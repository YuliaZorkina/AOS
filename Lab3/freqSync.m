clc; clear; close all;

Nfft = 2048;
M_data = 16;
M_pilot = 4;
SNR_dB = 0:1:15;
CP = 100;
NumOfdmSymbols = 10;

Naverage = 1000;
useMaxFrames = true;     % Флаг ограничения максимального числа кадров
useMinErrors = true;     % Флаг минимального количества ошибок
maxFrames = 100000;       % Максимальное количество кадров
minErrors = 1000;        % Минимальное суммарное количество ошибок
minFrames = 100;  

NumSC = 1600;
guard = (Nfft - NumSC)/2;
active_idx = (guard+1):(guard+NumSC);
pilot_step = 10;
pilot_idx_active = 1:pilot_step:NumSC;
pilot_idx = active_idx(pilot_idx_active);
data_idx_active = setdiff(1:NumSC, pilot_idx_active);
data_idx = active_idx(data_idx_active);
nPilots = length(pilot_idx);
nDataSymbols = length(data_idx);
kDataBits = log2(M_data);
kPilotBits = log2(M_pilot);
Ps = mean(abs(qammod(0:M_data-1, M_data, "UnitAveragePower", true)).^2) / log2(M_data) / Nfft;
rng(239)

%% ПУНКТ 1: BER при разных частотных расстройках, с/без эквалайзера

fn_list = [0.001 0.01 0.1 1];
maxFrames = 1000;
BER_noEQ = zeros(length(fn_list), length(SNR_dB));
BER_EQ = zeros(length(fn_list), length(SNR_dB));

for i_fn = 1:length(fn_list)
    fn = fn_list(i_fn);
    fprintf('Частотная расстройка fn = %.1f\n', fn);

    for i_snr = 1:length(SNR_dB)
        SNR = SNR_dB(i_snr);
        N0 = Ps / 10^(SNR/10);
        sigma = sqrt(N0/2);
        errorBits_noEQ = 0;
        errorBits_EQ = 0;
        totalBits = 0;
        frameCount = 0;

        while (errorBits_noEQ < Naverage ||errorBits_EQ < Naverage) && frameCount < maxFrames

            frameCount = frameCount + 1;

            tx_bits = randi([0 1], nDataSymbols*kDataBits, NumOfdmSymbols);
            tx_data = qammod(tx_bits(:), M_data,'InputType', 'bit', 'UnitAveragePower', true);
            tx_data = reshape(tx_data, nDataSymbols, NumOfdmSymbols);
            ofdm_symbols = zeros(Nfft, NumOfdmSymbols);
            ofdm_symbols(data_idx,:) = tx_data;

            pilot_bits = randi([0 1], nPilots*kPilotBits, NumOfdmSymbols);
            pilot_symbols = qammod(pilot_bits(:), M_pilot,'InputType', 'bit', 'UnitAveragePower', true);
            pilot_symbols = reshape(pilot_symbols, nPilots, NumOfdmSymbols);

            ofdm_symbols(pilot_idx,:) = pilot_symbols;
            ofdm_time = ifft(ifftshift(ofdm_symbols,1), Nfft, 1);
            cp = ofdm_time(end-CP+1:end,:);
            OFDM_vector = reshape([cp; ofdm_time], [], 1);

            n = (0:length(OFDM_vector)-1).';
            OFDM_vector_CFO = OFDM_vector .*  exp(1j*2*pi*fn*n/Nfft);
            noise = sigma*(randn(size(OFDM_vector_CFO)) + 1j*randn(size(OFDM_vector_CFO)));
  

            OFDM_vector_rx = OFDM_vector_CFO + noise;

            rx_matrix = reshape(OFDM_vector_rx,Nfft+CP, NumOfdmSymbols);
            rx_noCP = rx_matrix(CP+1:end,:);
            rx_symbols = fftshift(fft(rx_noCP, Nfft, 1),1);

            rx_data_noEQ = rx_symbols(data_idx,:);

            rx_bits_noEQ = qamdemod(rx_data_noEQ(:), M_data, 'OutputType', 'bit', 'UnitAveragePower', true);
            H_pilot = rx_symbols(pilot_idx,:) ./ pilot_symbols;
            H_est = interp1(pilot_idx(:), H_pilot, data_idx(:), 'linear', 'extrap');
            rx_data_EQ = rx_symbols(data_idx,:) ./ H_est;

            rx_bits_EQ = qamdemod(rx_data_EQ(:), M_data, 'OutputType', 'bit', 'UnitAveragePower', true);

            errorBits_noEQ = errorBits_noEQ + sum(tx_bits(:) ~= rx_bits_noEQ);
            errorBits_EQ = errorBits_EQ +  sum(tx_bits(:) ~= rx_bits_EQ);
            totalBits = totalBits + numel(tx_bits);
        end
        BER_noEQ(i_fn,i_snr) = errorBits_noEQ / totalBits;
        BER_EQ(i_fn,i_snr) = errorBits_EQ / totalBits;

        fprintf(['fn = %.3f | SNR = %2d dB | ''BER без EQ = %.6e | BER с EQ = %.6e | '  'Реализаций = %d\n'], fn, SNR, BER_noEQ(i_fn,i_snr),BER_EQ(i_fn,i_snr), frameCount);
    end
end

figure;
hold on;
grid on;

for i_fn = 1:length(fn_list)
    semilogy(SNR_dB, BER_noEQ(i_fn, :), '-o','LineWidth', 1.5, 'DisplayName', sprintf('f_n = %.3f', fn_list(i_fn)));
end
semilogy(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k','LineWidth',2, 'DisplayName','Теория AWGN');

xlabel('SNR, дБ');
ylabel('BER');
set(gca,'YScale','log');
title('BER при дробной частотной расстройке без эквалайзера');
legend('Location','southwest');
ylim([1e-5 1]);

figure;
hold on;
grid on;

for i_fn = 1:length(fn_list)
    semilogy(SNR_dB, BER_EQ(i_fn, :), '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('f_n = %.3f', fn_list(i_fn)));
end

semilogy(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k', 'LineWidth',2,'DisplayName','Теория AWGN');
xlabel('SNR, дБ');
ylabel('BER');
set(gca,'YScale','log');
title('BER при дробной частотной расстройке с эквалайзером');
legend('Location','southwest');
ylim([1e-5 1]);



%% ПУНКТ 2: Влияние количества OFDM-символов для когерентного накопления на BER при частотной синхронизации по CP

N_symb_list = 10;
N_acc_list = [1 2 5 10];
fn_list = [0.48 0.49 0.495];
BER_CFO = zeros(length(fn_list),length(N_acc_list),length(SNR_dB));

for i_fn = 1:length(fn_list)
    fn = fn_list(i_fn);
    for j = 1:length(N_symb_list)
        N_symb = N_symb_list(j);
        for j_acc = 1:length(N_acc_list)
            N_acc = N_acc_list(j_acc);
            for k = 1:length(SNR_dB)
                    errorBits = 0;
                    totalBits = 0;
                    frameCount = 0;

                while (~useMaxFrames || frameCount < maxFrames) && (~useMinErrors || errorBits < minErrors || frameCount < minFrames)
                    frameCount = frameCount+1;
                    data_bits = randi([0 1], nDataSymbols*kDataBits, N_symb);
                    mod_data = qammod(data_bits, M_data, 'InputType','bit','UnitAveragePower',true);
                    pilot_bits = randi([0 1], nPilots*kPilotBits, N_symb);
                    mod_pilots = qammod(pilot_bits, M_pilot, 'InputType','bit','UnitAveragePower',true);
                    
                    ofdm_symbols = zeros(Nfft, N_symb);
                    ofdm_symbols(pilot_idx,:) = mod_pilots;
                    ofdm_symbols(data_idx,:) = mod_data;
                    ofdm_time = ifft(ifftshift(ofdm_symbols,1), Nfft, 1);
                    cp = ofdm_time(end-CP+1:end,:);
                    tx_vector = reshape([cp; ofdm_time], [], 1);
                    % Частотная расстройка
                    n = (0:length(tx_vector)-1).';
                    tx_CFO = tx_vector.*exp(1j*2*pi*fn*n/Nfft);

                    % Добавление шума
                    N0 = Ps/10^(SNR_dB(k)/10);
                    sigma = sqrt(N0/2);
                    noise = sigma*randn(size(tx_CFO))+1j*sigma*randn(size(tx_CFO));
                    rx_signal = tx_CFO+noise;

                    % Оценка частотной расстройки по ЦП
                    rx_matrix = reshape(rx_signal,Nfft+CP,N_symb);
                    rx_cp = rx_matrix(1:CP,1:N_acc);
                    rx_tail = rx_matrix(Nfft+1:Nfft+CP,1:N_acc);

                    R = sum(conj(rx_cp(:)) .* rx_tail(:));
                    fn_est = angle(R)/(2*pi);

                    % Компенсация частотной расстройки
                    n_rx = (0:length(rx_signal)-1).';
                    rx_corrected = rx_signal.*exp(-1j*2*pi*fn_est*n_rx/Nfft);

                    % Приём OFDM-сигнала
                    rx_matrix_corrected = reshape(rx_corrected,Nfft+CP,N_symb);
                    rx_noCP = rx_matrix_corrected(CP+1:end,:);
                    rx_freq = fftshift(fft(rx_noCP,Nfft,1),1);
                    rx_data = rx_freq(data_idx,:);

                    % Эквалайзер по пилотам
                    H_pilot = rx_freq(pilot_idx,:)./mod_pilots;
                    H_est = interp1(pilot_idx(:),H_pilot,data_idx(:),'linear','extrap');
                    rx_data_eq = rx_data./H_est;

                    % Демодуляция и подсчёт BER
                    rx_bits = qamdemod(rx_data_eq(:),M_data,'OutputType','bit','UnitAveragePower',true);
                    errorBits = errorBits+sum(data_bits(:)~=rx_bits);
                    totalBits = totalBits+numel(data_bits);
                end
                BER_CFO(i_fn,j_acc,k) = errorBits/totalBits;
                fprintf('fn = %.2f, N_acc = %d, SNR = %d dB, BER = %.6g\n',fn,N_acc,SNR_dB(k),BER_CFO(i_fn,j_acc,k));
            end
        end
    end
end

styles = {'-','--',':','-.','-','--'};
for i_fn = 1:length(fn_list)
    figure;
    hold on;
    grid on;
    for j_acc = 1:length(N_acc_list) 
        BER_plot = squeeze(BER_CFO(i_fn,j_acc,:));
        semilogy(SNR_dB,BER_plot, 'LineStyle',styles{j_acc},'LineWidth',1.5, 'DisplayName',sprintf('N_{acc} = %d',N_acc_list(j_acc)));
    end
    semilogy(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k','LineWidth',2, 'DisplayName','Теория AWGN');
    set(gca,'YScale','log');
    xlabel('SNR, дБ');
    ylabel('BER');
    title(sprintf('BER при частотной синхронизации по ЦП, f_n = %.3f',fn_list(i_fn)));
    legend('Location','best');
    ylim([1e-5 1]);
    hold off;
end



%% Пунтк 3  ВЛИЯНИЯ ДЛИНЫ ЦИКЛИЧЕСКОГО ПРЕФИКСА 
CP_list = [4 32 512];
fn_list = [0.48 0.49];
N_acc = NumOfdmSymbols;

BER_CP = zeros(length(fn_list),length(CP_list),length(SNR_dB));

for i_fn = 1:length(fn_list)
    fn = fn_list(i_fn);
    for j_cp = 1:length(CP_list)
        CP_current = CP_list(j_cp);
        for k = 1:length(SNR_dB)

             errorBits = 0;
             totalBits = 0;
             frameCount = 0;

             while (~useMaxFrames || frameCount < maxFrames) && (~useMinErrors || errorBits < minErrors || frameCount < minFrames)
                    frameCount = frameCount+1;

                % Генерация информационных битов
                data_bits = randi([0 1],nDataSymbols*kDataBits,NumOfdmSymbols);
                mod_data = qammod(data_bits(:),M_data,'InputType','bit','UnitAveragePower',true);
                mod_data = reshape(mod_data,nDataSymbols,NumOfdmSymbols);

                % Генерация пилотных символов
                pilot_bits = randi([0 1],nPilots*kPilotBits,NumOfdmSymbols);
                mod_pilots = qammod(pilot_bits(:),M_pilot,'InputType','bit','UnitAveragePower',true);
                mod_pilots = reshape(mod_pilots,nPilots,NumOfdmSymbols);

                % Формирование частотного OFDM-сигнала
                ofdm_symbols = zeros(Nfft,NumOfdmSymbols);
                ofdm_symbols(data_idx,:) = mod_data;
                ofdm_symbols(pilot_idx,:) = mod_pilots;

                % OFDM-модуляция
                ofdm_time = ifft(ifftshift(ofdm_symbols,1),Nfft,1);
                cp = ofdm_time(end-CP_current+1:end,:);
                tx_vector = reshape([cp;ofdm_time],[],1);

                % Внесение частотной расстройки
                n = (0:length(tx_vector)-1).';
                tx_CFO = tx_vector.*exp(1j*2*pi*fn*n/Nfft);

                % Добавление АБГШ
                N0 = Ps/10^(SNR_dB(k)/10);
                sigma = sqrt(N0/2);
                noise = sigma*randn(size(tx_CFO))+1j*sigma*randn(size(tx_CFO));
                rx_signal = tx_CFO+noise;

                % Оценка частотной расстройки по ЦП
                rx_matrix = reshape(rx_signal,Nfft+CP_current,NumOfdmSymbols);
                rx_cp = rx_matrix(1:CP_current,1:N_acc);
                rx_tail = rx_matrix(Nfft+1:Nfft+CP_current,1:N_acc);

                R = sum(conj(rx_cp(:)).*rx_tail(:));
                fn_est = angle(R)/(2*pi);

                % Компенсация частотной расстройки
                n_rx = (0:length(rx_signal)-1).';
                rx_corrected = rx_signal.*exp(-1j*2*pi*fn_est*n_rx/Nfft);

                % Удаление ЦП и OFDM-демодуляция
                rx_matrix_corrected = reshape(rx_corrected,Nfft+CP_current,NumOfdmSymbols);
                rx_noCP = rx_matrix_corrected(CP_current+1:end,:);
                rx_freq = fftshift(fft(rx_noCP,Nfft,1),1);

                % Оценка канала по пилотам и эквализация
                rx_data = rx_freq(data_idx,:);
                H_pilot = rx_freq(pilot_idx,:)./mod_pilots;
                H_est = interp1(pilot_idx(:),H_pilot,data_idx(:),'linear','extrap');
                rx_data_eq = rx_data./H_est;

                % Демодуляция и подсчёт ошибок
                rx_bits = qamdemod(rx_data_eq(:),M_data,'OutputType','bit','UnitAveragePower',true);
                errorBits = errorBits+sum(data_bits(:)~=rx_bits);
                totalBits = totalBits+numel(data_bits);
            end

            BER_CP(i_fn,j_cp,k) = errorBits/totalBits;
            fprintf('fn = %.3f, CP = %d, SNR = %d dB, BER = %.6g\n',fn,CP_current,SNR_dB(k),BER_CP(i_fn,j_cp,k));
        end
    end
end

colors = lines(length(CP_list));
for i_fn = 1:length(fn_list)
    figure;
    hold on;
    grid on;
    for j_cp = 1:length(CP_list)
        BER_plot = squeeze(BER_CP(i_fn,j_cp,:));
        semilogy(SNR_dB,BER_plot,'Color',colors(j_cp,:),'LineWidth',1.5,'DisplayName',sprintf('CP = %d',CP_list(j_cp)));
    end
    semilogy(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k','LineWidth',2, 'DisplayName','Теория AWGN');
    set(gca,'YScale','log');
    xlabel('SNR, дБ');
    ylabel('BER');
    title(sprintf('BER при частотной синхронизации по ЦП, f_n = %.3f',fn_list(i_fn)));
    legend('Location','best');
    ylim([1e-5 1]);
    hold off;
end



%% ЧАСТОТНАЯ СИНХРОНИЗАЦИЯ МЕТОДОМ МИННА 

fn_sets = {-0.5:0.25:0.5,-1.5:0.5:1.5};
range_names = {'(-0.5, 0.5)','(-2, 2)'};

useMaxFrames = true;
useMinErrors = true;
maxFrames = 10000;
minErrors = Naverage;
minFrames = 100;

L = Nfft/4;
pn = comm.PNSequence('Polynomial',[9 5 0], ...
    'InitialConditions',ones(1,9), ...
    'SamplesPerFrame',511);
mseq_bits = pn();
release(pn);

% Преобразование М-последовательности в биполярную
mseq = 2*double(mseq_bits)-1;

% Дополнение до длины 512
mseq = [mseq;mseq(1)];
% Нормировка мощности преамбулы
A = mseq;
MinnBody = [A;A;-A;-A];

MinnPreamble = [MinnBody(end-CP+1:end);MinnBody];

BER_Minn = cell(1,length(fn_sets));

for i_range = 1:length(fn_sets)

    fn_list = fn_sets{i_range};
    BER_Minn{i_range} = zeros(length(fn_list),length(SNR_dB));

    for i_fn = 1:length(fn_list)

        fn = fn_list(i_fn);

        for k = 1:length(SNR_dB)

            errorBits = 0;
            totalBits = 0;
            frameCount = 0;

            while (~useMaxFrames || frameCount < maxFrames) && (~useMinErrors || errorBits < minErrors || frameCount < minFrames)

                frameCount = frameCount+1;

                data_bits = randi([0 1],nDataSymbols*kDataBits,NumOfdmSymbols);
                mod_data = qammod(data_bits(:),M_data,'InputType','bit','UnitAveragePower',true);
                mod_data = reshape(mod_data,nDataSymbols,NumOfdmSymbols);
                pilot_bits = randi([0 1],nPilots*kPilotBits,NumOfdmSymbols);
                mod_pilots = qammod(pilot_bits(:),M_pilot,'InputType','bit','UnitAveragePower',true);
                mod_pilots = reshape(mod_pilots,nPilots,NumOfdmSymbols);
                ofdm_symbols = zeros(Nfft,NumOfdmSymbols);
                ofdm_symbols(data_idx,:) = mod_data;
                ofdm_symbols(pilot_idx,:) = mod_pilots;

                ofdm_time = ifft(ifftshift(ofdm_symbols,1),Nfft,1);

                % Добавление циклического префикса
                cp_data = ofdm_time(end-CP+1:end,:);
                tx_data = reshape([cp_data;ofdm_time],[],1);

                tx_frame = [MinnPreamble;tx_data];

                n = (0:length(tx_frame)-1).';
                tx_CFO = tx_frame.*exp(1j*2*pi*fn*n/Nfft);

                N0 = Ps/10^(SNR_dB(k)/10);
                sigma = sqrt(N0/2);
                noise = sigma*randn(size(tx_CFO))+1j*sigma*randn(size(tx_CFO));
                rx_signal = tx_CFO+noise;

                % Извлечение принятой преамбулы
                rx_preamble = rx_signal(1:length(MinnPreamble));

                rx_minn = rx_preamble(CP+1:end);

                rx_A1 = rx_minn(1:L);
                rx_A2 = rx_minn(L+1:2*L);
                rx_A3 = rx_minn(2*L+1:3*L);
                rx_A4 = rx_minn(3*L+1:4*L);

                P1 = sum(conj(rx_A1).*rx_A2);
                P2 = sum(conj(rx_A3).*rx_A4);
                P = P1+P2;
                fn_est = (2/pi) * angle(P);

                n_rx = (0:length(rx_signal)-1).';
                rx_corrected = rx_signal.*exp(-1j*2*pi*fn_est*n_rx/Nfft);
                rx_data_frame = rx_corrected(length(MinnPreamble)+1:end);

                rx_matrix = reshape(rx_data_frame,Nfft+CP,NumOfdmSymbols);
                rx_noCP = rx_matrix(CP+1:end,:);
                rx_freq = fftshift(fft(rx_noCP,Nfft,1),1);

                rx_data = rx_freq(data_idx,:);
                H_pilot = rx_freq(pilot_idx,:)./mod_pilots;
                H_est = interp1(pilot_idx(:),H_pilot,data_idx(:),'linear','extrap');
                rx_data_eq = rx_data./H_est;

                rx_bits = qamdemod(rx_data_eq(:),M_data,'OutputType','bit','UnitAveragePower',true);

                errorBits = errorBits+sum(data_bits(:)~=rx_bits);
                totalBits = totalBits+numel(data_bits);

            end
            BER_Minn{i_range}(i_fn,k) = errorBits/totalBits;
            fprintf('Диапазон %s | fn = %.2f | SNR = %d dB | Кадров = %d | Ошибок = %d | Битов = %d | BER = %.6g\n', ...
                range_names{i_range},fn,SNR_dB(k),frameCount,errorBits,totalBits,BER_Minn{i_range}(i_fn,k));

        end
    end
end
for i_range = 1:length(fn_sets)
    fn_list = fn_sets{i_range};
    figure;
    hold on;
    grid on;
    colors = lines(length(fn_list));
    for i_fn = 1:length(fn_list)

        BER_plot = BER_Minn{i_range}(i_fn,:);

        semilogy(SNR_dB,BER_plot, 'Color',colors(i_fn,:),'LineStyle','-', 'LineWidth',1.5, 'DisplayName',sprintf('f_n = %.2f',fn_list(i_fn)));
    end
    semilogy(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k','LineWidth',2, 'DisplayName','Теория AWGN');
    set(gca,'YScale','log');
    xlabel('SNR, дБ');
    ylabel('BER');
    title(sprintf('Частотная синхронизация методом Минна, f_n \\in %s',range_names{i_range}));
    legend('Location','best');
    ylim([1e-5 1]);
    hold off;

end
