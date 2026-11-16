clc; clear; close all;

Nfft = 2048;
M_data = 16;
M_pilot = 4;
SNR_dB = 0:1:15;
CP = 100;
NumOfdmSymbols = 10;
Naverage = 100;

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

Ps = mean(abs(qammod(0:M_data-1, M_data, "UnitAveragePower", true)).^2)/log2(M_data)/Nfft;

%% ПУНКТ 1: BER при разных временных сдвигах, с/без эквалайзера

tau_cases = [-150, -50, 50, 0];
tau_labels = {'δ < -NCP', '-NCP ≤ δ < 0', 'δ > 0', 'δ = 0'};
eq_cases = [0, 1];
eq_labels = {'без экв.', 'с экв.'};

figure(1); clf; hold on; grid on;
colors = {'r','g','b','m','c','k','y','b'};
markers = {'o','s','d','^','v','>','<','p'};

for i_tau = 1:length(tau_cases)
    tau = tau_cases(i_tau);
    for i_eq = 1:length(eq_cases)
        use_eq = eq_cases(i_eq);
        BER = zeros(1, length(SNR_dB));
        
        for k = 1:length(SNR_dB)
            errorBits = 0; totalBits = 0;
            while errorBits < Naverage
                data_bits = randi([0 1], nDataSymbols*kDataBits, NumOfdmSymbols);
                mod_data = qammod(data_bits, M_data, 'InputType','bit','UnitAveragePower',true);
                pilot_bits = randi([0 1], nPilots*kPilotBits, NumOfdmSymbols);
                mod_pilots = qammod(pilot_bits, M_pilot, 'InputType','bit','UnitAveragePower',true);
                
                ofdm_symbols = zeros(Nfft, NumOfdmSymbols);
                ofdm_symbols(pilot_idx,:) = mod_pilots;
                ofdm_symbols(data_idx,:) = mod_data;
                ofdm_time = ifft(ifftshift(ofdm_symbols,1), Nfft, 1);
                
                cp = ofdm_time(end-CP+1:end,:);
                tx_vector = reshape([cp; ofdm_time], [], 1);
                
                % Временной сдвиг (циклический)
                shifted_signal = circshift(tx_vector, tau);
                
                % AWGN
                N0 = Ps / 10^(SNR_dB(k)/10);
                sigma = sqrt(N0/2);
                noise = sigma*(randn(size(shifted_signal)) + 1j*randn(size(shifted_signal)));
                rx_signal = shifted_signal + noise;
                
                % Приём БЕЗ синхронизации
                rx_matrix = reshape(rx_signal, Nfft+CP, NumOfdmSymbols);
                rx_no_cp = rx_matrix(CP+1:end,:);
                rx_freq = fftshift(fft(rx_no_cp, Nfft, 1), 1);
                
                rx_data = rx_freq(data_idx,:);
                if use_eq
                    H_est = mean(rx_freq(pilot_idx,:)./mod_pilots, 2);
                    H_interp = interp1(pilot_idx, H_est, data_idx, 'linear','extrap');
                    rx_data_eq = rx_data ./ H_interp.';
                else
                    rx_data_eq = rx_data;
                end
                
                rx_bits = qamdemod(rx_data_eq, M_data, 'OutputType','bit','UnitAveragePower',true);
                errorBits = errorBits + biterr(rx_bits, data_bits);
                totalBits = totalBits + numel(data_bits);
            end
            BER(k) = errorBits/totalBits;
        end
        idx_line = (i_tau-1)*2 + i_eq;
        semilogy(SNR_dB, BER, '-o', 'Color', colors{idx_line}, 'Marker', markers{idx_line}, ...
            'LineWidth',1.5, 'MarkerSize',5, ...
            'DisplayName', sprintf('τ=%d, %s', tau, eq_labels{i_eq}));
        fprintf('τ=%d, %s: BER@10дБ=%.4e\n', tau, eq_labels{i_eq}, BER(11));
    end
end

semilogy(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k', 'LineWidth',2, 'DisplayName','Теория AWGN');
xlabel('SNR (дБ)'); ylabel('BER');
title('Пункт 1: BER при разных временных сдвигах, с/без эквалайзера (CP=100)');
legend('Location','southwest'); 
set(gca, 'YScale', 'log');  
ylim([1e-5 1]);
hold off;

%% ПУНКТ 2: АКФ при разных N_symb (фиксированный SNR) с усред

SNR_fixed = 10;
N_symb_list = [3, 10, 30, 100];
tau_true = 50;
N_realizations = 50;  % Число реализаций для усреднения АКФ

figure(2); clf; hold on; grid on;

for j = 1:length(N_symb_list)
    N_symb = N_symb_list(j);
    
    AllDelay = 0:(CP-1);
    ACF_avg = zeros(length(AllDelay), 1);  
    
    for r = 1:N_realizations
        % генерация данных (новая для каждой реализации)
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
        
        % временной сдвиг
        shifted_signal = circshift(tx_vector, tau_true);
        N0 = Ps / 10^(SNR_fixed/10);
        sigma = sqrt(N0/2);
        noise = sigma*(randn(size(shifted_signal)) + 1j*randn(size(shifted_signal)));
        rx_signal = shifted_signal + noise;
        
        % АКФ по ЦП для данной реализации
        ACF = zeros(length(AllDelay), 1);
        
        for j_delay = AllDelay
            for j_symb = 0:(N_symb-2)  % отбрасываем последний символ
                sym_start = j_symb * (CP + Nfft) + 1;
                left_idx = sym_start + j_delay + (0:CP-1);
                right_idx = left_idx + Nfft;
                
                if max(right_idx) <= length(rx_signal)
                    left = rx_signal(left_idx);
                    right = rx_signal(right_idx);
                    ACF(j_delay+1) = ACF(j_delay+1) + sum(left .* conj(right));
                end
            end
        end
        
        ACF = ACF / (N_symb - 1);
        ACF_avg = ACF_avg + ACF;  % Накопление
    end
    
    ACF_avg = ACF_avg / N_realizations;  % Усреднение
    
    plot(AllDelay, abs(ACF_avg), 'LineWidth', 1.3, 'DisplayName', sprintf('N_symb = %d', N_symb));
    
    [~, peak_idx] = max(abs(ACF_avg));
    fprintf('N_symb=%d: пик АКФ на задержке %d (истинная τ=%d)\n', ...
        N_symb, AllDelay(peak_idx), tau_true);
end

xlabel('Отсчёт задержки'); ylabel('|АКФ|');
title(sprintf('Пункт 2: АКФ при разных N_symb (SNR=%d дБ, τ=%d, усреднено по %d реализациям)', SNR_fixed, tau_true, N_realizations));
legend('Location','best'); hold off;

%% ПУНКТ 3: BER при разных N_symb (синхронизация по ЦП)

N_symb_list = [5, 10, 20, 40];
tau_true = 50;
eq_cases = [0,1]; 
eq_labels = { 'без экв.','с экв.'};

figure(3); clf; hold on; grid on;

for j = 1:length(N_symb_list)
    N_symb = N_symb_list(j);
    for i_eq = 1:length(eq_cases)
        use_eq = eq_cases(i_eq);
        BER = zeros(1, length(SNR_dB));
        
        for k = 1:length(SNR_dB)
            errorBits = 0; totalBits = 0;
            while errorBits < Naverage
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
                
                shifted_signal = circshift(tx_vector, tau_true);
                N0 = Ps / 10^(SNR_dB(k)/10);
                sigma = sqrt(N0/2);
                noise = sigma*(randn(size(shifted_signal)) + 1j*randn(size(shifted_signal)));
                rx_signal = shifted_signal + noise;
                
                % Синхронизация по ЦП
                AllDelay = 0:(CP-1);
                ACF = zeros(length(AllDelay), 1);
                
                % Используем только N_symb-1 символов 
                for j_delay = AllDelay
                    for j_symb = 0:(N_symb-2)  % отбрасываем последний символ
                        sym_start = j_symb * (CP + Nfft) + 1;
                        left_idx = sym_start + j_delay + (0:CP-1);
                        right_idx = left_idx + Nfft;
                        
                        left = rx_signal(left_idx);
                        right = rx_signal(right_idx);
                        ACF(j_delay+1) = ACF(j_delay+1) + sum(left .* conj(right));
                    end
                end
                ACF = ACF / (N_symb - 1);
                
                [~, Delay] = max(abs(ACF));
                
                % приём (используем N_symb-1 символов)
                Signal_sync = rx_signal(Delay:end);
                N_use = N_symb - 1;
                Signal_sync = Signal_sync(1:(Nfft+CP)*N_use);
                symAWGN = reshape(Signal_sync, Nfft+CP, N_use);
                ncpAWGN = symAWGN(CP+1:end,:);
                rx_freq = fftshift(fft(ncpAWGN, Nfft, 1), 1);
                
                rx_data = rx_freq(data_idx,:);
                if use_eq
                    H_est = mean(rx_freq(pilot_idx,:)./mod_pilots(:,1:N_use), 2);
                    H_interp = interp1(pilot_idx, H_est, data_idx, 'linear','extrap');
                    rx_data_eq = rx_data ./ H_interp.';
                else
                    rx_data_eq = rx_data;
                end
                
                rx_bits = qamdemod(rx_data_eq, M_data, 'OutputType','bit','UnitAveragePower',true);
                errorBits = errorBits + biterr(rx_bits, data_bits(:,1:N_use));
                totalBits = totalBits + numel(data_bits(:,1:N_use));
            end
            BER(k) = errorBits/totalBits;
        end
        
        semilogy(SNR_dB, BER, '-o', 'LineWidth',1.5, 'MarkerSize',5,'DisplayName', sprintf('N_symb=%d, %s', N_symb, eq_labels{i_eq}));
        fprintf('N_symb=%d, %s: BER@10дБ=%.4e\n', N_symb, eq_labels{i_eq}, BER(11));
    end
end

semilogy(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k', 'LineWidth',2, 'DisplayName','Теория AWGN');
xlabel('SNR (дБ)'); ylabel('BER');
title(sprintf('Пункт 3: BER при разных N_symb (синхр. по ЦП, τ=%d)', tau_true));
legend('Location','southwest'); 
set(gca, 'YScale', 'log'); 
ylim([1e-5 1]);
hold off;

%% ПУНКТ 4: МЕТРИКА МИННА с усред

L = Nfft/4; % длина одного сегмента
tau_true = 50; % истинная временная задержка
MaxDelay = 200; % максимальная исследуемая задержка
SNR_fixed = 17; % SNR для 4 пункта
N_realizations_minn = 50;  % Число реализаций для усреднения метрики

pn = comm.PNSequence( ...
    'Polynomial', [9 5 0], ...
    'InitialConditions', ones(1,9), ...
    'SamplesPerFrame', L);

A_bits = pn();
A = 2*A_bits - 1;            
sync_symbol = [A; A; -A; -A];

AllDelay = 0:MaxDelay;
Metric_avg = zeros(length(AllDelay), 1);  

for r = 1:N_realizations_minn
    % Генерация информационных данных (новая для каждой реализации)
    data_bits = randi([0 1], nDataSymbols*kDataBits, N_symb);
    mod_data = qammod(data_bits, M_data,'InputType','bit','UnitAveragePower',true);
    pilot_bits = randi([0 1], nPilots*kPilotBits, N_symb);
    mod_pilots = qammod(pilot_bits, M_pilot,'InputType','bit','UnitAveragePower',true);

    % Формирование OFDM-символов
    ofdm_symbols = zeros(Nfft, N_symb);
    ofdm_symbols(pilot_idx,:) = mod_pilots;
    ofdm_symbols(data_idx,:) = mod_data;
    ofdm_time = ifft(ifftshift(ofdm_symbols,1), Nfft, 1);

    % Добавление циклического префикса
    cp = ofdm_time(end-CP+1:end,:);
    tx_payload = reshape([cp; ofdm_time], [], 1);

    % Добавление преамбулы перед информационными символами
    tx_signal = [sync_symbol; tx_payload];

    % Внесение временной задержки
    shifted_signal = [zeros(tau_true,1); tx_signal];

    % Добавление AWGN
    N0 = Ps / 10^(SNR_fixed/10);
    sigma = sqrt(N0/2);
    noise = sigma*(randn(size(shifted_signal)) + 1j*randn(size(shifted_signal)));
    rx_signal = shifted_signal + noise;

    % Вычисление метрики Минна для данной реализации
    Metric = zeros(length(AllDelay), 1);

    for j = 1:length(AllDelay)
        d = AllDelay(j);

        % Первая половина синхросимвола [A A]
        left1 = rx_signal(d + (1:L)).';
        right1 = rx_signal(d + L + (1:L)).';

        % Вторая половина синхросимвола [-A -A]
        left2 = rx_signal(d + 2*L + (1:L)).';
        right2 = rx_signal(d + 3*L + (1:L)).';

        % Корреляция между сегментами первой половины
        P1 = conv(left1, fliplr(conj(right1)), 'valid');

        % Корреляция между сегментами второй половины
        P2 = conv(left2, fliplr(conj(right2)), 'valid');

        % Суммарная корреляция
        P = P1 + P2;

        % Энергии
        E1 = sum(abs(left1).^2) + sum(abs(right1).^2);
        E2 = sum(abs(left2).^2) + sum(abs(right2).^2);

        % Нормированная метрика Минна
        Metric(j) = abs(P)^2 / (E1*E2);
    end
    
    Metric_avg = Metric_avg + Metric;  % Накопление
end

Metric_avg = Metric_avg / N_realizations_minn;  % Усреднение

% Определение максимума усреднённой метрики
[~, peak_idx] = max(Metric_avg);
Delay = AllDelay(peak_idx);

figure(4); clf;
plot(AllDelay, Metric_avg, 'b', 'LineWidth', 1.5);
hold on;
xline(tau_true, '--r', 'Истинная задержка', 'LineWidth', 1.5);
xline(Delay, '--g', 'Оцененная задержка', 'LineWidth', 1.5);
grid on;
xlabel('Временная задержка d, отсчёты');
ylabel('M(d)');
title(sprintf('Метрика Минна (SNR = %d дБ, усреднено по %d реализациям)', SNR_fixed, N_realizations_minn));
legend('Метрика Минна', 'Истинная задержка', 'Оцененная задержка', 'Location', 'best');
hold off;

eq_cases = [0 1];
eq_labels = {'Без эквалайзера', 'С эквалайзером'};
figure(5); clf;
hold on;
grid on;
for i_eq = 1:length(eq_cases)

    use_eq = eq_cases(i_eq);
    BER = zeros(1,length(SNR_dB));
    for k = 1:length(SNR_dB)
        errorBits = 0;
        totalBits = 0;
        while errorBits < Naverage

            data_bits = randi([0 1], nDataSymbols*kDataBits, NumOfdmSymbols);
            mod_data = qammod(data_bits, M_data, 'InputType','bit','UnitAveragePower',true);
            pilot_bits = randi([0 1], nPilots*kPilotBits, NumOfdmSymbols);
            mod_pilots = qammod(pilot_bits, M_pilot, 'InputType','bit','UnitAveragePower',true);

            ofdm_symbols = zeros(Nfft, NumOfdmSymbols);
            ofdm_symbols(pilot_idx,:) = mod_pilots;
            ofdm_symbols(data_idx,:) = mod_data;
            ofdm_time = ifft(ifftshift(ofdm_symbols,1), Nfft, 1);
            cp = ofdm_time(end-CP+1:end,:);

            tx_payload = reshape([cp; ofdm_time], [], 1);

            tx_signal = [sync_symbol; tx_payload];

            shifted_signal = [zeros(tau_true,1); tx_signal];

            N0 = Ps / 10^(SNR_dB(k)/10);
            sigma = sqrt(N0/2);
            noise = sigma*(randn(size(shifted_signal)) + 1j*randn(size(shifted_signal)));
            rx_signal = shifted_signal + noise;

            Metric = zeros(length(AllDelay),1);

            for j = 1:length(AllDelay)

                d = AllDelay(j);

                left1 = rx_signal(d + (1:L)).';
                right1 = rx_signal(d + L + (1:L)).';

                left2 = rx_signal(d + 2*L + (1:L)).';
                right2 = rx_signal(d + 3*L + (1:L)).';

                % Корреляции
                P1 = conv(left1, fliplr(conj(right1)), 'valid');
                P2 = conv(left2, fliplr(conj(right2)), 'valid');
                P = P1 + P2;

                % Энергии
                E1 = sum(abs(left1).^2) + sum(abs(right1).^2);
                E2 = sum(abs(left2).^2) + sum(abs(right2).^2);
                % Метрика Минна
                Metric(j) = abs(P)^2 / (E1*E2);

            end
            % Оценка временной задержки
            [~, peak_idx] = max(Metric);
            Delay = AllDelay(peak_idx);
            % Начало полезной OFDM-последовательности
            payload_start = Delay + Nfft + 1;
            N_total = (Nfft+CP)*NumOfdmSymbols;
            % Проверка границ массива
            if payload_start + N_total - 1 > length(rx_signal)
                continue;  
            end

            % Извлечение полезного сигнала
            Signal_sync = rx_signal(payload_start : payload_start+N_total-1);

            rx_matrix = reshape(Signal_sync, Nfft+CP, NumOfdmSymbols);

            rx_no_cp = rx_matrix(CP+1:end,:);

            rx_freq = fftshift(fft(rx_no_cp,Nfft,1),1);

            rx_data = rx_freq(data_idx,:);

            if use_eq
                H_est = mean(rx_freq(pilot_idx,:) ./ mod_pilots, 2);
                H_interp = interp1(pilot_idx, H_est, data_idx, 'linear','extrap');
                rx_data_eq = rx_data ./ H_interp.';
            else
                rx_data_eq = rx_data;
            end

            rx_bits = qamdemod(rx_data_eq, M_data, 'OutputType','bit', 'UnitAveragePower',true);
            errorBits = errorBits + biterr(rx_bits, data_bits);
            totalBits = totalBits + numel(data_bits);
        end
        BER(k) = errorBits / totalBits;
    end
    semilogy(SNR_dB, BER, '-o', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName',eq_labels{i_eq});
    fprintf('%s: BER@10дБ = %.4e\n', eq_labels{i_eq}, BER(11));
end

semilogy(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k', 'LineWidth',2,'DisplayName','Теория AWGN');
xlabel('SNR (дБ)');
ylabel('BER');

title(sprintf('Помехоустойчивость OFDM при временной синхронизации по алгоритму Минна (CP=%d)',CP));
legend('Location','southwest');
set(gca,'YScale','log');
ylim([1e-5 1]);
grid on;
hold off;
