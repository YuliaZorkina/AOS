clc; clear; close all;

Nfft = 2048;                       % Размер FFT
M_data = 16;                       % Модуляция данных: 16-QAM
M_pilot = 4;                       % Модуляция пилотов: QPSK
CP_len = Nfft/4;                   % Длина CP = 512 отсчётов
T_sym = 100e-6;                    % Длительность OFDM-символа, 100 мкс
Fs = Nfft/T_sym;                   % Частота дискретизации, 20.48 МГц
SNR_dB = 1:2:30;                   % Диапазон ОСШ, дБ
NumOfdmSymbols = 200;              % Количество OFDM-символов
NumSC_fad = 1200;

Es = mean(abs(qammod(0:M_data-1, M_data)).^2) / Nfft;
%% ПРОФИЛИ КАНАЛОВ EPA, EVA, ETU

% Задержки лучей, нс
EPA_delays_ns = [0, 30, 70, 80, 110, 190, 410];
EVA_delays_ns = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510];
ETU_delays_ns = [0, 50, 120, 200, 230, 500, 1600, 2300, 5000];

% Средние относительные мощности лучей, дБ
EPA_gains_dB = [0, -1, -2, -3, -8, -17.2, -20.8];
EVA_gains_dB = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7, -12, -16.9];
ETU_gains_dB = [-1, -1, -1, 0, 0, 0, -3, -5, -7];

% Перевод задержек в секунды
EPA_delays = EPA_delays_ns*1e-9;
EVA_delays = EVA_delays_ns*1e-9;
ETU_delays = ETU_delays_ns*1e-9;

% Перевод мощностей из дБ в линейный масштаб
EPA_gains_lin = 10.^(EPA_gains_dB/10);
EVA_gains_lin = 10.^(EVA_gains_dB/10);
ETU_gains_lin = 10.^(ETU_gains_dB/10);

% Нормировка суммарной мощности каждого профиля к единице
EPA_gains_lin = EPA_gains_lin/sum(EPA_gains_lin);
EVA_gains_lin = EVA_gains_lin/sum(EVA_gains_lin);
ETU_gains_lin = ETU_gains_lin/sum(ETU_gains_lin);

channel_names = {'EPA', 'EVA', 'ETU'};
colors_ch = {'b', 'r', 'g'};

%% ДОПОЛНИТЕЛЬНОЕ ЗАДАНИЕ 1: BER БЕЗ КАНАЛА

fprintf('\nДополнительное задание 1: BER без канала\n');

NumSC_list_dop = [200 800 1200 1600];
BER_no_channel = zeros(length(NumSC_list_dop),1);

for idxSC = 1:length(NumSC_list_dop)

    NumSC_dop = NumSC_list_dop(idxSC);
    NumGuardSC_dop = round((Nfft-NumSC_dop)/2);

    N_bits = NumSC_dop*log2(M_data);
    totalBits = NumOfdmSymbols*N_bits;
    errorBits = 0;

    for sym = 1:NumOfdmSymbols

        % Генерация битов и модуляция
        data_bits = randi([0 1], N_bits, 1);

        mod_data = qammod(data_bits, M_data,'InputType', 'bit', 'UnitAveragePower', true);

        % Формирование OFDM-символа
        Sp = zeros(Nfft,1);
        Sp(NumGuardSC_dop+(1:NumSC_dop)) = mod_data;

        Sig = ifft(ifftshift(Sp));
        Sig_CP = [Sig(end-CP_len+1:end); Sig];

        % Идеальный канал: без шума и искажений
        Sig_rec = Sig_CP(CP_len+1:CP_len+Nfft);
        Sp_rec = fftshift(fft(Sig_rec));

        % Демодуляция в биты
        demod_bits = qamdemod(Sp_rec(NumGuardSC_dop+(1:NumSC_dop)),M_data, 'OutputType', 'bit', 'UnitAveragePower', true);

        errorBits = errorBits + sum(demod_bits(:) ~= data_bits(:));
    end

    BER_no_channel(idxSC) = errorBits/totalBits;

    fprintf('NumSC = %4d: BER = %e\n', NumSC_dop, BER_no_channel(idxSC));
end

figure(10);
clf;
bar(NumSC_list_dop, BER_no_channel, 0.6);
set(gca,'YScale','log');
xlabel('Количество поднесущих (NumSC)');
ylabel('BER');
title('BER системы без канала (идеальный случай)');
grid on;
ylim([1e-10 1]);


%% ДОПОЛНИТЕЛЬНОЕ ЗАДАНИЕ 2: СПЕКТР OFDM-СИГНАЛА МЕТОДОМ УЭЛЧА

fprintf('\nДополнительное задание 2: Спектр OFDM-сигнала\n');

NumSC_spectr = 1200;
NumGuardSC_spectr = round((Nfft-NumSC_spectr)/2);
N_symbols_psd = 100;

sig_len_noCP = Nfft*N_symbols_psd;
sig_len_withCP = (Nfft+CP_len)*N_symbols_psd;

Sig_noCP_long = zeros(sig_len_noCP,1);
Sig_withCP_long = zeros(sig_len_withCP,1);

for sym = 1:N_symbols_psd

    % Генерация OFDM-данных
    data = randi([0 M_data-1], NumSC_spectr, 1);

    mod_data = qammod(data, M_data, ...
        'UnitAveragePower', true);

    Sp = zeros(Nfft,1);
    Sp(NumGuardSC_spectr+(1:NumSC_spectr)) = mod_data;

    Sig = ifft(ifftshift(Sp));

    % Сигнал без CP
    idx_noCP = (sym-1)*Nfft+1 : sym*Nfft;
    Sig_noCP_long(idx_noCP) = Sig;

    % Сигнал с CP
    idx_withCP = (sym-1)*(Nfft+CP_len)+1 : sym*(Nfft+CP_len);
    Sig_withCP_long(idx_withCP) = ...
        [Sig(end-CP_len+1:end); Sig];
end

% Частота дискретизации нормирована к единице
Fs_psd = 1;
nfft_psd = 4*Nfft;

[Pxx_noCP, f_noCP] = pwelch(Sig_noCP_long, [], [], nfft_psd, Fs_psd, 'centered');

[Pxx_withCP, f_withCP] = pwelch(Sig_withCP_long, [], [], nfft_psd, Fs_psd, 'centered');

figure(11);
clf;
hold on;
grid on;

plot(f_noCP, 10*log10(Pxx_noCP/max(Pxx_noCP)),'b-', 'LineWidth', 1.5, 'DisplayName', 'Без CP');

plot(f_withCP, 10*log10(Pxx_withCP/max(Pxx_withCP)), 'r--', 'LineWidth', 1.5,'DisplayName', sprintf('С CP (CP = %d)', CP_len));

xlabel('Нормированная частота (f/Fs)');
ylabel('Нормированная PSD, дБ');
title('Спектр OFDM-сигнала (метод Уэлча)');
legend('Location','best');
xlim([-0.5 0.5]);
ylim([-80 5]);
hold off;


%% ДОПОЛНИТЕЛЬНОЕ ЗАДАНИЕ 2: СПЕКТРАЛЬНЫЙ СОСТАВ ПОДНЕСУЩИХ

figure(12);
clf;
hold on;
grid on;

data = randi([0 M_data-1], NumSC_spectr, 1);
mod_data = qammod(data, M_data, 'UnitAveragePower', true);

Sp = zeros(Nfft,1);
Sp(NumGuardSC_spectr+(1:NumSC_spectr)) = mod_data;

amplitudes = abs(Sp);
f_subcarriers = (-Nfft/2:Nfft/2-1)/Nfft;

stem(f_subcarriers, amplitudes, 'b.', 'MarkerSize', 4);

xlabel('Нормированная частота (f/Fs)');
ylabel('|X(k)|');
title('Поднесущие OFDM-сигнала в частотной области');
xlim([-0.5 0.5]);
hold off;


%% ДОПОЛНИТЕЛЬНОЕ ЗАДАНИЕ 2: АКТИВНЫЕ ПОДНЕСУЩИЕ

figure(13);
clf;
hold on;
grid on;

active_idx_spectr = NumGuardSC_spectr+(1:NumSC_spectr);

stem(f_subcarriers(active_idx_spectr),amplitudes(active_idx_spectr), 'b.', 'MarkerSize', 3);

xlabel('Нормированная частота (f/Fs)');
ylabel('|X(k)|');
title(sprintf('Активные поднесущие (NumSC = %d)', NumSC_spectr));
hold off;


%% ДОПОЛНИТЕЛЬНОЕ ЗАДАНИЕ 2.1: ИХ И ЧХ КАНАЛОВ EPA, EVA, ETU

fprintf('\nДополнительное задание 2.1: ИХ и ЧХ каналов EPA, EVA, ETU\n');

channels_delays = {EPA_delays, EVA_delays, ETU_delays};
channels_gains = {EPA_gains_lin, EVA_gains_lin, ETU_gains_lin};

figure(14);
clf;

for ch_idx = 1:3

    tau = channels_delays{ch_idx};
    P_lin = channels_gains{ch_idx};

    tau_max_samples = round(tau(end)*Fs);
    h_len = max(tau_max_samples+1, 10);
    h_ch = zeros(1,h_len);

    for p = 1:length(tau)
        idx = round(tau(p)*Fs)+1;

        if idx <= h_len
            h_ch(idx) = h_ch(idx)+sqrt(P_lin(p));
        end
    end

    % Средняя и среднеквадратичная задержки
    tau_avg = sum(P_lin.*tau)/sum(P_lin);

    tau_rms = sqrt( ...
        sum(P_lin.*(tau-tau_avg).^2)/sum(P_lin));

    % Приближённая полоса когерентности
    Bc = 1/tau_rms;

    subplot(3,1,ch_idx);

    stem((0:h_len-1)/Fs*1e9, abs(h_ch),[colors_ch{ch_idx}, '.'],'MarkerSize', 6, 'LineWidth', 1.5);

    title(sprintf('Импульсная характеристика %s: ''\\tau_{rms}=%.2f нс, B_c=%.2f кГц',channel_names{ch_idx}, tau_rms*1e9, Bc/1e3));

    xlabel('Задержка, нс');
    ylabel('|h(t)|');
    xlim([0 tau(end)*1e9*1.2]);
    grid on;

    fprintf('Канал %s: tau_rms = %.2f нс, Bc = %.2f кГц\n',channel_names{ch_idx}, tau_rms*1e9, Bc/1e3);
end


%% ЧАСТОТНЫЕ ХАРАКТЕРИСТИКИ КАНАЛОВ EPA, EVA, ETU

figure(15);
clf;
hold on;
grid on;

for ch_idx = 1:3

    tau = channels_delays{ch_idx};
    P_lin = channels_gains{ch_idx};

    tau_max_samples = round(tau(end)*Fs);
    h_len = max(tau_max_samples+1, 10);
    h_ch = zeros(1,h_len);

    for p = 1:length(tau)
        idx = round(tau(p)*Fs)+1;

        if idx <= h_len
            h_ch(idx) = h_ch(idx)+sqrt(P_lin(p));
        end
    end

    Nfft_ch = 2^14;
    H = fftshift(fft(h_ch,Nfft_ch));

    f_ch = (-Nfft_ch/2:Nfft_ch/2-1)*(Fs/Nfft_ch);
    H_db = 20*log10(abs(H)+eps);

    tau_avg = sum(P_lin.*tau)/sum(P_lin);
    tau_rms = sqrt(sum(P_lin.*(tau-tau_avg).^2)/sum(P_lin));
    Bc = 1/tau_rms;

    plot(f_ch/1e6,H_db,[colors_ch{ch_idx},'-'],'LineWidth',1.5,'DisplayName',sprintf('%s: tau_rms=%.1f нс, Bc=%.1f кГц',channel_names{ch_idx},tau_rms*1e9,Bc/1e3));
end

xlabel('Частота, МГц');
ylabel('|H(f)|, дБ');
title('Частотные характеристики каналов EPA, EVA, ETU');
legend('Location','best');
xlim([-Fs/2/1e6 Fs/2/1e6]);
ylim([-60 5]);
hold off;

%% ДОП ЗАДАНИЕ 2.2: АНАЛИТИЧЕСКАЯ ОЦЕНКА СЕЛЕКТИВНОСТИ
fprintf('\nДОП ЗАДАНИЕ 2.2: Анализ селективности\n');

% Шаг между поднесущими
Delta_f = Fs / Nfft;       % = 1/T_sym
B_subcarrier = Delta_f;

% Полоса активных поднесущих
B_total = NumSC_fad * Delta_f;

% Длительность циклического префикса
T_CP = CP_len / Fs;

fprintf('Шаг между поднесущими: %.2f кГц\n', Delta_f/1e3);
fprintf('Полоса активных поднесущих: %.2f МГц\n', B_total/1e6);
fprintf('Длительность CP: %.2f мкс\n\n', T_CP*1e6);

for ch_idx = 1:3

    tau = channels_delays{ch_idx};
    P_lin = channels_gains{ch_idx};

    % Средняя задержка и среднеквадратичное рассеяние задержек
    tau_avg = sum(P_lin .* tau) / sum(P_lin);
    tau_rms = sqrt(sum(P_lin .* (tau - tau_avg).^2) / sum(P_lin));

    % Приближённая оценка полосы когерентности
    Bc = 1 / tau_rms;

    % Максимальная задержка профиля
    tau_max = max(tau);

    fprintf('Канал %s:\n', channel_names{ch_idx});
    fprintf('  Среднеквадратичное рассеяние задержек: %.2f нс\n', tau_rms*1e9);
    fprintf('  Оценка полосы когерентности B_c: %.2f МГц\n', Bc/1e6);
    fprintf('  Максимальная задержка: %.2f мкс\n', tau_max*1e6);

    % Частотная селективность по полной полосе OFDM
    if B_total > Bc
        fprintf('  B_total > B_c: канал частотно-селективный по полосе OFDM\n');
    else
        fprintf('  B_total <= B_c: канал приблизительно плоский по полосе OFDM\n');
    end

    % Приближённое условие плоского замирания на одной поднесущей
    if B_subcarrier < Bc
        fprintf('  Delta_f < B_c: замирания на отдельной поднесущей приблизительно плоские\n');
    else
        fprintf('  Delta_f >= B_c: возможна частотная селективность внутри поднесущей\n');
    end

    % Проверка достаточности циклического префикса
    if T_CP >= tau_max
        fprintf('  T_CP >= tau_max: CP достаточен для защиты от МСИ по задержкам профиля\n');
    else
        fprintf('  T_CP < tau_max: CP может быть недостаточен для устранения МСИ\n');
    end

    fprintf('\n');
end

%% ДОП ЗАДАНИЕ 2.3-2.4: BER В КАНАЛАХ С ЗАМИРАНИЯМИ (EPA, EVA, ETU)
fprintf('ДОП ЗАДАНИЕ 2.3-2.4: BER в каналах с замираниями\n');
NumGuardSC_fad = round((Nfft - NumSC_fad)/2);

pilot_step_fad = 16;
pilot_idx_fad = 1:pilot_step_fad:NumSC_fad;
data_idx_fad = setdiff(1:NumSC_fad, pilot_idx_fad);

nDataSymbols_fad = length(data_idx_fad);
nPilots_fad = length(pilot_idx_fad);

rayleighchan_EPA = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'PathDelays', EPA_delays, ...
    'AveragePathGains', EPA_gains_dB, ...
    'MaximumDopplerShift', 0, ...
    'NormalizePathGains', true);

rayleighchan_EVA = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'PathDelays', EVA_delays, ...
    'AveragePathGains', EVA_gains_dB, ...
    'MaximumDopplerShift', 0, ...
    'NormalizePathGains', true);

rayleighchan_ETU = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'PathDelays', ETU_delays, ...
    'AveragePathGains', ETU_gains_dB, ...
    'MaximumDopplerShift', 0, ...
    'NormalizePathGains', true);

channels_fad = {rayleighchan_EPA, rayleighchan_EVA, rayleighchan_ETU};

% Количество независимых реализаций канала
N_realizations = 10;

BER_fading = zeros(length(channels_fad), length(SNR_dB));

for ch_idx = 1:length(channels_fad)

    fprintf('Канал %s:\n', channel_names{ch_idx});

    for k = 1:length(SNR_dB)

        errorBits = 0;
        totalBits = 0;

        SNR_linear = 10^(SNR_dB(k)/10);

        % Усреднение по реализациям канала
        for realization = 1:N_realizations

            % Новое случайное состояние для каждой реализации
            rng(1000 + realization);

            % Сброс состояния канала перед новой реализацией
            reset(channels_fad{ch_idx});

            for sym = 1:NumOfdmSymbols

                % Информационные биты и QAM
                data_bits = randi([0 1], nDataSymbols_fad * log2(M_data), 1);
                mod_data = qammod(data_bits, M_data,'InputType', 'bit','UnitAveragePower', true);

                % Пилоты QPSK
                pilot_bits = randi([0 1],nPilots_fad * log2(M_pilot), 1);
                mod_pilots = qammod(pilot_bits, M_pilot, 'InputType', 'bit','UnitAveragePower', true);

                % Формирование OFDM-символа
                Sp = zeros(Nfft, 1);

                Sp(NumGuardSC_fad + data_idx_fad) = mod_data;
                Sp(NumGuardSC_fad + pilot_idx_fad) = mod_pilots;

                Sig = ifft(ifftshift(Sp));
                Sig_CP = [Sig(end-CP_len+1:end); Sig];

                % Передача через канал
                Rx_Signal = channels_fad{ch_idx}(Sig_CP);

                % Добавление АБГШ

                % Мощность принятого сигнала до добавления шума
                Ps = mean(abs(Rx_Signal).^2);

                N0 = Ps / SNR_linear;
                sigma_noise = sqrt(N0/2);
                noise = sigma_noise *(randn(size(Rx_Signal)) + 1j*randn(size(Rx_Signal)));
                Rx_Signal = Rx_Signal + noise;

                % Приём OFDM-символа
                Sig_rec = Rx_Signal(CP_len+1:CP_len+Nfft);
                Sp_rec = fftshift(fft(Sig_rec));

                % Оценка канала по пилотам
                H_pilots = Sp_rec(NumGuardSC_fad + pilot_idx_fad)./ mod_pilots;
                H_est = interp1(pilot_idx_fad, H_pilots, (1:NumSC_fad)', 'linear', 'extrap');
                Sp_eq = Sp_rec(NumGuardSC_fad + data_idx_fad) ./ H_est(data_idx_fad);

                % Демодуляция
                bits_rx = qamdemod(Sp_eq, M_data,'OutputType', 'bit','UnitAveragePower', true);

                % Накопление ошибок
                errorBits = errorBits + sum(data_bits ~= bits_rx);
                totalBits = totalBits + length(data_bits);

            end
        end

        % BER по всем реализациям
        BER_fading(ch_idx, k) = errorBits / totalBits;

    end

    fprintf('  BER при SNR = %g дБ: %e\n', SNR_dB(6), BER_fading(ch_idx, 6));
end

% Построение BER
figure(16);
clf;
hold on;
grid on;

markers_ch = {'o-', 's-', 'd-'};

for ch_idx = 1:length(channels_fad)

    semilogy(SNR_dB, BER_fading(ch_idx,:),[colors_ch{ch_idx} markers_ch{ch_idx}],'LineWidth', 1.5,'MarkerSize', 6,'DisplayName', sprintf('%s (с эквалайзером)', channel_names{ch_idx}));

end

semilogy(SNR_dB, berawgn(SNR_dB, 'qam', M_data), '--k', 'LineWidth', 2,'DisplayName', 'Теория AWGN (berawgn)');

xlabel('SNR (дБ)');
ylabel('BER');
title('Помехоустойчивость OFDM в каналах с замираниями (EPA, EVA, ETU)');
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
ylim([1e-5 1]);
hold off;
%% ДОП ЗАДАНИЕ 2.5: СОЗВЕЗДИЕ ДО И ПОСЛЕ ЭКВАЛАЙЗЕРА
fprintf('\nСозвездие до и после эквалайзера (канал EPA, SNR=15 дБ)\n');

SNR_test = 25;
SNR_linear_test = 10^(SNR_test/10);
release(channels_fad{1});
totalBits_const = NumOfdmSymbols * nDataSymbols_fad * log2(M_data); % общее количество переданных информационных битов 

data_bits = randi([0 1], nDataSymbols_fad * log2(M_data), NumOfdmSymbols);
mod_data = qammod(data_bits, M_data, 'InputType', 'bit', 'UnitAveragePower', true);

pilot_bits = randi([0 1], nPilots_fad * log2(M_pilot), NumOfdmSymbols);
mod_pilots = qammod(pilot_bits, M_pilot, 'InputType', 'bit', 'UnitAveragePower', true);
        
% Матрица OFDM-символов в частотной области
Sp = zeros(Nfft, NumOfdmSymbols);

Sp(NumGuardSC_fad + data_idx_fad, :) = mod_data;
Sp(NumGuardSC_fad + pilot_idx_fad, :) = mod_pilots;

Sig = ifft(ifftshift(Sp, 1), [], 1);
Sig_CP = [Sig(end-CP_len+1:end, :); Sig];

tx_signal = Sig_CP(:);
release(channels_fad{1});

Rx_Signal = channels_fad{1}(tx_signal);
N0 = Es / SNR_linear_test;
sigma_noise = sqrt(N0/2);
noise = sigma_noise * (randn(size(Rx_Signal)) + 1j*randn(size(Rx_Signal)));

Rx_Signal = Rx_Signal + noise;

Rx_matrix = reshape(Rx_Signal, Nfft + CP_len, NumOfdmSymbols);
Sig_rec = Rx_matrix(CP_len+1:CP_len+Nfft, :);

% Переход в частотную область
Sp_rec = fftshift(fft(Sig_rec, [], 1), 1);

constellation_before_eq = Sp_rec(NumGuardSC_fad + data_idx_fad, 1:min(3, NumOfdmSymbols));

% Оценка канала по пилотам для каждого OFDM-символа
H_pilots = Sp_rec(NumGuardSC_fad + pilot_idx_fad, :) ./ mod_pilots;

H_est = interp1(pilot_idx_fad, H_pilots, (1:NumSC_fad)', 'linear', 'extrap');

received_data = Sp_rec(NumGuardSC_fad + data_idx_fad, :);
H_data = H_est(data_idx_fad, :);

Sp_eq = received_data ./ H_data;

constellation_after_eq = Sp_eq(:, 1:min(3, NumOfdmSymbols));

bits_rx_matrix = qamdemod(Sp_eq, M_data, 'OutputType', 'bit', 'UnitAveragePower', true);

% Подсчёт BER
errorBits_const = sum(data_bits(:) ~= bits_rx_matrix(:));
BER_const = errorBits_const / totalBits_const;

fprintf('BER при SNR=%d дБ в канале EPA = %e\n', SNR_test, BER_const);

scatterplot(constellation_before_eq(:), 1, 0, 'b.');
title('Созвездие ДО эквалайзера (канал EPA, SNR=15 дБ)');
xlim([-3, 3]); ylim([-3, 3]);
grid on;

scatterplot(constellation_after_eq(:), 1, 0, 'r.');
title('Созвездие ПОСЛЕ эквалайзера (канал EPA, SNR=15 дБ)');
xlim([-3, 3]); ylim([-3, 3]);
grid on;

%% ДОП ЗАДАНИЕ 2.6: ДОПЛЕРОВСКИЕ СПЕКТРЫ EPA5, EVA70, ETU300
fprintf('\nДОП ЗАДАНИЕ 2.6: Доплеровские спектры\n');

fD_values_dop = [5 70 300];
N_symbols_doppler = [100000, 30000, 20000];

channel_names_dop = {'EPA5','EVA70','ETU300'};
delays_dop = {EPA_delays, EVA_delays, ETU_delays};
gains_dop = {EPA_gains_dB, EVA_gains_dB, ETU_gains_dB};
colors_dop = {'b','r','g'};

% Частота дискретизации процесса замираний:
% один отсчёт на OFDM-символ
Fs_doppler = 1 / T_sym;

% Параметры оценки спектра методом Welch
Nfft_doppler = 2^16;
window_len = 2048;
win = hamming(window_len);
noverlap = window_len / 2;

figure(19);
clf;
hold on;
grid on;

for idx = 1:length(fD_values_dop)

    fD = fD_values_dop(idx);
    N_sym = N_symbols_doppler(idx);

    % Канал с моделью спектра Джейкса
    chan = comm.RayleighChannel( ...
        'SampleRate', Fs_doppler, ...
        'PathDelays', delays_dop{idx}, ...
        'AveragePathGains', gains_dop{idx}, ...
        'MaximumDopplerShift', fD, ...
        'NormalizePathGains', true, ...
        'DopplerSpectrum', doppler('Jakes'), ...
        'PathGainsOutputPort', true);

    % Генерация независимой реализации процесса замираний
    reset(chan);
    [~, pathGains] = chan(ones(N_sym,1));

    % Коэффициент замираний первого пути
    h_dop = pathGains(:,1);

    % Оценка PSD методом Welch
    [Pxx, f] = pwelch(h_dop, win, noverlap, Nfft_doppler, Fs_doppler, 'centered');

    % Нормировка спектра в дБ
    Pxx_dB = 10*log10(Pxx / max(Pxx) + eps);

    % Отображаем область вокруг доплеровского спектра
    mask = abs(f) <= 1.5*fD;

    plot(f(mask), Pxx_dB(mask), colors_dop{idx},'LineWidth', 1.5,'DisplayName', sprintf('%s, f_D=%d Гц',channel_names_dop{idx}, fD));

end

xlabel('Доплеровская частота, Гц');
ylabel('Нормированная PSD, дБ');
title('Доплеровские спектры каналов (модель Джейкса)');
legend('Location','best');
ylim([-50 5]);
hold off;

fprintf('Доплеровские спектры построены.\n');

%% ДОП ЗАДАНИЕ 2.7: ВРЕМЯ КОГЕРЕНТНОСТИ
fprintf('\nДОП ЗАДАНИЕ 2.7: Время когерентности\n');

% Максимальные доплеровские частоты для профилей каналов
fD_values = [5 70 300];
channel_names_Tc = {'EPA5', 'EVA70', 'ETU300'};

% Время когерентности по формуле из лекции:
Tc = 1 ./ fD_values;

% Длительность полезной части OFDM-символа
T_sym = Nfft / Fs;

fprintf('Длительность полезной части OFDM-символа: %.3f ms\n\n', T_sym * 1000);

for idx = 1:length(fD_values)

    fprintf('%s:\n', channel_names_Tc{idx});
    fprintf('  Максимальная доплеровская частота: %d Hz\n', fD_values(idx));
    fprintf('  Время когерентности: %.3f ms\n', Tc(idx) * 1000);

    ratio = T_sym / Tc(idx);

    fprintf('  T_sym / t_corr = %.5f\n', ratio);

    if ratio < 0.1
        fprintf('  Канал можно считать медленно меняющимся за OFDM-символ.\n\n');
    else
        fprintf('  Изменение канала за OFDM-символ становится заметным.\n\n');
    end
end

%% ДОП ЗАДАНИЕ 2.8: MSE ОЦЕНКИ КАНАЛА ОТ ОСШ (EPA)

fprintf('\nДОП ЗАДАНИЕ 2.8: MSE оценки канала\n');

% Параметры OFDM
NumSC_mse = 1200;
NumGuardSC_mse = round((Nfft - NumSC_mse) / 2);

% Расстановка пилотов
pilot_step_mse = 16;
pilot_idx_mse = 1:pilot_step_mse:NumSC_mse;
data_idx_mse = setdiff(1:NumSC_mse, pilot_idx_mse);

nPilots_mse = length(pilot_idx_mse);

% Индексы активных поднесущих после fftshift
active_idx_mse = NumGuardSC_mse + (1:NumSC_mse);

% Количество OFDM-символов для усреднения
N_MSE_symbols = 500;

% Диапазон ОСШ
SNR_MSE_dB = 0:2:50;

% Массивы для хранения MSE
MSE_linear = zeros(length(SNR_MSE_dB), 1);
MSE_nearest = zeros(length(SNR_MSE_dB), 1);

% Формирование истинной импульсной характеристики канала EPA

% Перевод задержек путей в отсчёты
delay_samples_mse = round(EPA_delays * Fs);

% Импульсная характеристика канала
h_EPA_true = zeros(max(delay_samples_mse) + 1, 1);

% Фиксируем случайные фазы путей
rng(1);
path_phases_mse = 2*pi*rand(size(EPA_gains_lin));

for p = 1:length(delay_samples_mse)

    idx_h = delay_samples_mse(p) + 1;

    % Комплексная амплитуда пути
    path_gain = sqrt(EPA_gains_lin(p)) * exp(1j * path_phases_mse(p));

    % Суммирование путей, попавших в один отсчёт
    h_EPA_true(idx_h) = h_EPA_true(idx_h) + path_gain;

end

% Нормировка энергии импульсной характеристики
h_EPA_true = h_EPA_true / sqrt(sum(abs(h_EPA_true).^2));

% Истинный частотный отклик канала EPA
H_true_freq = fftshift(fft(h_EPA_true, Nfft));

% Истинный канал на активных поднесущих
H_true_active = H_true_freq(active_idx_mse);

% Расчёт MSE для разных значений ОСШ

for k = 1:length(SNR_MSE_dB)

    MSE_sum_linear = 0;
    MSE_sum_nearest = 0;

    SNR_linear = 10^(SNR_MSE_dB(k) / 10);

    for sym = 1:N_MSE_symbols

        % Формирование пилотного OFDM-символа
        pilot_bits = randi([0 1], nPilots_mse * log2(M_pilot), 1);
        mod_pilots = qammod(pilot_bits, M_pilot,'InputType', 'bit', 'UnitAveragePower', true);

        % OFDM-символ содержит пилоты на пилотных поднесущих
        Sp = zeros(Nfft, 1);
        Sp(active_idx_mse(pilot_idx_mse)) = mod_pilots;

        % Передача через канал EPA

        % Применение истинного частотного отклика канала
        Sp_freq_rx = Sp .* H_true_freq;

        % Переход во временную область
        Sig_rx = ifft(ifftshift(Sp_freq_rx));

        % Добавление циклического префикса
        Sig_rx_CP = [Sig_rx(end-CP_len+1:end); Sig_rx];

        % Добавление АБГШ
        Ps = mean(abs(Sig_rx_CP).^2);
        N0 = Ps / SNR_linear;
        sigma_noise = sqrt(N0 / 2);
        noise = sigma_noise *(randn(size(Sig_rx_CP)) + 1j * randn(size(Sig_rx_CP)));
        Rx_noisy = Sig_rx_CP + noise;

        % Приём OFDM-символа

        % Удаление циклического префикса
        Sig_rec = Rx_noisy(CP_len+1:CP_len+Nfft);

        % FFT и переход к порядку поднесущих после fftshift
        Sp_rec_noisy = fftshift(fft(Sig_rec));

        % LS-оценка канала на пилотах

        % Принятые пилоты делим на переданные
        H_est_pilots = Sp_rec_noisy(active_idx_mse(pilot_idx_mse))./ mod_pilots;

        % Интерполяция оценки канала

        % Линейная интерполяция
        H_est_linear = interp1(pilot_idx_mse, H_est_pilots,(1:NumSC_mse)', 'linear', 'extrap');

        % Интерполяция макима
        H_est_nearest = interp1(pilot_idx_mse, H_est_pilots,(1:NumSC_mse)', 'makima', 'extrap');

        % Расчёт MSE на информационных поднесущих
        H_true_data = H_true_active(data_idx_mse);

        H_linear_data = H_est_linear(data_idx_mse);
        H_nearest_data = H_est_nearest(data_idx_mse);

        MSE_sum_linear = MSE_sum_linear + mean(abs(H_linear_data - H_true_data).^2);

        MSE_sum_nearest = MSE_sum_nearest + mean(abs(H_nearest_data - H_true_data).^2);

    end

    % Усреднение по OFDM-символам
    MSE_linear(k) = MSE_sum_linear / N_MSE_symbols;
    MSE_nearest(k) = MSE_sum_nearest / N_MSE_symbols;

end

% Построение графика MSE
figure(20);
plot(SNR_MSE_dB, MSE_linear,'b-o', 'LineWidth', 1.5, 'MarkerSize', 6,'DisplayName', 'Линейная интерполяция');
hold on;

plot(SNR_MSE_dB, MSE_nearest, 'r-s', 'LineWidth', 1.5,'MarkerSize', 6,'DisplayName', 'Интерполяция Makima');

set(gca, 'YScale', 'log'); 
grid on;
xlabel('SNR, дБ');
ylabel('MSE');
legend('Location', 'best');
title('Зависимость MSE от SNR');

fprintf('MSE оценки канала для EPA рассчитаны.\n');
