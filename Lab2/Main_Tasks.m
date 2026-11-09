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

%% ПУНКТ 1: BER В КАНАЛЕ АБГШ ДЛЯ РАЗНОГО КОЛИЧЕСТВА ПОДНЕСУЩИХ

fprintf('ПУНКТ 1: BER в канале АБГШ\n');

NumSC_list = [200 800 1200 1600];
BER_AWGN = zeros(length(NumSC_list), length(SNR_dB));

for idxSC = 1:length(NumSC_list)

    NumSC = NumSC_list(idxSC);
    NumGuardSC = round((Nfft - NumSC)/2);

    N_bits = NumSC*log2(M_data);
    totalBits = NumOfdmSymbols*N_bits;

    for k = 1:length(SNR_dB)

        errorBits = 0;
        SNR_linear = 10^(SNR_dB(k)/10);

        for sym = 1:NumOfdmSymbols

            % Генерация и QAM-модуляция данных
            data_bits = randi([0 1], N_bits, 1);
            mod_data = qammod(data_bits, M_data,'InputType', 'bit', 'UnitAveragePower', true);

            % Размещение данных на активных поднесущих
            Sp = zeros(Nfft, 1);
            Sp(NumGuardSC + (1:NumSC)) = mod_data;

            % Формирование OFDM-символа
            Sig = ifft(ifftshift(Sp));

            % Добавление циклического префикса
            Sig_CP = [Sig(end-CP_len+1:end); Sig];

            % Расчёт мощности сигнала и дисперсии шума
            Ps = mean(abs(Sig_CP).^2);
            Pn = Ps/SNR_linear;

            noise = sqrt(Pn/2)*(randn(size(Sig_CP)) + 1j*randn(size(Sig_CP)));
            Rx_Signal = Sig_CP + noise;

            % Удаление CP и переход в частотную область
            Sig_rec = Rx_Signal(CP_len+1:CP_len+Nfft);
            Sp_rec = fftshift(fft(Sig_rec));

            % Демодуляция принятых данных
            demod_bits = qamdemod(Sp_rec(NumGuardSC+(1:NumSC)), M_data,'OutputType', 'bit', 'UnitAveragePower', true);

            errorBits = errorBits + sum(demod_bits(:) ~= data_bits(:));
        end

        BER_AWGN(idxSC, k) = errorBits/totalBits;
    end

    fprintf('NumSC = %4d: BER при SNR = 10 дБ: %e\n', NumSC, BER_AWGN(idxSC, 5));
end

% Теоретическая BER для 16-QAM в канале АБГШ
ber_theory = berawgn(SNR_dB, 'qam', M_data);

figure(1);
clf;
hold on;
grid on;

markers = {'o-', 's-', 'd-', '^-'};
colors = {'b', 'r', 'g', 'm'};

for idxSC = 1:length(NumSC_list)
    plot(SNR_dB, BER_AWGN(idxSC,:), ...
        [colors{idxSC} markers{idxSC}], ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf( ...
        'NumSC = %d (моделирование)', NumSC_list(idxSC)));
end

plot(SNR_dB, ber_theory, '--k','LineWidth', 2.5,'DisplayName', 'Теория (berawgn)');

xlabel('ОСШ, дБ');
ylabel('BER');
title('Пункт 1: Помехоустойчивость OFDM в канале АБГШ');
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
ylim([1e-5 1]);
hold off;

fprintf('\n');

%% ПУНКТ 2: BER БЕЗ ЗАЩИТНОГО ИНТЕРВАЛА ПРИ РАЗНОМ ШАГЕ ПИЛОТОВ

fprintf('ПУНКТ 2: BER без защитного интервала, разный шаг пилотов\n');

CP_len_p2 = 0;
pilot_steps_p2 = [8 16 32 64];

NumSC_p2 = 1200;
NumGuardSC_p2 = round((Nfft-NumSC_p2)/2);

% Канал EVA со статическими замираниями (fD = 0)
rayleighchan_p2 = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'PathDelays', EVA_delays, ...
    'AveragePathGains', EVA_gains_dB, ...
    'MaximumDopplerShift', 0, ...
    'NormalizePathGains', true);

BER_p2 = zeros(length(pilot_steps_p2), length(SNR_dB));

for idxP = 1:length(pilot_steps_p2)

    pilot_step_p2 = pilot_steps_p2(idxP);
    pilot_idx_p2 = 1:pilot_step_p2:NumSC_p2;
    data_idx_p2 = setdiff(1:NumSC_p2, pilot_idx_p2);

    nDataSymbols_p2 = length(data_idx_p2);
    nPilots_p2 = length(pilot_idx_p2);

    fprintf('Шаг пилотов = %d:\n', pilot_step_p2);

    for k = 1:length(SNR_dB)

        errorBits = 0;
        totalBits = NumOfdmSymbols*nDataSymbols_p2*log2(M_data);

        SNR_linear = 10^(SNR_dB(k)/10);

        for sym = 1:NumOfdmSymbols

            % Генерация информационных битов и QAM-модуляция
            data_bits = randi([0 1],nDataSymbols_p2*log2(M_data), 1);

            mod_data = qammod(data_bits, M_data,'InputType', 'bit', 'UnitAveragePower', true);

            % Генерация пилотов QPSK
            pilot_bits = randi([0 1],nPilots_p2*log2(M_pilot), 1);

            mod_pilots = qammod(pilot_bits, M_pilot,'InputType', 'bit', 'UnitAveragePower', true);

            % Формирование частотного OFDM-символа
            Sp = zeros(Nfft, 1);
            Sp(NumGuardSC_p2 + data_idx_p2) = mod_data;
            Sp(NumGuardSC_p2 + pilot_idx_p2) = mod_pilots;

            % Переход во временную область
            Sig = ifft(ifftshift(Sp));

            % CP = 0
            Sig_CP = Sig;

            % Передача через канал EVA
            release(rayleighchan_p2);
            Rx_Signal = rayleighchan_p2(Sig_CP);

            % Добавление АБГШ
            Ps = mean(abs(Rx_Signal).^2);
            Pn = Ps/SNR_linear;

            noise = sqrt(Pn/2)*(randn(size(Rx_Signal)) + 1j*randn(size(Rx_Signal)));

            Rx_Signal = Rx_Signal + noise;

            % Приём OFDM-символа
            Sig_rec = Rx_Signal(1:Nfft);
            Sp_rec = fftshift(fft(Sig_rec));

            % Оценка канала по пилотам
            H_pilots = Sp_rec(NumGuardSC_p2 + pilot_idx_p2)./ mod_pilots;

            % Линейная интерполяция оценки канала
            H_est = interp1(pilot_idx_p2, H_pilots,(1:NumSC_p2)', 'linear', 'extrap');

            % Эквализация информационных поднесущих
            Sp_eq = Sp_rec(NumGuardSC_p2 + data_idx_p2)./ H_est(data_idx_p2);

            % Демодуляция в биты
            bits_rx = qamdemod(Sp_eq, M_data,'OutputType', 'bit', 'UnitAveragePower', true);

            errorBits = errorBits + sum(data_bits ~= bits_rx);
        end

        BER_p2(idxP, k) = errorBits/totalBits;
    end

    fprintf('  BER при SNR = 10 дБ: %e\n', BER_p2(idxP, 5));
end

figure(2);
clf;
hold on;
grid on;

markers_p2 = {'o-', 's-', 'd-', '^-'};
colors_p2 = {'b', 'r', 'g', 'm'};

for idxP = 1:length(pilot_steps_p2)
    plot(SNR_dB, BER_p2(idxP,:), ...
        [colors_p2{idxP} markers_p2{idxP}], ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('Шаг пилотов = %d', ...
        pilot_steps_p2(idxP)));
end

plot(SNR_dB, berawgn(SNR_dB, 'qam', M_data), '--k','LineWidth', 2,'DisplayName', 'Теория AWGN');

xlabel('ОСШ, дБ');
ylabel('BER');
title('Пункт 2: BER без защитного интервала (CP = 0)');
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
ylim([1e-5 1]);
hold off;

%% ПУНКТ 3: СРАВНЕНИЕ ZP И CP

fprintf('\nПУНКТ 3: Сравнение ZP и CP\n');

CP_len_p3 = 100;% Длина защитного интервала
guard_types = {'CP', 'ZP'};

NumSC_p3 = 1200;
NumGuardSC_p3 = round((Nfft-NumSC_p3)/2);

pilot_step_p3 = 16;
pilot_idx_p3 = 1:pilot_step_p3:NumSC_p3;
data_idx_p3 = setdiff(1:NumSC_p3, pilot_idx_p3);

nDataSymbols_p3 = length(data_idx_p3);
nPilots_p3 = length(pilot_idx_p3);

% Канал EPA со статическими замираниями (fD = 0)
rayleighchan_p3 = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'PathDelays', EPA_delays, ...
    'AveragePathGains', EPA_gains_dB, ...
    'MaximumDopplerShift', 0, ...
    'NormalizePathGains', true);

BER_p3 = zeros(length(guard_types), length(SNR_dB));

for g_idx = 1:length(guard_types)

    fprintf('Тип защитного интервала: %s\n', guard_types{g_idx});

    for k = 1:length(SNR_dB)

        errorBits = 0;
        totalBits = NumOfdmSymbols*nDataSymbols_p3*log2(M_data);
        SNR_linear = 10^(SNR_dB(k)/10);

        for sym = 1:NumOfdmSymbols

            % Генерация информационных битов
            data_bits = randi([0 1],nDataSymbols_p3*log2(M_data), 1);

            mod_data = qammod(data_bits, M_data,'InputType', 'bit', 'UnitAveragePower', true);

            % Генерация пилотов QPSK
            pilot_bits = randi([0 1],nPilots_p3*log2(M_pilot), 1);

            mod_pilots = qammod(pilot_bits, M_pilot,'InputType', 'bit', 'UnitAveragePower', true);

            % Формирование частотного OFDM-символа
            Sp = zeros(Nfft, 1);
            Sp(NumGuardSC_p3 + data_idx_p3) = mod_data;
            Sp(NumGuardSC_p3 + pilot_idx_p3) = mod_pilots;

            % Переход во временную область
            Sig = ifft(ifftshift(Sp));

            % Добавление защитного интервала: CP или ZP
            if strcmp(guard_types{g_idx}, 'CP')
                Sig_guard = [Sig(end-CP_len_p3+1:end); Sig];
            else
                Sig_guard = [zeros(CP_len_p3, 1); Sig];
            end

            % Передача через канал EPA
            release(rayleighchan_p3);
            Rx_Signal = rayleighchan_p3(Sig_guard);

            % Добавление АБГШ
            Ps = mean(abs(Rx_Signal).^2);
            Pn = Ps/SNR_linear;

            noise = sqrt(Pn/2)*(randn(size(Rx_Signal)) + 1j*randn(size(Rx_Signal)));

            Rx_Signal = Rx_Signal + noise;

            % Удаление защитного интервала
            Sig_rec = Rx_Signal(CP_len_p3+1:CP_len_p3+Nfft);
            Sp_rec = fftshift(fft(Sig_rec));

            % Оценка канала по пилотам
            H_pilots = Sp_rec(NumGuardSC_p3 + pilot_idx_p3)./ mod_pilots;

            H_est = interp1(pilot_idx_p3, H_pilots,(1:NumSC_p3)', 'linear', 'extrap');

            % Эквализация информационных поднесущих
            Sp_eq = Sp_rec(NumGuardSC_p3 + data_idx_p3)./ H_est(data_idx_p3);

            % Демодуляция в биты
            bits_rx = qamdemod(Sp_eq, M_data,'OutputType', 'bit', 'UnitAveragePower', true);

            errorBits = errorBits + sum(data_bits ~= bits_rx);
        end

        BER_p3(g_idx, k) = errorBits/totalBits;
    end

    fprintf('  BER при SNR = 10 дБ: %e\n', BER_p3(g_idx, 5));
end

figure(3);
clf;
hold on;
grid on;

colors_p3 = {'b', 'r'};
markers_p3 = {'o-', 's-'};

for g_idx = 1:length(guard_types)
    plot(SNR_dB, BER_p3(g_idx,:), ...
        [colors_p3{g_idx} markers_p3{g_idx}], ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('%s (длина = %d)', ...
        guard_types{g_idx}, CP_len_p3));
end

plot(SNR_dB, berawgn(SNR_dB, 'qam', M_data), '--k','LineWidth', 2,'DisplayName', 'Теория AWGN');

xlabel('ОСШ, дБ');
ylabel('BER');
title(sprintf('Пункт 3: Сравнение CP и ZP (канал EPA, длина = %d)', ...
    CP_len_p3));
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
ylim([1e-5 1]);
hold off;
%% ПУНКТ 4: BER ДЛЯ РАЗНЫХ ДЛИН CP (МАТРИЦА Gi)

fprintf('\nПУНКТ 4: Разные длины CP\n');

CP_lengths_p4 = [0 25 50 100 150 200];

NumSC_p4 = 1200;
NumGuardSC_p4 = round((Nfft-NumSC_p4)/2);

pilot_step_p4 = 4;
pilot_idx_p4 = 1:pilot_step_p4:NumSC_p4;
data_idx_p4 = setdiff(1:NumSC_p4, pilot_idx_p4);

nDataSymbols_p4 = length(data_idx_p4);
nPilots_p4 = length(pilot_idx_p4);

NumOfdmSymbols_p4 = 500;

BER_p4 = zeros(length(CP_lengths_p4), length(SNR_dB));

% Импульсная характеристика канала ETU
delay_samples_p4 = round(ETU_delays*Fs);
path_power_p4 = 10.^(ETU_gains_dB/10);
path_power_p4 = path_power_p4/sum(path_power_p4);

rng(100);

path_coeff_p4 = sqrt(path_power_p4(:)/2).*(randn(length(path_power_p4),1) + 1j*randn(length(path_power_p4),1));

path_coeff_p4 = path_coeff_p4/norm(path_coeff_p4);

Lh_p4 = max(delay_samples_p4)+1;
h_p4 = zeros(Lh_p4,1);

for p = 1:length(delay_samples_p4)
    h_p4(delay_samples_p4(p)+1) = h_p4(delay_samples_p4(p)+1) + path_coeff_p4(p);
end

% Моделирование для разных длин CP
for ci = 1:length(CP_lengths_p4)

    CP_len_p4 = CP_lengths_p4(ci);
    block_len_p4 = Nfft + CP_len_p4;

    % Матрица свёртки Gi для канала ETU
    Gi = zeros(block_len_p4+Lh_p4-1, block_len_p4);

    for r = 1:size(Gi,1)
        for c = 1:block_len_p4
            d = r-c+1;

            if d >= 1 && d <= Lh_p4
                Gi(r,c) = h_p4(d);
            end
        end
    end

    fprintf('CP = %d\n', CP_len_p4);

    for k = 1:length(SNR_dB)

        errors = 0;
        totalBits = NumOfdmSymbols_p4*nDataSymbols_p4*log2(M_data);
        SNRlin = 10^(SNR_dB(k)/10);

        % Состояние канала сохраняется между OFDM-символами
        zi = zeros(Lh_p4-1,1);

        for sym = 1:NumOfdmSymbols_p4

            % Генерация информационных битов
            data_bits = randi([0 1],nDataSymbols_p4*log2(M_data), 1);

            mod_data = qammod(data_bits, M_data,'InputType', 'bit', 'UnitAveragePower', true);

            % Генерация пилотов
            pilot_bits = randi([0 1], nPilots_p4*log2(M_pilot), 1);

            mod_pilots = qammod(pilot_bits, M_pilot,'InputType', 'bit', 'UnitAveragePower', true);

            % Формирование частотного OFDM-символа
            Sp = zeros(Nfft,1);
            Sp(NumGuardSC_p4+data_idx_p4) = mod_data;
            Sp(NumGuardSC_p4+pilot_idx_p4) = mod_pilots;

            Sig = ifft(ifftshift(Sp));

            % Добавление CP
            if CP_len_p4 > 0
                Sig_CP = [Sig(end-CP_len_p4+1:end); Sig];
            else
                Sig_CP = Sig;
            end

            % Передача через канал ETU с памятью
            [Rx_clean, zi] = filter(h_p4, 1, Sig_CP, zi);

            % Добавление АБГШ согласно заданному ОСШ
            N0 = Es/SNRlin;
            noise = sqrt(N0/2)*(randn(size(Rx_clean)) + 1j*randn(size(Rx_clean)));

            Rx = Rx_clean + noise;

            % Удаление CP и FFT
            if CP_len_p4 > 0
                Sig_rec = Rx(CP_len_p4+1:CP_len_p4+Nfft);
            else
                Sig_rec = Rx(1:Nfft);
            end

            Sp_rec = fftshift(fft(Sig_rec));

            % Оценка канала по пилотам
            H_pilots = Sp_rec(NumGuardSC_p4+pilot_idx_p4)./ mod_pilots;

            H_est = interp1(pilot_idx_p4, H_pilots,(1:NumSC_p4)', 'linear', 'extrap');

            Y_data = Sp_rec(NumGuardSC_p4+data_idx_p4);
            H_data = H_est(data_idx_p4);

            Sp_eq = Y_data.*conj(H_data)./(abs(H_data).^2 + Nfft*N0);

            % Демодуляция в биты
            bits_rx = qamdemod(Sp_eq, M_data,'OutputType', 'bit', 'UnitAveragePower', true);

            errors = errors + sum(data_bits ~= bits_rx);
        end

        BER_p4(ci,k) = errors/totalBits;
    end

    % Вывод ближайшего значения ОСШ к 10 дБ
    [~, idx10] = min(abs(SNR_dB-10));
    fprintf('BER при SNR = %.1f дБ: %e\n', SNR_dB(idx10), BER_p4(ci,idx10));
end

% График BER для разных длин CP
figure(4);
clf;
hold on;
grid on;

colors_p4 = {'b','r','g','m','c','k'};
markers_p4 = {'o-','s-','d-','^-','x-','p-'};

for ci = 1:length(CP_lengths_p4)

    b = BER_p4(ci,:);
    b(b==0) = NaN;

    plot(SNR_dB, b, [colors_p4{ci} markers_p4{ci}],'LineWidth', 1.5,'MarkerSize', 6,'DisplayName', sprintf('CP = %d', CP_lengths_p4(ci)));
end

plot(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k', 'LineWidth', 2,'DisplayName', 'Теория AWGN');

xlabel('ОСШ, дБ');
ylabel('BER');
title('Пункт 4: BER для разных длин CP (канал ETU)');
legend('Location','southwest');
set(gca,'YScale','log');
ylim([1e-5 1]);
hold off;


%% ПУНКТ 5: РАЗНЫЕ ЗНАЧЕНИЯ fD ПРИ ФИКСИРОВАННОЙ ДЛИНЕ CP

fprintf('\nПУНКТ 5: Разные значения fD при фиксированной длине CP\n');

fD_values_p5 = [0 50 100 200 500 1000];
CP_len_p5 = 512;

NumSC_p5 = 1200;
NumGuardSC_p5 = round((Nfft-NumSC_p5)/2);

pilot_step_p5 = 8;
pilot_idx_p5 = 1:pilot_step_p5:NumSC_p5;
data_idx_p5 = setdiff(1:NumSC_p5, pilot_idx_p5);

nDataSymbols_p5 = length(data_idx_p5);
nPilots_p5 = length(pilot_idx_p5);

NumOfdmSymbols_p5 = 500;

BER_p5 = zeros(length(fD_values_p5), length(SNR_dB));

for fD_idx = 1:length(fD_values_p5)

    fD_p5 = fD_values_p5(fD_idx);
    fprintf('fD = %d Гц:\n', fD_p5);

    % Канал EVA с заданной доплеровской частотой
    rayleighchan_p5 = comm.RayleighChannel( ...
        'SampleRate', Fs, ...
        'PathDelays', EVA_delays, ...
        'AveragePathGains', EVA_gains_dB, ...
        'MaximumDopplerShift', fD_p5, ...
        'NormalizePathGains', true);

    for k = 1:length(SNR_dB)

        errorBits = 0;
        totalBits = NumOfdmSymbols_p5*nDataSymbols_p5*log2(M_data);
        SNR_linear = 10^(SNR_dB(k)/10);

        for sym = 1:NumOfdmSymbols_p5

            % Информационные биты
            data_bits = randi([0 1],nDataSymbols_p5*log2(M_data), 1);

            mod_data = qammod(data_bits, M_data,'InputType', 'bit', 'UnitAveragePower', true);

            % Пилоты QPSK
            pilot_bits = randi([0 1],nPilots_p5*log2(M_pilot), 1);

            mod_pilots = qammod(pilot_bits, M_pilot,'InputType', 'bit', 'UnitAveragePower', true);

            % Формирование OFDM-символа
            Sp = zeros(Nfft,1);
            Sp(NumGuardSC_p5+data_idx_p5) = mod_data;
            Sp(NumGuardSC_p5+pilot_idx_p5) = mod_pilots;

            Sig = ifft(ifftshift(Sp));
            Sig_CP = [Sig(end-CP_len_p5+1:end); Sig];

            % Передача через канал EVA
            Rx_Signal = rayleighchan_p5(Sig_CP);

            % Добавление АБГШ
            Ps = mean(abs(Rx_Signal).^2);
            N0 = Ps/SNR_linear;

            noise = sqrt(N0/2)*(randn(size(Rx_Signal)) + 1j*randn(size(Rx_Signal)));
            Rx_Signal = Rx_Signal + noise;

            % Приём OFDM-символа
            Sig_rec = Rx_Signal(CP_len_p5+1:CP_len_p5+Nfft);
            Sp_rec = fftshift(fft(Sig_rec));

            % Оценка канала по пилотам
            H_pilots = Sp_rec(NumGuardSC_p5+pilot_idx_p5)./ mod_pilots;

            H_est = interp1(pilot_idx_p5, H_pilots,(1:NumSC_p5)', 'linear', 'extrap');

            % Эквализация
            Sp_eq = Sp_rec(NumGuardSC_p5+data_idx_p5) ./ H_est(data_idx_p5);

            % Демодуляция в биты
            bits_rx = qamdemod(Sp_eq, M_data,'OutputType', 'bit', 'UnitAveragePower', true);

            errorBits = errorBits + sum(data_bits ~= bits_rx);
        end

        BER_p5(fD_idx,k) = errorBits/totalBits;
    end

    [~, idx10] = min(abs(SNR_dB-10));
    fprintf('  BER при SNR = %.1f дБ: %e\n', SNR_dB(idx10), BER_p5(fD_idx,idx10));

    release(rayleighchan_p5);
end

figure(5);
clf;
hold on;
grid on;

colors_p5 = {'b','r','g','m','c','y'};
markers_p5 = {'o-','s-','d-','^-','x-','p-'};

for fD_idx = 1:length(fD_values_p5)
    plot(SNR_dB, BER_p5(fD_idx,:), ...
        [colors_p5{fD_idx} markers_p5{fD_idx}], ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('fD = %d Гц', fD_values_p5(fD_idx)));
end

plot(SNR_dB, berawgn(SNR_dB,'qam',M_data), '--k','LineWidth', 2,'DisplayName', 'Теория AWGN');

xlabel('ОСШ, дБ');
ylabel('BER');
title(sprintf('Пункт 5: BER для разных fD (CP = %d, канал EVA)', CP_len_p5));
legend('Location','southwest');
set(gca,'YScale','log');
ylim([1e-5 1]);
hold off;