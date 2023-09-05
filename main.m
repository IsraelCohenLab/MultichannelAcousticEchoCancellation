%% Clear previous data.
clear all;
close all;
clc;
%% General Parameters.
fs= 16000; % If this is changed or the .wav files are changed, might need to change and downsample / sample more of all signals. Also, might need to change simulink files.
Ts= 1 / fs;
overlap= 0.75;
win_length= 512;
lambda= 0.35; 
beta= 5;
Short_Time_Average_Frames= 100;
STFT_jump= (1 - overlap) * win_length;
freqs= win_length / 2 + 1;
f_vec= (linspace(0, fs / 2, freqs))';
Segment_time= 10;
segment_length= Segment_time * fs;
Forefront_frames= overlap / (1 - overlap);
Correlation_Estimation_Last_Frame= Segment_time * fs / STFT_jump - Forefront_frames;
%% Simulated Room Parameters.
mex -setup C++
mex rir_generator.cpp rir_generator_core.cpp
%% Create Signals for Simulated Room Experiment.
[far_sig_1, ~]= audioread('SA1.wav'); 
[far_sig_2, ~]= audioread('SA2_FEMALE.wav');
[far_sig_3, ~]= audioread('SI1271.wav');
[near_sig_1, ~]= audioread('SA2_MALE.wav'); 
[near_sig_2, ~]= audioread('SI564.wav');
[near_sig_3, ~]= audioread('SI1194.wav');
[near_sig_4, ~]= audioread('SI1824.wav'); % If fs changes might need to decimate signals.
x= [far_sig_1; far_sig_2; far_sig_3];
s= [near_sig_1; near_sig_2; near_sig_3; near_sig_4];
x= x(1 : segment_length);
s= s(1 : segment_length);
audiowrite('X.wav', x, fs); 
audiowrite('S.wav', s, fs);
silence= zeros(segment_length, 1);
Simulink_sim_time= Segment_time - Ts;
Simulink_time_vec= transpose(0 : Ts : Simulink_sim_time);
x_forsimulink=  [Simulink_time_vec, x];
open('acoustical_lib.slx');
sim('NonLinear_System.slx', Simulink_sim_time);
x_NL(:, 1)= ans.x_NL(1, 1, 1 : segment_length);
x_NL_First_Loudspeaker= [x_NL; x_NL ; x_NL; silence ; silence];
x_NL_Second_Loudspeaker= [silence; silence; silence; x_NL; x_NL];
s_First_Speaker=[silence; s; silence; s; silence];
s_Second_Speaker= [silence; silence; s; silence; s];
x_NL= x_NL_First_Loudspeaker + x_NL_Second_Loudspeaker;
x= [x; x; x; x; x];
%% Plot Signals for M=L= 4.
M= 4;
L= 4;
[y, u, v, d] = Generate_Signals(M, x_NL_First_Loudspeaker, x_NL_Second_Loudspeaker, s_First_Speaker, s_Second_Speaker, fs, segment_length);
[D, Y, U, X, V, frames]= Generate_STFT_Signals(d, y, u, x, v, win_length, STFT_jump, beta, freqs, M);
Last_Time= (overlap * win_length + STFT_jump * frames)/ fs;
SER_First_Loudspeaker_First_Speaker= mag2db(rms(u((segment_length + 1) : (2 * segment_length), 1)) / rms(y((segment_length + 1) : (2 * segment_length), 1)))
SER_First_Loudspeaker_Second_Speaker= mag2db(rms(u((2 * segment_length + 1) : (3 * segment_length), 1)) / rms(y((2 * segment_length + 1) : (3 * segment_length), 1)))
SER_Second_Loudspeaker_First_Speaker= mag2db(rms(u((3 * segment_length + 1) : (4 * segment_length), 1)) / rms(y((3 * segment_length + 1) : (4 * segment_length), 1)))
SER_Second_Loudspeaker_Second_Speaker= mag2db(rms(u((4 * segment_length + 1) : (5 * segment_length), 1)) / rms(y((4 * segment_length + 1) : (5 * segment_length), 1)))
[~, ~, ~, ~, u_Hat_Spatial_LCMV]= Spatial_LCMV_Method(fs, frames, win_length, STFT_jump, beta, M, L, X, D, U, Y, V); 
audiowrite('U_Hat_Proposed_M=4_L=4.wav', u_Hat_Spatial_LCMV, fs);
figure;
subplot(4, 1, 1);
plot(linspace(0, Last_Time, length(d(:, 1))) - Segment_time, d(:, 1), 'LineWidth', 0.001);
xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
set(gca, 'FontSize', 7.5);
title('(a)', 'Interpreter', 'latex');
xline(Segment_time, 'c', 'LineWidth', 2);
xline(2 * Segment_time, 'c', 'LineWidth', 2);
xline(3 * Segment_time, 'c', 'LineWidth', 2);
xlim([0, Last_Time - Segment_time]);
ylim([-0.15, 0.15]);
yticks([-0.15, 0.15]);
subplot(4, 1, 2);
plot(linspace(0, Last_Time, length(y(:, 1))) - Segment_time, y(:, 1), 'LineWidth', 0.001);
xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
set(gca, 'FontSize', 7.5);
title('(b)', 'Interpreter', 'latex');
xline(Segment_time, 'c', 'LineWidth', 2);
xline(2 * Segment_time, 'c', 'LineWidth', 2);
xline(3 * Segment_time, 'c', 'LineWidth', 2);
xlim([0, Last_Time - Segment_time]);
ylim([-0.15, 0.15]);
yticks([-0.15, 0.15]);
subplot(4, 1, 3);
plot(linspace(0, Last_Time, length(u(:, 1))) - Segment_time, u(:, 1), 'LineWidth', 0.001);
xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
set(gca, 'FontSize', 7.5);
title('(c)', 'Interpreter', 'latex');
xline(Segment_time, 'c', 'LineWidth', 2);
xline(2 * Segment_time, 'c', 'LineWidth', 2);
xline(3 * Segment_time, 'c', 'LineWidth', 2);
xlim([0, Last_Time - Segment_time]);
ylim([-0.15, 0.15]);
yticks([-0.15, 0.15]);
subplot(4, 1, 4);
plot(linspace(0, Last_Time, length(u_Hat_Spatial_LCMV)) - Segment_time, u_Hat_Spatial_LCMV, 'LineWidth', 0.001);
xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
set(gca, 'FontSize', 7.5);
title('(d)', 'Interpreter', 'latex');
xline(Segment_time, 'c', 'LineWidth', 2);
xline(2 * Segment_time, 'c', 'LineWidth', 2);
xline(3 * Segment_time, 'c', 'LineWidth', 2);
xlim([0, Last_Time - Segment_time]);
ylim([-0.15, 0.15]);
yticks([-0.15, 0.15]);
set(gcf, 'PaperUnits', 'inches', 'Units', 'inches');
set(gcf, 'PaperPosition', [0 0 7.9 5]); 
set(gcf, 'PaperSize', [7.9 5]); 
saveas(gcf, 'Signals_Time_Domain', 'pdf');
%% Run as function of M.
figure;
L= 4;
for M= 2 : 5
    [y, u, v, d] = Generate_Signals(M, x_NL_First_Loudspeaker, x_NL_Second_Loudspeaker, s_First_Speaker, s_Second_Speaker, fs, segment_length);
    [D, Y, U, X, V, frames]= Generate_STFT_Signals(d, y, u, x, v, win_length, STFT_jump, beta, freqs, M);
    Last_Time= (overlap * win_length + STFT_jump * frames)/ fs;
    [ERLE_Spatial_LCMV, DI_Spatial_LCMV, ~, ~, u_Hat_Spatial_LCMV]= Spatial_LCMV_Method(fs, frames, win_length, STFT_jump, beta, ...
        M, L, X, D, U, Y, V); 
    audiowrite(strcat('U_Hat_Proposed_M=', num2str(M), '_L=4.wav'), u_Hat_Spatial_LCMV, fs);
    subplot(2, 1, 1);
    plot(linspace(0, Last_Time, length(ERLE_Spatial_LCMV)) - Segment_time, pow2db(ERLE_Spatial_LCMV),  ...
        'LineWidth', 0.001);
    hold on;
    xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
    ylabel('[dB]', 'Interpreter', 'latex');
    set(gca, 'FontSize', 7.5);
    title('(a)', 'Interpreter', 'latex');
    xline(Segment_time, 'c', 'LineWidth', 2);
    xline(2 * Segment_time, 'c', 'LineWidth', 2);
    xline(3 * Segment_time, 'c', 'LineWidth', 2);
    xlim([0, Last_Time - Segment_time]);
    ylim([0, 40]);
    subplot(2, 1, 2);
    plot(linspace(0, Last_Time, length(DI_Spatial_LCMV)) - Segment_time, pow2db(DI_Spatial_LCMV), 'LineWidth', 0.001);
    hold on;
    xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
    ylabel('[dB]', 'Interpreter', 'latex');
    set(gca, 'FontSize', 7.5);
    title('(b)', 'Interpreter', 'latex');
    xline(Segment_time, 'c', 'LineWidth', 2);
    xline(2 * Segment_time, 'c', 'LineWidth', 2);
    xline(3 * Segment_time, 'c', 'LineWidth', 2);
    xlim([0, Last_Time - Segment_time]);
    ylim([-20, 15]);
end
set(gcf, 'PaperUnits', 'inches', 'Units', 'inches');
set(gcf, 'PaperPosition', [0 0 7.9 5]); 
set(gcf, 'PaperSize', [7.9 5]); 
saveas(gcf, 'ERLE_DI_Varying_M', 'pdf');
%% Run as Function of L.
figure;
M= 4;
[y, u, v, d] = Generate_Signals(M, x_NL_First_Loudspeaker, x_NL_Second_Loudspeaker, s_First_Speaker, s_Second_Speaker, fs, segment_length);
[D, Y, U, X, V, frames]= Generate_STFT_Signals(d, y, u, x, v, win_length, STFT_jump, beta, freqs, M);
Last_Time= (overlap * win_length + STFT_jump * frames)/ fs;
for L= 2 : 5 
    [ERLE_Spatial_LCMV, DI_Spatial_LCMV, ~, ~, u_Hat_Spatial_LCMV]= Spatial_LCMV_Method(fs, frames, win_length, STFT_jump, beta, ...
        M, L, X, D, U, Y, V); 
    audiowrite(strcat('U_Hat_Proposed_M=4_L=', num2str(L), '.wav'), u_Hat_Spatial_LCMV, fs);
    subplot(2, 1, 1);
    plot(linspace(0, Last_Time, length(ERLE_Spatial_LCMV)) - Segment_time, pow2db(ERLE_Spatial_LCMV),  ...
        'LineWidth', 0.001);
    hold on;
    xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
    ylabel('[dB]', 'Interpreter', 'latex');
    set(gca, 'FontSize', 7.5);
    title('(a)', 'Interpreter', 'latex');
    xline(Segment_time, 'c', 'LineWidth', 2);
    xline(2 * Segment_time, 'c', 'LineWidth', 2);
    xline(3 * Segment_time, 'c', 'LineWidth', 2);
    xlim([0, Last_Time - Segment_time]);
    ylim([0, 40]);
    subplot(2, 1, 2);
    plot(linspace(0, Last_Time, length(DI_Spatial_LCMV)) - Segment_time, pow2db(DI_Spatial_LCMV), 'LineWidth', 0.001);
    hold on;
    xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
    ylabel('[dB]', 'Interpreter', 'latex');
    set(gca, 'FontSize', 7.5);
    title('(b)', 'Interpreter', 'latex');
    xline(Segment_time, 'c', 'LineWidth', 2);
    xline(2 * Segment_time, 'c', 'LineWidth', 2);
    xline(3 * Segment_time, 'c', 'LineWidth', 2);
    xlim([0, Last_Time - Segment_time]);
    ylim([-20, 15]);
end
set(gcf, 'PaperUnits', 'inches', 'Units', 'inches');
set(gcf, 'PaperPosition', [0 0 7.9 5]); 
set(gcf, 'PaperSize', [7.9 5]); 
saveas(gcf, 'ERLE_DI_Varying_L', 'pdf');
%% Compare with Baseline.
figure;
L=4;
M=4;
[y, u, v, d] = Generate_Signals(M, x_NL_First_Loudspeaker, x_NL_Second_Loudspeaker, s_First_Speaker, s_Second_Speaker, fs, segment_length);
[D, Y, U, X, V, frames]= Generate_STFT_Signals(d, y, u, x, v, win_length, STFT_jump, beta, freqs, M);
Last_Time= (overlap * win_length + STFT_jump * frames)/ fs;
[ERLE_Spatial_LCMV, DI_Spatial_LCMV, ~, ~, u_Hat_Spatial_LCMV]= Spatial_LCMV_Method(fs, frames, win_length, STFT_jump, beta, M, ...
    L, X, D, U, Y, V); 
Y_tilde=Get_Y_Tilde_by_NLMS(d, fs, beta, win_length, overlap, freqs, X, D(:, :, 1), frames, Correlation_Estimation_Last_Frame);
[ERLE_Baseline, DI_Baseline, ~, ~, u_Hat_Baseline]= Baseline(fs, frames, win_length, STFT_jump, beta, L, D(:, :, 1), U(:, :, 1), Y(:, :, 1), ...
    V(:, :, 1), Y_tilde, lambda, Short_Time_Average_Frames);
audiowrite(strcat('U_Hat_Baseline.wav'), u_Hat_Baseline, fs);
audiowrite(strcat('U_Hat_Proposed.wav'), u_Hat_Spatial_LCMV, fs);
subplot(2, 1, 1);
plot(linspace(0, Last_Time, length(ERLE_Spatial_LCMV))- Segment_time, pow2db(ERLE_Spatial_LCMV),  ...
    'LineWidth', 0.001);
hold on;
plot(linspace(0, Last_Time, length(ERLE_Baseline))- Segment_time, pow2db(ERLE_Baseline),  ...
    'LineWidth', 0.001);
xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
ylabel('[dB]', 'Interpreter', 'latex');
set(gca, 'FontSize', 7.5);
title('(a)', 'Interpreter', 'latex');
xline(Segment_time, 'c', 'LineWidth', 2);
xline(2 * Segment_time, 'c', 'LineWidth', 2);
xline(3 * Segment_time, 'c', 'LineWidth', 2);
xlim([0, Last_Time - Segment_time]);
ylim([0, 40]);
subplot(2, 1, 2);
plot(linspace(0, Last_Time, length(DI_Spatial_LCMV)) - Segment_time, pow2db(DI_Spatial_LCMV), 'LineWidth', 0.001);
hold on;
plot(linspace(0, Last_Time, length(DI_Baseline))- Segment_time, pow2db(DI_Baseline), 'LineWidth', 0.001);
hold on;
xlabel('$t\left(s\right)$', 'Interpreter', 'latex');
ylabel('[dB]', 'Interpreter', 'latex');
set(gca, 'FontSize', 7.5);
title('(b)', 'Interpreter', 'latex');
xline(Segment_time, 'c', 'LineWidth', 2);
xline(2 * Segment_time, 'c', 'LineWidth', 2);
xline(3 * Segment_time, 'c', 'LineWidth', 2);
xlim([0, Last_Time - Segment_time]);
ylim([-20, 15]);
set(gcf, 'PaperUnits', 'inches', 'Units', 'inches');
set(gcf, 'PaperPosition', [0 0 7.9 5]); 
set(gcf, 'PaperSize', [7.9 5]); 
saveas(gcf, 'ERLE_DI_Competing_Methods', 'pdf');