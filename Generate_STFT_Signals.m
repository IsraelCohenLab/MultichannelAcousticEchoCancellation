function [D, Y, U, X, V, frames] = Generate_STFT_Signals(d, y, u, x, v, win_length, STFT_jump, beta, freqs, Mics)
%% This function finds All the STFT Signals.
%% Calculate STFTs.
X= stft_changed_by_Yuval(x, win_length, STFT_jump, 1, beta);
for m= 1 : Mics
    D(:, :, m)= stft_changed_by_Yuval(d(:, m), win_length, STFT_jump, 1, beta);
    Y(:, :, m)= stft_changed_by_Yuval(y(:, m), win_length, STFT_jump, 1, beta);
    U(:, :, m)= stft_changed_by_Yuval(u(:, m), win_length, STFT_jump, 1, beta);
    V(:, :, m)= stft_changed_by_Yuval(v(:, m), win_length, STFT_jump, 1, beta);
end
X= X(1 : freqs, :);
D= D(1 : freqs, :, :);
Y= Y(1 : freqs, :, :);
U= U(1 : freqs, :, :);
V= V(1 : freqs, :, :);
%% Find the Number of Frames.
frames= size(D, 2);
end