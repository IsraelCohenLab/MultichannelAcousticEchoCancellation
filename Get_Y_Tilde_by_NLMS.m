function Y_Tilde = Get_Y_Tilde_by_NLMS(d, fs, beta, win_length, overlap, freqs, X, D, frames, Last_NLMS_Adaptation_Frame)
%% This Function Finds the initila guess for U With the NLMS Algorithm.
%% Create Y_tilde.
Y_Tilde= zeros(freqs, frames);
audiowrite('temp_in.wav', d, fs);
fin= 'temp_in.wav';
fout= 'temp_out';
[~, ~, sigma_v_squared]= omlsa_changed_by_Yuval(fin, fout, beta, fs, win_length, overlap);
sigma_e_squared= zeros(freqs, 1);
sigma_x_squared= zeros(freqs, 1);
K= 2;
Q_tag= 20;
L_NLMS= (2  * K + 1) * Q_tag;
h_NLMS= zeros(L_NLMS, freqs);
for n= Q_tag : frames
    for k= 1 : freqs
        freq_indeces= mod((k - K - 1 : k + K - 1), freqs) + 1;
        frame_indeces= n : -1 : (n - Q_tag + 1);
        x_NLMS= reshape(X(freq_indeces, frame_indeces), L_NLMS, 1);
        if(n <= Last_NLMS_Adaptation_Frame)
            [h_NLMS(:, k), sigma_e_squared(k), sigma_x_squared(k)]= Update_Freq_Domain_NLMS_Filter(h_NLMS(:, k), sigma_e_squared(k), sigma_x_squared(k), x_NLMS, ...
            sigma_v_squared(k, n), D(k, n), Q_tag, K);
        else
            [~, sigma_e_squared(k), sigma_x_squared(k)]= Update_Freq_Domain_NLMS_Filter(h_NLMS(:, k), sigma_e_squared(k), sigma_x_squared(k), x_NLMS, sigma_v_squared(k, n), D(k, n), Q_tag, K);
        end
        Y_Tilde(k, n)= h_NLMS(:, k)' * x_NLMS;
    end
end
