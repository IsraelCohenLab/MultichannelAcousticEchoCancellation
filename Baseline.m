function [ERLE, DI, NRF, U_Hat, u_hat] = Baseline(fs, frames, win_length, STFT_jump, beta, L, D, U, Y, V, Y_tilde, ...
    lambda, Short_Time_Average_Frames)
%% This function implements the temporal MVDR filter minus the corresponding NLMS output.
%% Define Parameters.
freqs= win_length / 2 + 1;
Summation_Last_Frame= Short_Time_Average_Frames + L - 1;
%% Implements the total MVDR beamformer.
Current_phi_d= zeros(freqs, L, L);
Current_phi_y_tilde= zeros(freqs, L, L);
Phi_d= zeros(freqs, L, L);
Phi_y_tilde= zeros(freqs, L, L);
Phi_u_c= zeros(freqs, L, L);
U_f= zeros(freqs, frames);
Y_re= zeros(freqs, frames);
V_re= zeros(freqs, frames);
U_Hat= zeros(freqs, frames);
for n= L : frames
    First_Frame= n - L + 1;  
    d_vec= transpose(D(:, n : -1 : First_Frame));
    y_tilde_vec= transpose(Y_tilde(:, n : -1 : First_Frame));
    y_vec= transpose(Y(:, n : -1 : First_Frame));
    u_vec= transpose(U(:, n : -1 : First_Frame));
    v_vec= transpose(V(:, n : -1 : First_Frame));
    for k= 1 : freqs
        Current_phi_d(k, :, :)= d_vec(:, k) * d_vec(:, k)';
        Current_phi_y_tilde(k, :, :)= y_tilde_vec(:, k) * y_tilde_vec(:, k)';
    end
    if(n < Summation_Last_Frame)
        Phi_d= Phi_d + Current_phi_d;
        Phi_y_tilde= Phi_y_tilde + Current_phi_y_tilde;
    else 
        if(n==Summation_Last_Frame)
            Phi_d= (Phi_d + Current_phi_d) / Short_Time_Average_Frames;
            Phi_y_tilde= (Phi_y_tilde + Current_phi_y_tilde) / Short_Time_Average_Frames;
        else
            Phi_d= lambda * Phi_d + (1 - lambda) * Current_phi_d;
            Phi_y_tilde= lambda * Phi_y_tilde + (1 - lambda) * Current_phi_y_tilde;
        end
    end
    Phi_u= Phi_d - Phi_y_tilde;
    phi_u= Phi_u(:, 1, 1);
    gamma_u= Phi_u(:, :, 1) ./ phi_u;
    for k=1 : freqs
        phi_u_freq= phi_u(k);
        gamma_u_freq(:, 1)= gamma_u(k, :);
        Phi_u_c(k, :, :)= phi_u_freq * (gamma_u_freq * gamma_u_freq'); 
    end
    Phi_in= Phi_y_tilde + Phi_u - Phi_u_c;
    for k= 1 : freqs 
        Phi_in_freq(:, :)= Phi_in(k, :, :);
        gamma_u_freq(:, 1)= gamma_u(k, :);
        Expression_freq= (Phi_in_freq \ gamma_u_freq); 
        h= Expression_freq / (gamma_u_freq' * Expression_freq);   
        U_f(k, n)= h' * u_vec(:, k);
        Y_re(k, n)= h' * y_vec(:, k);
        U_Hat(k, n)= h' * d_vec(:, k);
        V_re(k, n)= h' * v_vec(:, k);
    end
end
%% Find ISTFTs.
U_f((freqs + 1) : win_length, :)= conj(U_f(((freqs - 1) : -1 : 2), :));
u_f= istft_changed_by_Yuval(U_f, win_length, STFT_jump, 1, beta);
U((freqs + 1) : win_length, :, :)= conj(U(((freqs - 1) : -1 : 2), :, :));
u_reference= istft_changed_by_Yuval(U(:, :, 1), win_length, STFT_jump, 1, beta);
U_Hat((freqs + 1) : win_length, :)= conj(U_Hat(((freqs - 1) : -1 : 2), :));
u_hat= istft_changed_by_Yuval(U_Hat, win_length, STFT_jump, 1, beta);
Y((freqs + 1) : win_length, :, :)= conj(Y(((freqs - 1) : -1 : 2), :, :));
y_reference= istft_changed_by_Yuval(Y(:, :, 1), win_length, STFT_jump, 1, beta);
V((freqs + 1) : win_length, :, :)= conj(V(((freqs - 1) : -1 : 2), :, :));
v_reference= istft_changed_by_Yuval(V(:, :, 1), win_length, STFT_jump, 1, beta);
Y_re((freqs + 1) : win_length, :)= conj(Y_re(((freqs - 1) : -1 : 2), :));
y_re= istft_changed_by_Yuval(Y_re, win_length, STFT_jump, 1, beta);
V_re((freqs + 1) : win_length, :)= conj(V_re(((freqs - 1) : -1 : 2), :));
v_re= istft_changed_by_Yuval(V_re, win_length, STFT_jump, 1, beta);
%% Find ERFs and Distortion indeces.
ERLE= LPF_Function(y_reference.^2)./LPF_Function((y_re).^2);
DI= LPF_Function((u_reference - u_f).^2)./ LPF_Function(u_reference.^2);
NRF= LPF_Function(v_reference.^2)./LPF_Function((v_re).^2);
%% Return Positive frequencies for U_Hat.
U_Hat= U_Hat(1 : freqs, :);
audiowrite('U_Hat_Baseline.wav', u_hat, fs);
end