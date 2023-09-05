function [ERLE, DI, NRF, U_Hat, u_hat] = Spatial_LCMV_Method(fs, frames, win_length, STFT_jump, beta, M, L, X, D, U, Y, V)
%% This function implements the spatial MVDR filter minus the corresponding NLMS output.
%% Define Parameters.
freqs= win_length / 2 + 1;
%% Implements the total MVDR beamformer.
U_f= zeros(freqs, frames);
Y_re= zeros(freqs, frames);
V_re= zeros(freqs, frames);
U_Hat= zeros(freqs, frames);
Y_tilde= zeros(freqs, frames, M);
ms= nchoosek(1 : M, 2);
ls= nchoosek(1 : L, 2);
m_couples= size(ms, 1);
l_couples= size(ls, 1);
Equations= m_couples * l_couples;
res_vec= zeros(freqs, Equations);
Mat_g= zeros(freqs, Equations, M);
Current_res_vec= zeros(freqs, Equations);
Current_Mat_g= zeros(freqs, Equations, M);
g_vec= zeros(M, freqs);
Current_u_tilde_expectation_vec= zeros(freqs, M);
for n= L : frames
    for k= 1 : freqs
        for o= 1 : m_couples
            m1= ms(o, 1);
            m2= ms(o, 2);
            for p= 1 : l_couples
                l1= ls(p, 1);
                l2= ls(p, 2);
                equation= (o - 1) * l_couples + p;
                Current_res_vec(k, equation)= D(k, n - l1 + 1, m1) * D(k, n - l2 + 1, m2) - D(k, n - l2 + 1, m1) * D(k, n - l1 + 1, m2);
                Current_Mat_g(k, equation, m1)= X(k, n - l1 + 1) * D(k, n - l2 + 1, m2) - X(k, n - l2 + 1) * D(k, n - l1 + 1, m2);
                Current_Mat_g(k, equation, m2)= X(k, n - l2 + 1) * D(k, n - l1 + 1, m1) - X(k, n - l1 + 1) * D(k, n - l2 + 1, m1);
            end
        end
    end
    res_vec(:, :)= Current_res_vec;
    Mat_g(:, :, :)=Current_Mat_g(:, :, :);
    d_vec= transpose(reshape(D(:, n, :), freqs, M));
    y_vec= transpose(reshape(Y(:, n, :), freqs, M));
    u_vec= transpose(reshape(U(:, n, :), freqs, M));
    v_vec= transpose(reshape(V(:, n, :), freqs, M));
    for k= 1 : freqs
        Mat_g_freq(:, :)= Mat_g(k, :, :);
        res_vec_freq(:, 1)= res_vec(k, :);
        g_vec(:, k)= pinv(Mat_g_freq) * res_vec_freq;
        Y_tilde(k, n, :)= g_vec(:, k) * X(k, n);
        Current_u_tilde_expectation_vec(k, :)= D(k, n, :) - Y_tilde(k, n, :);
    end
    U_tilde_expectation_vec= Current_u_tilde_expectation_vec;
    q_vec= transpose(U_tilde_expectation_vec ./ U_tilde_expectation_vec(:, 1));
    for k=1 : freqs
        C= [q_vec(:, k), g_vec(:, k)];
        h= C * ((C' * C) \ [1 ; 0]);
        U_f(k, n)= h' * u_vec(:, k);
        Y_re(k, n)= h' * y_vec(:, k);
        V_re(k, n)= h' * v_vec(:, k);
        U_Hat(k, n)= h' * d_vec(:, k);
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
%% Take Positive Frequencies of U_Hat.
U_Hat= U_Hat(1 : freqs, :);
audiowrite(strcat('U_Hat_M=', num2str(M), '_L=', num2str(L), '.wav'), u_hat, fs);
end