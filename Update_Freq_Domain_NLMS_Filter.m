function [h_NLMS, sigma_e_squared, sigma_x_squared] = Update_Freq_Domain_NLMS_Filter(h_NLMS, sigma_e_squared, sigma_x_squared, x_NLMS, sigma_v_squared, D, Q_tag, K)
%% Implements the Frequency Domain NLMS Filter.
epsilon= 1e-8;
constant= 20;
beta= 1 - 1 / (6 * Q_tag);
e= D - h_NLMS' * x_NLMS;
sigma_x_squared= beta * sigma_x_squared + (1 - beta) * abs(x_NLMS(K + 1))^2;
sigma_e_squared= beta * sigma_e_squared + (1 - beta) * abs(e)^2;
miu= 1 / (epsilon + constant * sigma_x_squared + x_NLMS' * x_NLMS) * (1 - sqrt(sigma_v_squared) / (epsilon + sqrt(sigma_e_squared)));
h_NLMS= h_NLMS + miu * x_NLMS * e';
end

