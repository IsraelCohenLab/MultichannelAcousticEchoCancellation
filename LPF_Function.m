function LPF_Function = LPF_Function(sig)
%% Implements a LPF.
filter_length= 4096;
sig_length= length(sig);
extended_sig= [zeros(filter_length - 1, 1); sig];
LPF_Function= zeros(sig_length, 1);
for i= 1 : sig_length
    LPF_Function(i)= sum(extended_sig(i : (i + filter_length - 1)));
end
end