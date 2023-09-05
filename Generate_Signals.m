function [y, u, v, d] = Generate_Signals(M, x_First_Loudspeaker, x_Second_Loudspeaker, s_First_Speaker, s_Second_Speaker, fs, segment_length)
%% This function generates the received signals by the microphones.
c= 340;
Radius= 0.075;
T_60= 0.3;
SNR_dB= 30;
Room_Dims= [6, 6, 4.5];
First_Loudspeaker_Location= [3, 3, 0.1];
Second_Loudspeaker_Location= [3, 3, 0.5];
First_Speaker_Location= [3.5, 3, 0.5];
Second_Speaker_Location= [2.5, 3, 0.5];
theta_vec= (linspace(0, 2 * pi * (M - 1) / M, M))';
for m= 1: M
    Sensor_First_Array= First_Loudspeaker_Location + [Radius * cos(theta_vec(m)), Radius * sin(theta_vec(m)), 0];
    Sensor_Second_Array= Second_Loudspeaker_Location + [Radius * cos(theta_vec(m)), Radius * sin(theta_vec(m)), 0];
    g_First_Array_First_Loudspeaker= rir_generator(c, fs, Sensor_First_Array, First_Loudspeaker_Location, Room_Dims, T_60);
    g_First_Array_Second_Loudspeaker= rir_generator(c, fs, Sensor_First_Array, Second_Loudspeaker_Location, Room_Dims, T_60);
    g_Second_Array_First_Loudspeaker= rir_generator(c, fs, Sensor_Second_Array, First_Loudspeaker_Location, Room_Dims, T_60);
    g_Second_Array_Second_Loudspeaker= rir_generator(c, fs, Sensor_Second_Array, Second_Loudspeaker_Location, Room_Dims, T_60);
    y_First_Array= fftfilt(g_First_Array_First_Loudspeaker, x_First_Loudspeaker) + fftfilt(g_First_Array_Second_Loudspeaker, x_Second_Loudspeaker);
    y_Second_Array= fftfilt(g_Second_Array_First_Loudspeaker, x_First_Loudspeaker) + fftfilt(g_Second_Array_Second_Loudspeaker, x_Second_Loudspeaker);
    y(:, m)= [y_First_Array(1 : (3 * segment_length)) ; y_Second_Array((3 * segment_length + 1) : (5 * segment_length))];
    audiowrite(strcat('Y_m=', num2str(m), '_M=', num2str(M), '.wav'), y(:, m), fs);
    g_First_Array_First_Speaker= rir_generator(c, fs, Sensor_First_Array, First_Speaker_Location, Room_Dims, T_60);
    g_First_Array_Second_Speaker= rir_generator(c, fs, Sensor_First_Array, Second_Speaker_Location, Room_Dims, T_60);
    g_Second_Array_First_Speaker= rir_generator(c, fs, Sensor_Second_Array, First_Speaker_Location, Room_Dims, T_60);
    g_Second_Array_Second_Speaker= rir_generator(c, fs, Sensor_Second_Array, Second_Speaker_Location, Room_Dims, T_60);
    u_First_Array= fftfilt(g_First_Array_First_Speaker, s_First_Speaker) + fftfilt(g_First_Array_Second_Speaker, s_Second_Speaker);
    u_Second_Array= fftfilt(g_Second_Array_First_Speaker, s_First_Speaker) + fftfilt(g_Second_Array_Second_Speaker, s_Second_Speaker);
    u(:, m)= [u_First_Array(1 : (3 * segment_length)) ; u_Second_Array((3 * segment_length + 1) : (5 * segment_length))];
    audiowrite(strcat('U_m=', num2str(m), '_M=', num2str(M), '.wav'), u(:, m), fs);
    clean_d= u(:, m) + y(:, m);
    v(:, m)= rms(clean_d) / db2mag(SNR_dB) * randn(size(clean_d)); 
    d(:, m)= clean_d + v(:, m);
    audiowrite(strcat('D_m=', num2str(m), '_M=', num2str(M), '.wav'), d(:, m), fs); 
end
end