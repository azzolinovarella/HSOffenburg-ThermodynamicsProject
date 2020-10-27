clear
clc

PIi = 1;
PIf = 200;
T3i = 850;
T3f = 2273;
k = 1.4;

n = 50;  % Numero de intervalos!

PI = PIi: (PIf - PIi)/n :PIf;
T3 = T3i: (T3f - T3i)/n :T3f;

[PI, T3] = meshgrid(PI, T3);

w = -1*(297.40 .* (PI.^0.28 - 1) + 1.02 .* T3 .* (PI.^-0.28 - 1));

figure(1)
surf(PI, T3, w)
zlabel ('w_p_r_o_d_u_c_e_d (kJ/kg)')
xlabel ('\Pi')
ylabel('T_3 (K)')
hold on
grid on

PI = 42;
T3 = 850 :1: 2273;
w = -1*(297.40 .* (PI.^0.28 - 1) + 1.02 .* T3 .* (PI.^-0.28 - 1));
figure(2)
plot(T3, w)
grid on
ylabel ('w_p_r_o_d_u_c_e_d (kJ/kg)')
xlabel('T_3 (K)')

T3 = 1273;
PI = 1:200;
w = -1*(297.40 .* (PI.^0.28 - 1) + 1.02 .* T3 .* (PI.^-0.28 - 1));
figure(3)
plot(PI, w)
grid on
ylabel ('w_p_r_o_d_u_c_e_d (kJ/kg)')
xlabel('\Pi')

ef = 1 - PI.^-((k-1)/k);
figure(4)
plot(PI, ef)
grid on
ylabel ('\eta_t_h')
xlabel('\Pi')