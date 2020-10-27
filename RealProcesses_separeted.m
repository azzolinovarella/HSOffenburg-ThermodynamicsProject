opengl software
clear
clc

% We know:
p1 = 100;
p4 = p1;
Ro = 8.314;
M = 28.86; 
k = 1.4;
T1 = 20 + 273;
T3 = 1000 + 273;
PI = 42;
nc = 0.85;
nt = 0.82;
R = Ro/M;
cp = k*R/(k-1);
cv = cp/k;


% So, for the pv-diagram...
R = Ro/M;
T2_perfect = T1*PI^((k-1)/k);
T4_perfect = T3*PI^((1-k)/k);
T2 = (T2_perfect - T1)/nc + T1;
T4 = nt*(T4_perfect - T3) + T3;

p2 = PI*p1;
p3 = p2 - 100;
p3_perfect = p2;

v1 = R*T1/p1;
v2 = R*T2/p2;
v3 = R*T3/p3;
v4 = R*T4/p4;
v3_perfect = R*T3/p3_perfect;

v12 = v1:-0.00001:v2;
v23 = v2:0.000001:v3;
v34 = v3:0.001:v4;
v41 = v4:-0.0001:v1;

n12 = log(p2/p1)/log(v1/v2);
n23 = log(p3/p2)/log(v2/v3);
n34 = log(p4/p3)/log(v3/v4);
n41 = log(p1/p4)/log(v4/v1);

c12 = (p1*v1^n12 + p2*v2^n12)/2; 
c23 = (p2*v2^n23 + p3*v3^n23)/2;
c34 = (p3*v3^n34 + p4*v4^n34)/2;
c41 = (p4*v4^n41 + p1*v1^n41)/2;

p12 = c12./v12.^n12;
p23 = c23./v23.^n23;
p34 = c34./v34.^n34;

cn12 = cv * ((n12-k)/(n12-1));
cn23 = cv * ((n23-k)/(n23-1));
cn34 = cv * ((n34-k)/(n34-1));
cn41 = cv * ((n41-k)/(n41-1));

s1 = 0;
s2 = cp*log(T2/T1) - R*log(p2/p1) + s1;
s2_perfect = s1;
s3 = cp*log(T3/T1) - R*log(p3/p1) + s1;
s4 = cp*log(T4/T1) - R*log(p4/p1) + s1;
s4_perfect = s3;

Ds12 = s1:0.001:s2;  
Ds23 = s2:0.001:s3;
Ds34 = s3:0.001:s4;  

T12 = T1*exp((Ds12 - s1)/cn12);
T23 = T2*exp((Ds23 - s2)/cn23);
T34 = T3*exp((Ds34 - s3)/cn34);
T23_perfect = T2*exp((Ds23 - s2)/cp);



figure(1)
hold on
Ds = -0.05:0.01:0.15;
p2const = T2*exp((Ds - s2)/cp);
p1const = T1*exp((Ds - s1)/cp);
plot(Ds, p1const, 'r-.')
plot(Ds, p2const, 'g-.')
line([s1 s2_perfect], [T1 T2_perfect],'color', [0 0 1])
plot(Ds12, T12, 'color', [0 0 0])
grid on
ylabel ('T (K)')
xlabel ('s (J/K*kg)')
legend('p1 = const', 'p2 = const', 'Ideal process', 'Real Process','Location','SouthEast')
hold off



figure(2)
hold on
line([v2 v2], [4000 4400], 'color', [0 1 0], 'LineStyle', '-.')
line([v3_perfect v3_perfect], [4000 4400], 'color', [1 0 1], 'LineStyle', '-.')
line([v3 v3], [4000 4400], 'color', [1 0 0], 'LineStyle', '-.')
line([v2 v3_perfect], [p2 p3_perfect], 'color', [0 0 0])
plot(v23, p23, 'color', [0 0 1])
grid on
ylabel ('p (kPa)')
xlabel ('v (m^3/kg)')
legend('v2 = const', 'v3_p_e_r_f_e_c_t = const', 'v3_r_e_a_l_ = const', 'Ideal process', 'Real Process','Location','NorthWest')
hold off


figure(3)
hold on
Ds = 0.35:0.01:0.75;
p3const = T2*exp((Ds - s2)/cp);
p4const = T1*exp((Ds - s1)/cp);
plot(Ds, p3const, 'g-.')
plot(Ds, p4const, 'r-.')
line([s3 s4_perfect], [T3 T4_perfect], 'color', [0 0 0])
plot(Ds34, T34, 'color', [0 0 1])
grid on
ylabel ('T (K)')
xlabel ('s (J/K*kg)')
legend('p3 = const', 'p4 = const', 'Ideal process', 'Real Process','Location','NorthWest')

