clear
opengl software
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

v1 = R*T1/p1;
v2 = R*T2/p2;
v3 = R*T3/p3;
v4 = R*T4/p4;

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
% p41 = c41./v41.^n41;  Not considering process from 4 --> 1

figure(1)
hold on
grid on
xlabel('v (m^3/kg)')
ylabel('p (kPa)')
plot(v12, p12, 'color', [0 0 0])
plot(v23, p23, 'color', [0 0 1])
plot(v34, p34, 'color', [0 1 0])
% plot(v41, p41, 'color', [1 0 1])  Not considering process from 4 --> 1
plot(v1,p1,'r*')
text(v1, p1, '1','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(v2, p2 ,'r*')
text(v2, p2, '2','HorizontalAlignment', 'right','VerticalAlignment','bottom')
plot(v3,p3,'r*')
text(v3, p3, '3','HorizontalAlignment', 'left','VerticalAlignment','bottom')
plot(v4, p4 ,'r*')
text(v4, p4, '4','HorizontalAlignment', 'left','VerticalAlignment','top')
% legend('1->2', '2->3', '3->4', '4->1', 'Location','NorthEast')  Not considering process from 4 --> 1
legend('1->2', '2->3', '3->4', 'Location','NorthEast')


% So, for the T,s-diagram...
cn12 = cv * ((n12-k)/(n12-1));
cn23 = cv * ((n23-k)/(n23-1));
cn34 = cv * ((n34-k)/(n34-1));
% cn41 = cv * ((n41-k)/(n41-1));  Not considering process from 4 --> 1

% So, for the T,s-diagram:

s1 = 0;
s2 = cp*log(T2/T1) - R*log(p2/p1) + s1;  
s3 = cp*log(T3/T1) - R*log(p3/p1) + s1;
s4 = cp*log(T4/T1) - R*log(p4/p1) + s1; 

Ds12 = s1:0.001:s2;  
Ds23 = s2:0.001:s3;
Ds34 = s3:0.001:s4;  
% Ds41 = s4:-0.001:s1;  Not considering process from 4 --> 1

T12 = T1*exp((Ds12 - s1)/cn12);
T23 = T2*exp((Ds23 - s2)/cn23);
T34 = T3*exp((Ds34 - s3)/cn34);
% T41 = T4*exp((Ds41 - s4)/cn41);  Not considering process from 4 --> 1

% Finaly:
figure(2)
hold on
plot(Ds12, T12, 'color', [0 0 0])
plot(Ds23, T23, 'color', [0 0 1])
plot(Ds34, T34, 'color', [0 1 0])
% plot(Ds41, T41, 'color', [1 0 1])  Not considering process from 4 --> 1
plot(s1,T1,'r*')
text(s1, T1, '1','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(s2, T2 ,'r*')
text(s2, T2, '2','HorizontalAlignment', 'right','VerticalAlignment','bottom')
plot(s3,T3,'r*')
text(s3, T3, '3','HorizontalAlignment', 'left','VerticalAlignment','bottom')
plot(s4, T4 ,'r*')
text(s4, T4, '4','HorizontalAlignment', 'left','VerticalAlignment','top')
grid on
ylabel ('T (K)')
xlabel ('s(J/(K*kg))')
% legend('1->2','2->3', '3->4', '4->1', 'Location','NorthWest')  Not considering process from 4 --> 1
legend('1->2','2->3', '3->4', 'Location','NorthWest')
hold off

disp(cn12)
disp(cn23)
disp(cn34)




% To estimate the amount of work produced we can use numerical integration...
w12 = 0;  % This will change!
w23 = 0;  % This will continue 0 --> Isobaric! 
w34 = 0;  % This will change!
w41 = 0;  % This will continue 0 --> Isobaric!

n = 1000;  % Numb of intervals 

p12 = p1 : (p2-p1)/n : p2;
p23 = p2 : (p3-p2)/n : p3;
p34 = p3 : (p4-p3)/n : p4;
p41 = p4 : (p1-p4)/n : p1;

v12 = (c12./p12).^(1./n12);
v23 = (c23./p23).^(1./n23);
v34 = (c34./p34).^(1./n34);

dp = (p2 - p1)/n;
for i=1:length(v12)
   w12 = w12 + v12(i).*dp; 
end

dp = (p3 - p2)/n;
for i=1:length(v23)
   w23 = w23 + v23(i).*dp; 
end

dp = (p4 - p3)/n;
for i=1:length(v34)
   w34 = w34 + v34(i).*dp; 
end

wt = w12 + w23 + w34 + w41;
qinv = cn12*(T2-T1) + cn23*(T3-T2) + cn34*(T4-T3);
efi = -wt/qinv;

disp("The amount of specific net work produced in the real cycle is " +(-1*wt)+ " kJ/kg")
disp("The efficiency is " +efi*100+ "%")

