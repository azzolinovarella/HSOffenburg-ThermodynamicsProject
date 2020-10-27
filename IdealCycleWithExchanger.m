clear
opengl software
clc

% We know:
p1 = 100;
p4a = 100;
p4b = p4a;
Ro = 8.314;
M = 28.86; 
k = 1.4;
T1 = 20 + 273;
T2b = 644.99;
T4b = T2b;
T3 = 1000 + 273;
PI = 10;



% So, for the p,v-diagram:
R = Ro/M;
p2a = p1*PI;
p2b = p2a;
p3 = p2b;
cp = k*R/(k-1);

v1 = R.*T1./p1;  
v2a = v1./(PI.^(1./k));
v2b = R.*T2b./p2b;
v3 = R.*T3./p3;
v4a = v3.*PI^(1./k);
v4b = R.*T4b./p4b;

T2a = p2a*v2a/R;
T2b = p2b*v2b/R;
T4a = p4a*v4a/R;
T4b = p4b*v4b/R;

v12a = v1:-0.00001:v2a;
v2a2b = v2a:0.001:v2b;
v32b = v3:0.001:v2b;
v34a = v3:0.00001:v4a;
v4a4b = v4a:-0.001:v4b;
v2b4b = (v2b-0.1):0.001:(v4b+0.1);

c12a = (p1.*v1.^k + p2a.*v2a.^k)/2;
c34a = (p3.*v3.^k + p4a.*v4a.^k)/2;

p12a = c12a./(v12a.^k);
p34a = c34a./(v34a.^k);

% Finally
figure(1)
hold on
plot(v12a, p12a, 'color', [0 0 0]); 
line([v2a v2b], [p2a p2b], 'color', [0 0 1])
line([v2b v3], [p2b p3], 'color', [0 1 0])
plot(v34a, p34a, 'color', [0 1 1]); 
line([v4a v4b], [p4a p4b], 'color', [1 0 1])
% line([v4b v1], [p4b p1], 'color', [1 1 0])  Not considering process from 4 --> 1
p2b4b = R.*T2b./v2b4b;
plot(v2b4b, p2b4b, 'k-.')
grid on
ylabel ('p (kPa)')
xlabel ('v (m^3/kg)')
plot(v1,p1,'r*')
text(v1, p1, '1','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(v2a, p2a ,'r*')
text(v2a, p2a, '2a','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(v2b, p2b ,'r*')
text(v2b, p2b, '2b','HorizontalAlignment', 'left','VerticalAlignment','bottom')
plot(v3,p3,'r*')
text(v3, p3, '3','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(v4a, p4a ,'r*')
text(v4a, p4a, '4a','HorizontalAlignment', 'left','VerticalAlignment','top')
plot(v4b, p4b ,'r*')
text(v4b, p4b, '4b','HorizontalAlignment', 'right','VerticalAlignment','bottom')
% legend('1  --> 2a','2a --> 2b', '2b --> 3', '3  --> 4a', '4a --> 4b', '4b --> 1','Isothermal','Location','NorthEast')  Not considering process from 4 --> 1
legend('1  --> 2a','2a --> 2b', '2b --> 3', '3  --> 4a', '4a --> 4b', 'Isothermal','Location','NorthEast')
hold off


% So, for the T,s-diagram:
s1 = 0;
s2a = s1;  % Isentropic!
s2b = cp*log(T2b/T2a) + s2a;
s3 = cp*log(T3/T2b) + s2b;
s4a = cp*log(T4a/T1) + s1;
s4b = cp*log(T4b/T1) + s1;


Ds12a = s1:0.001:s2a;  % Should be zero -> Isentropic!
Ds2a2b = s2a:0.001:s2b;
Ds2b3 = s2b:0.001:s3;
Ds34a = s3:0.001:s4a;  % Should be zero -> Isentropic!
Ds4a4b = s4a:-0.001:s4b;
Ds4b1 = s4b + 0.2:-0.001:s1 - 0.2;  % Not considering process from 4 --> 1

T2a2b = T2a*exp((Ds2a2b - s2a)/cp);
T2b3 = T2b*exp((Ds2b3 - s2b)/cp);
T4a4b = T4a*exp((Ds4a4b - s4a)/cp);
T4b1 = T1*exp((Ds4b1 - s1)/cp);  % Not considering process from 4 --> 1  

% Finaly:
figure(2)
hold on
line([s1 s2a], [T1 T2a], 'color', [0 0 0])
plot(Ds2a2b, T2a2b, 'color', [0 0 1])
plot(Ds2b3, T2b3, 'color', [0 1 0])
line([s3 s4a], [T3 T4a], 'color', [0 1 1])
plot(Ds4a4b, T4a4b, 'color', [1 0 1])
plot(Ds4b1, T4b1, 'color', [0 0 0], 'LineStyle', '-.')  % Not considering process from 4 --> 1
plot(s1,T1,'r*')
text(s1, T1, '1','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(s2a, T2a ,'r*')
text(s2a, T2a, '2a','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(s2b, T2b ,'r*')
text(s2b, T2b, '2b','HorizontalAlignment', 'left','VerticalAlignment','bottom')
plot(s3,T3,'r*')
text(s3, T3, '3','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(s4a, T4a ,'r*')
text(s4a, T4a, '4a','HorizontalAlignment', 'left','VerticalAlignment','top')
plot(s4b, T4b ,'r*')
text(s4b, T4b, '4b','HorizontalAlignment', 'right','VerticalAlignment','bottom')
legend('1  --> 2a','2a --> 2b', '2b --> 3', '3  --> 4a', '4a --> 4b', 'Isobaric','Location','NorthWest')
grid on
hold off


%{



s1 = cp*log(T1) - R*log(p1);
s2a = cp*log(T2a) - R*log(p2a);
s2b = cp*log(T2b) - R*log(p2b);
s3 = cp*log(T3) - R*log(p3);
s4a = cp*log(T4a) - R*log(p4a);
s4b = cp*log(T4b) - R*log(p4b);

Ds12a = 0;  % Isentropic!
Ds2a2b = s2a:0.001:s2b;  
Ds2b3 = s2b:0.001:s3;
Ds34a = 0;  % Isentropic!
Ds4a4b = s4a:-0.001:s4b;
Ds4b1 = s4b:-0.001:s1;

T2a2b = exp((Ds2a2b + R*log(p2a))/cp);  % REVER!!!!!!!!!!!!!!!!!!!!!
T2b3 = exp((Ds2b3+R*log(p2b))/cp);
T4a4b = exp((Ds4a4b + R*log(p4a))/cp);
T4b1 = exp((Ds4b1+R*log(p4b))/cp);

% Finaly:
figure(2)
hold on
line([s1 s2a], [T1 T2a], 'color', [0 0 0])
plot(Ds2a2b, T2a2b, 'color', [0 0 1])
plot(Ds2b3, T2b3, 'color', [0 1 0])
line([s3 s4a], [T3 T4a], 'color', [0 1 1])
plot(Ds4a4b, T4a4b, 'color', [1 0 1])
% plot(Ds4b1, T4b1, 'color', [1 1 0])  Not considering process from 4 --> 1


grid on
ylabel ('T (K)')
xlabel ('s(J/(K*kg))')
plot(s1,T1,'r*')
text(T1, T1, '1','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(s2a, T2a ,'r*')
text(s2a, T2a, '2a','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(s2b, T2b ,'r*')
text(s2b, T2b, '2b','HorizontalAlignment', 'left','VerticalAlignment','bottom')
plot(s3,T3,'r*')
text(s3, T3, '3','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(s4a, T4a ,'r*')
text(s4a, T4a, '4a','HorizontalAlignment', 'left','VerticalAlignment','top')
plot(s4b, T4b ,'r*')
text(s4b, T4b, '4b','HorizontalAlignment', 'right','VerticalAlignment','bottom')
% legend('1  --> 2a','2a --> 2b', '2b --> 3', '3  --> 4a', '4a --> 4b', '4b --> 1','Location','NorthWest')  Not considering process from 4 --> 1  
hold off
%}

% To estimate the amount of work produced we can use numerical integration
w12a = 0;  % This will change!
w2a2b = 0;  % This will continue 0 --> Isobaric!
w2b3 = 0;  % This will continue 0 --> Isobaric! 
w34a = 0;  % This will change!
w4a4b = 0;  % This will continue 0 --> Isobaric!
% w4b1 = 0;  % This will continue 0 --> Isobaric! Not considering process from 4 --> 1
n = 1000;  % Numb of intervals 

p12a = p1 : (p2a-p1)/n : p2a;
p34a = p3 : (p4a-p3)/n : p4a;

v12a = (c12a./p12a).^(1./k);
v34a = (c34a./p34a).^(1./k);

dp = (p2a - p1)/n;
for i=1:length(v12a)
   w12a = w12a + v12a(i).*dp; 
end

dp = (p4a - p3)/n;
for i=1:length(v34a)
   w34a = w34a + v34a(i).*dp; 
end

% wt = w12a + w2a2b + w2b3 + w34a + w4a4b + w4b1;  Not considering process from 4 --> 1
wt = w12a + w2a2b + w2b3 + w34a + w4a4b;
qinv = cp*(T3 - T2b);
efi = -wt/qinv;

disp("The amount of specific net work produced is " +(-1*wt)+ " kJ/kg")
disp("The efficiency is " +efi*100+ "%")