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
R = Ro/M;
p2 = p1*PI;
p3 = p2;
v1 = R.*T1./p1;  
v2 = v1./(PI.^(1./k));
v3 = R.*T3./p3;
v4 = v3.*PI^(1./k);
cp = k*R/(k-1);
T2 = p2*v2/R;
T4 = p4*v4/R;



% So, for the p,v-diagram...

v12 = v1:-0.00001:v2;
v34 = v3:0.00001:v4;

c12 = (p1.*v1.^k + p2.*v2.^k)/2;
c34 = (p3.*v3.^k + p4.*v4.^k)/2;

p12 = c12./(v12.^k);
p34 = c34./(v34.^k);

% Finaly:
figure(1)
hold on
plot(v12, p12, 'color', [0 0 0]);
line([v2 v3], [p2 p3], 'color', [0 0 1])
plot(v34, p34, 'color', [0 1 0])
% line([v4 v1], [p4 p1], 'color', [1 0 1])  Not considering process from 4 --> 1
plot(v1,p1,'r*')
text(v1, p1, '1','HorizontalAlignment', 'right','VerticalAlignment','top')
plot(v2, p2 ,'r*')
text(v2, p2, '2','HorizontalAlignment', 'right','VerticalAlignment','bottom')
plot(v3,p3,'r*')
text(v3, p3, '3','HorizontalAlignment', 'left','VerticalAlignment','bottom')
plot(v4, p4 ,'r*')
text(v4, p4, '4','HorizontalAlignment', 'left','VerticalAlignment','top')
grid on
ylabel ('p (kPa)')
xlabel ('v (m^3/kg)')
legend('1->2','2->3', '3->4','Location','NorthEast') 
% legend('1->2','2->3', '3->4', '4->1','Location','NorthEast')   Not considering process from 4 --> 1
hold off



% So, for the T,s-diagram:

s1 = 0;
s2 = s1;  % Isentropic!
s3 = cp*log(T3/T2) + s2;
s4 = cp*log(T4/T1) + s1;

Ds12 = s1:0.001:s2;  % Should be zero -> Isentropic!
Ds23 = s2:0.001:s3;
Ds34 = s3:0.001:s4;  % Should be zero -> Isentropic!
% Ds41 = s4:-0.001:s1;  Not considering process from 4 --> 1

T23 = T2*exp((Ds23 - s2)/cp);
% T41 = T1*exp((Ds41 - s1)/cp);  Not considering process from 4 --> 1  

% Finaly:
figure(2)
hold on
line([s1 s2], [T1 T2], 'color', [0 0 0])
plot(Ds23, T23, 'color', [0 0 1])
line([s3 s4], [T3 T4], 'color', [0 1 0])
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
legend('1->2','2->3', '3->4', 'Location','NorthWest')
% legend('1->2','2->3', '3->4', '4->1','Location','NorthWest')  Not considering process from 4 --> 1
hold off



% To estimate the amount of work produced we can use numerical integration...
w12 = 0;  % This will change!
w23 = 0;  % This will continue 0 --> Isobaric! 
w34 = 0;  % This will change!
% w41 = 0;  % This will continue 0 --> Isobaric!  Not considering process from 4 --> 1  

n = 1000;  % Numb of intervals 

p12 = p1 : (p2-p1)/n : p2;
p23 = p2 : (p3-p2)/n : p3;
p34 = p3 : (p4-p3)/n : p4;
% p41 = p4 : (p1-p4)/n : p1;  Not considering process from 4 --> 1

v12 = (c12./p12).^(1./k);
v34 = (c34./p34).^(1./k);

dp = (p2 - p1)/n;
for i=1:length(v12)
   w12 = w12 + v12(i).*dp; 
end

dp = (p4 - p3)/n;
for i=1:length(v34)
   w34 = w34 + v34(i).*dp; 
end

% wt = w12 + w23 + w34 + w41  Not considering process from 4 --> 1
wt = w12 + w23 + w34
efi = 1 - PI^(-(k-1)/k);

disp("The amount of specific net work produced is " +(-1*wt)+ " kJ/kg")
disp("The efficiency is " +efi*100+ "%")
