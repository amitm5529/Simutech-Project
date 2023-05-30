clc, clear
Pt = 511.66;
A1 = 7.97010; B1 = 1521.230; C1 = 233.970;
A2 = 8.01770; B2 = 1715.700; C2 = 234.268;
T1sat = 338; T2sat = 373;

A12 = -87.5644; A21 = 802.2673;
V1 = 42; V2 = 18;
R = 8.314; 
T = zeros();
Psat1 = zeros();
Psat2 = zeros();
a12 = zeros();
a21 = zeros();
m = zeros();
n = zeros();
W = zeros();
x = zeros(); y = zeros(); 
W = zeros();
gamma1 = zeros();
gamma2 = zeros();
T = T2sat:((T1sat-T2sat)/100):T1sat;

 for i = 1:length(T)
     a12(i) = (V2/V1)*exp(-A12./T(i).*R);
     a21(i) = (V1/V2)*exp(-A21./T(i).*R);
    m(i) = -log(x(i)+a12(i).*(1-x(i)));
    n(i)  = -log((1-x(i))+a21(i).*x(i));

     W(i) = a12(i)./(x(i)+(1-x(i)).*a12(i)) - a21(i)./(1+x(i).*(a21(i)-1));
    
    gamma1(i) = exp(m(i)+(1-x(i)).*W(i));
    gamma2(i) = exp(n(i)+x(i).*W(i));
    Psat1(i) = exp(A1-B1/(C1+T(i)));
    Psat2(i) = exp(A2-B2/(C2+T(i)));    
 end

     x = (Pt-gamma2.*Psat2)./(gamma1.*Psat1-gamma2.*Psat2);
    y = x.*gamma1.*Psat1./Pt;
x = 0:0.01:1;
y = 0:0.01:1;

plot(x,T,y,T,'r');
hold on
title('Txy diagram for Methanol/Water Mixture');
xlabel('Mole fraction of Methanol');
ylabel('Temperature');










