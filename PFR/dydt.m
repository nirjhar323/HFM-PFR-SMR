function dydt = dydt(t,y,R,T,F_CH4_0,F_H2O_0,F_CO2_0,F_H2_0,F_CO_0,V)  %F_CH4_0 passed as F_CH4 in HFM.m
global outflow;
%R = 8.314; %J/kgK
%V = 0.5; %m^3 - unsure

%P_CH4_0 = 2; P_H2O_0 = 2; P_CO2_0 = 2; P_CO_0 = 2; P_H2_0 = 2; % bar - unsure

%F_CH4_0 = 0.001; 
%F_H2O_0 = 3*F_CH4_0; F_H2_0 = 1.25*F_CH4_0;
%F_CO2_0 = 0; F_CO_0 = 0;

K1=1.448e13*exp(-221901/8.314/T);
K2=2.151e-2*exp(35030/8.314/T);
K3=3.116e11*exp(-189360/8.314/T);

KCH4=6.65e-4*exp(38280/8.314/T);
KCO=8.23e-5*exp(70650/8.314/T);
KH2=6.12e-9*exp(82900/8.314/T);
KH2O=1.77e5*exp(-88680/8.314/T);

kin1=4.225e15*exp(-240100/8.314/T);
kin2=1.955e6*exp(-67130/8.314/T);
kin3=1.020e15*exp(-243900/8.314/T);

% DEN=1+KCO*P_CO_0+KH2*P_H2_0+KCH4*P_CH4_0+KH2O*P_H2O_0/P_H2_0;
% r1=((kin1/P_H2_0^2.5)/DEN^2)*(P_CH4_0*P_H2O_0-P_H2_0^3*P_CO_0/K1);
% r2=((kin2/P_H2_0)/DEN^2)*(P_CO_0*P_H2O_0-P_H2_0*P_CO2_0/K2);
% r3=((kin3/P_H2_0^3.5)/DEN^2)*(P_CH4_0*P_H2O_0^2-P_H2_0^4*P_CO2_0/K3);

% Solve ODE system - y(1) = P_CH4, y(2) = P_H2O, y(3) = P_CO2, y(4) = P_H2,
% y(5) = P_CO

P_total = sum(y);

DEN=1+KCO*y(5)+KH2*y(4)+KCH4*y(1)+KH2O*y(2)/y(4);
r1=((kin1/y(4)^2.5)/DEN^2)*(y(1)*y(2)-y(4)^3*y(5)/K1);
r2=((kin2/y(4))/DEN^2)*(y(5)*y(2)-y(4)*y(3)/K2);
r3=((kin3/y(4)^3.5)/DEN^2)*(y(1)*y(2)^2-y(4)^4*y(3)/K3);

%Conversion from kmol/kgcat.hr to mol/kgcat.hr - 2_5_20
r1 = r1*1000; r2 = r2*1000; r3 = r3*1000;


F_in = F_CH4_0 + F_H2O_0 + F_H2_0 + F_CO_0 + F_CO2_0;

%Conversion factor - Pa to bar - mol/hr to mol/s - 2_5_20
factor = 1/(101325);%3600/101325

dydt(1) = (R*T/V)*factor*(F_CH4_0 - ((r1 + r3)*V) - (y(1)/P_total)*(F_in + 2*V*r1 + 2*V*r3));
dydt(2) = (R*T/V)*factor*(F_H2O_0 - ((r1 + r2 + 2*r3)*V) -(y(2)/P_total)*(F_in + 2*V*r1 + 2*V*r3));
dydt(3) = (R*T/V)*factor*(F_CO2_0 + ((r2 + r3)*V)- (y(3)/P_total)*(F_in + 2*V*r1 + 2*V*r3));
dydt(4) = (R*T/V)*factor*(F_H2_0 + ((3*r1 + r2 + 4*r3)*V)- (y(4)/P_total)*(F_in + 2*V*r1 + 2*V*r3));
dydt(5) = (R*T/V)*factor*(F_CO_0 + ((r1 - r2)*V)- (y(5)/P_total)*(F_in + 2*V*r1 + 2*V*r3));

dydt = dydt';
outflow = (F_in + 2*V*r1 + 2*V*r3);
end