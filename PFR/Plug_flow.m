%Nirjhar - 5/4/20
%Plug_flow
clear all;
close;
%Load plot data files
global outflow;
load T773.mat;load T798.mat; load T823.mat; load T848.mat;

%Create inflow value from data files
F_CH4_773 = 0.4./T773(:,1); F_CH4_798 = 0.4./T798(:,1); F_CH4_823 = 0.4./T823(:,1); F_CH4_848 = 0.4./T848(:,1);

%Unit conversion from mol/(g.cat hr) to mol/(kg.cat hr)
F_CH4_773 = F_CH4_773*1000; F_CH4_798 = F_CH4_798*1000; F_CH4_823 = F_CH4_823*1000; F_CH4_848 = F_CH4_848*1000; %m mol/hr to mol/s


%This is redundnant, already done in dydt.m
F_H2O_773 = F_CH4_773*3; F_H2O_798 = F_CH4_798*3; 
F_H2O_823 = F_CH4_823*3; F_H2O_848 = F_CH4_848*3;

F_H2_773 = F_CH4_773*1.25; F_H2_798 = F_CH4_798*1.25;
F_H2_823 = F_CH4_823*1.25; F_H2_848 = F_CH4_848*1.25;

R = 8.314; %J/kgK
V = 0.5; %m^3 - unsure - set in dydt now, remove rendundancy


%Picking an inflow
F_CH4 = F_CH4_773(1);
F_CH4_0 = F_CH4; %Not sure about the conversion definition
%Setting inflow values for first plug
F_H2O = F_CH4*3; F_H2 =F_CH4*1.25; F_CO2 = 0; F_CO = 0;
%Plugs
no_plugs = 70;

%PFR for T=773K and inflow at F_CH4_773(4)
CH4_conv = zeros(no_plugs,1); %Store steady state conversion of each plug
CO2_conv = zeros(no_plugs,1);
pressure = zeros(no_plugs,5); %Store steady state pressure of each plug
for i = 1:no_plugs
T = 848;%848; %K
tspan = [0 5];
y0 = [2 2 2 2 2];

if (i>1)
F_CH4 = (pressure(i-1,1)/sum(pressure(i-1,:)))*outflow;
F_H2O = (pressure(i-1,2)/sum(pressure(i-1,:)))*outflow;
F_CO2 = (pressure(i-1,3)/sum(pressure(i-1,:)))*outflow;
F_H2 = (pressure(i-1,4)/sum(pressure(i-1,:)))*outflow;
F_CO = (pressure(i-1,5)/sum(pressure(i-1,:)))*outflow;
end

fprintf("Methane inflow: %f\n",F_CH4);
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode15s(@(t,y) dydt(t,y,R,T,F_CH4,F_H2O,F_CO2,F_H2,F_CO,V),tspan,y0,options);

% Steady state checker
check_start = 3;
for j = dsearchn(t,check_start):size(t,1)
    steady_state = false;
    difference = abs(sum(y(j,:) - y(j-1,:)));
    if (difference < 1e-8)
        steady_state = true;
    end
end

%Find conversion if steady state reached, otherwise output warning
if (steady_state) 
    %F_CH4 = 0.001;
    %Getting outflow 
    %[dydt,outflow] = dydt(t,y,R);
    CH4_conv(i) = (F_CH4_0 - ((y(size(y,1),1)/sum(y(size(y,1),:)))*outflow))/F_CH4_0; %(F_CH4 - ((y(size(y,1),1)/sum(y(size(y,1),:)))*outflow))/F_CH4;
    CO2_conv(i) = ((y(size(y,1),3)/sum(y(size(y,1),:)))*outflow)/F_CH4_0;
    fprintf("Methane conversion: %f\n",CH4_conv(i));
    pressure(i,:) = y(size(y,1),:);
else
    fprintf("Steady state not reached: check run time or input configuration:");
end

end


%Plot
figure;
plot(linspace(1,no_plugs,no_plugs),CO2_conv,'--');
xlabel('No of plugs');ylabel('CO2 conversion');
title('CO2 conversion');

figure;
plot(linspace(1,no_plugs,no_plugs),CH4_conv,'--');
xlabel('No of plugs');ylabel('CH4 conversion');
title('CH4 conversion');

figure;
plot(linspace(1,no_plugs,no_plugs),pressure);
xlabel('No of plugs'); ylabel('Pressures(bar)');
title('Pressure of each species');
