    %HFM - Nirjhar - 2.4.20
%Load plot data files
clear all;

global outflow;
load T773.mat;load T798.mat; load T823.mat; load T848.mat;

%Create inflow value from data files

F_CH4_773 = 0.4./T773(:,1); F_CH4_798 = 0.4./T798(:,1); F_CH4_823 = 0.4./T823(:,1); F_CH4_848 = 0.4./T848(:,1);

%Converting to kmol/hr. - 14.4.20
%F_CH4_773 = F_CH4_773/3600; F_CH4_798 = F_CH4_798/3600; F_CH4_823 = F_CH4_823/3600; F_CH4_848 = F_CH4_848/3600; %m mol/hr to mol/s
 
%conversion  from mol/(g.cat hr) to mol/(kg.cat hr)
F_CH4_773 = F_CH4_773*1000; F_CH4_798 = F_CH4_798*1000; F_CH4_823 = F_CH4_823*1000; F_CH4_848 = F_CH4_848*1000; %m mol/hr to mol/s




F_H2O_773 = F_CH4_773*3; F_H2O_798 = F_CH4_798*3; 
F_H2O_823 = F_CH4_823*3; F_H2O_848 = F_CH4_848*3;

F_H2_773 = F_CH4_773*1.25; F_H2_798 = F_CH4_798*1.25;
F_H2_823 = F_CH4_823*1.25; F_H2_848 = F_CH4_848*1.25;

R = 8.314; %J/(mol K)
V = 0.5; %m^3 - unsure - set in dydt now, remove rendundancy

P_CH4_0 = 2; P_H2O_0 = 2; P_CO2_0 = 2; P_CO_0 = 2; P_H2_0 = 2; % bar - unsure


%% T = 773
CH4_conv_773 = zeros(size(F_CH4_773,1),1);
for i = 1:size(F_CH4_773,1)
T = 773; %K
tspan = [0 5];
y0 = [2 2 2 2 2];
F_CH4 = F_CH4_773(i);
options = odeset('RelTol',1e-7, 'AbsTol',1e-9);
fprintf("Methane inflow: %f\n",F_CH4);
[t,y] = ode15s(@(t,y) dydt(t,y,R,T,F_CH4,V),tspan,y0,options);

% Steady state checker

check_start = 3;
for j = dsearchn(t,check_start):size(t,1)
    steady_state = false;
    difference = sum(y(j,:) - y(j-1,:));
    if (difference < 1e-8)
        steady_state = true;
    end
end

%Find conversion if steady state reached, otherwise output warning
if (steady_state) 
    %F_CH4 = 0.001;
    %Getting outflow 
    %[dydt,outflow] = dydt(t,y,R);
    CH4_conv_773(i) = (F_CH4 - ((y(size(y,1),1)/sum(y(size(y,1),:)))*outflow))/F_CH4;
    fprintf("Methane conversion: %f\n",CH4_conv_773(i));
else
    fprintf("Steady state not reached: check run time or input configuration:");
end

end

%Plot and compare
figure;
plot(T773(:,1),T773(:,2));
hold on;
plot(T773(:,1),CH4_conv_773);
title('T-773');

%% T = 848
CH4_conv_848 = zeros(size(F_CH4_848,1),1);
for i = 1:size(F_CH4_848,1)
T = 848; %K
tspan = [0 5];
y0 = [2 2 2 2 2];
F_CH4 = F_CH4_848(i);

fprintf("Methane inflow: %f\n",F_CH4);
[t,y] = ode15s(@(t,y) dydt(t,y,R,T,F_CH4,V),tspan,y0);

% Steady state checker

check_start = 3;
for j = dsearchn(t,check_start):size(t,1)
    steady_state = false;
    difference = sum(y(j,:) - y(j-1,:));
    if (difference < 1e-8)
        steady_state = true;
    end
end

%Find conversion if steady state reached, otherwise output warning
if (steady_state) 
    %F_CH4 = 0.001;
    %Getting outflow 
    %[dydt,outflow] = dydt(t,y,R);
    CH4_conv_848(i) = (F_CH4 - ((y(size(y,1),1)/sum(y(size(y,1),:)))*outflow))/F_CH4;
    fprintf("Methane conversion: %f\n",CH4_conv_848(i));
else
    fprintf("Steady state not reached: check run time or input configuration:");
end

end

%Plot and compare
figure;
plot(T848(:,1),T848(:,2));
hold on;
plot(T848(:,1),CH4_conv_848);
title('T-848');


%% T = 823
CH4_conv_823 = zeros(size(F_CH4_823,1),1);
for i = 1:size(F_CH4_823,1)
T = 823; %K
tspan = [0 5];
y0 = [2 2 2 2 2];
F_CH4 = F_CH4_823(i);

fprintf("Methane inflow: %f\n",F_CH4);
[t,y] = ode15s(@(t,y) dydt(t,y,R,T,F_CH4,V),tspan,y0);

% Steady state checker

check_start = 3;
for j = dsearchn(t,check_start):size(t,1)
    steady_state = false;
    difference = sum(y(j,:) - y(j-1,:));
    if (difference < 1e-8)
        steady_state = true;
    end
end

%Find conversion if steady state reached, otherwise output warning
if (steady_state) 
    %F_CH4 = 0.001;
    %Getting outflow 
    %[dydt,outflow] = dydt(t,y,R);
    CH4_conv_823(i) = (F_CH4 - ((y(size(y,1),1)/sum(y(size(y,1),:)))*outflow))/F_CH4;
    fprintf("Methane conversion: %f\n",CH4_conv_823(i));
else
    fprintf("Steady state not reached: check run time or input configuration:");
end

end

%Plot and compare
figure;
plot(T823(:,1),T823(:,2));
hold on;
plot(T823(:,1),CH4_conv_823);
title('T-823');


%% T = 798
CH4_conv_798 = zeros(size(F_CH4_798,1),1);
for i = 1:size(F_CH4_798,1)
T = 798; %K
tspan = [0 5];
y0 = [2 2 2 2 2];
F_CH4 = F_CH4_798(i);

fprintf("Methane inflow: %f\n",F_CH4);
[t,y] = ode15s(@(t,y) dydt(t,y,R,T,F_CH4,V),tspan,y0);

% Steady state checker

check_start = 3;
for j = dsearchn(t,check_start):size(t,1)
    steady_state = false;
    difference = sum(y(j,:) - y(j-1,:));
    if (difference < 1e-8)
        steady_state = true;
    end
end

%Find conversion if steady state reached, otherwise output warning
if (steady_state) 
    %F_CH4 = 0.001;
    %Getting outflow 
    %[dydt,outflow] = dydt(t,y,R);
    CH4_conv_798(i) = (F_CH4 - ((y(size(y,1),1)/sum(y(size(y,1),:)))*outflow))/F_CH4;
    fprintf("Methane conversion: %f\n",CH4_conv_798(i));
else
    fprintf("Steady state not reached: check run time or input configuration:");
end

end

%Plot and compare
figure;
plot(T798(:,1),T798(:,2));
hold on;
plot(T798(:,1),CH4_conv_798);
title('T-798');

    
    


%for i = 1:size(F_CH4_773,1)
    
