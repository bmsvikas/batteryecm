
clc
clear all
close all
Rep_Meatech='/Users/joel-louiskone/Documents/Meatech_code';
addpath('/Users/joel-louiskone/Documents/Meatech_code')


%User input
Qnom=2.4162; %Capacity, cycle 16 (600km)
SOC_init=50;

%% Experimental data
load ([Rep_Meatech,filesep, '11-13-19_05.34 1023_WLTP206c.mat']); %cycle 600


time=meas.Time; %Experimental time data
current=meas.Current; %Experimental current data
voltage=meas.Voltage;
T_init=meas.Battery_Temp_degC(1);

%% Load ECM parameters (LookUpTables)

LookUpTables_HPPC=load([Rep_Meatech,filesep,'LookUpTables.mat']);

%%C1,C2
LookUpTables.C1=LookUpTables.Tau1./LookUpTables.R1;
LookUpTables.C2=LookUpTables.Tau2./LookUpTables.R2;

%% Simulation ECM


Ouput=sim('ECM_Meatech.slx'); % Results of data simulating from simulink


% Time, Current, Voltage, SOC   
Time=Ouput.tout;
Voltage_simu=Ouput.Voltage.signals.values;

% Absolute error and mean Absolute error
Voltage_Fit=fit(Time,Voltage_simu,'pchip');
Voltage_simu_interp=Voltage_Fit(time);
Error=abs(Voltage_simu_interp-voltage) ;
mean_Error=mean(Error);


%% Plots

f=figure;
ax=gca;
ax.LineWidth=1;

subplot(2,1,1)
hold on
plot(time,voltage,'r','Markersize',1,'DisplayName','Experimental')
plot(time,Voltage_simu_interp,'b','DisplayName','Voltage\_JL')
grid on
xlabel('time(s)')
ylabel('voltage(V)')

subplot(2,1,2)
hold on
plot(time,Error,'b','DisplayName','Error')
line([time(1) time(end)] ,[mean_Error mean_Error],'Color','b','LineWidth',2,'DisplayName','MeanError')
grid on
xlabel('time(s)')
ylabel('Error(V)')


%% Save results and plots

mkdir ([Rep_Meatech,filesep, 'validation'])
saveas(f,[Rep_Meatech, filesep 'validation',filesep, 'Voltage_and_Error.fig'])




