clc
clear all
close all

Rep_Meatech='/Users/joel-louiskone/Documents/Meatech_code';
addpath('/Users/joel-louiskone/Documents/Meatech_code')
mkdir([Rep_Meatech, filesep, 'EKF'])

%% User input
cycle='0';

%% Load experimental data and ECM parameters

% Experimental data
load ('UDDS_25degC.mat')
time_meas=meas.Time;
current_meas=meas.Current;
V_meas=meas.Voltage;
Temp_meas=meas.Battery_Temp_degC;


% Load OCV

[OCV_SOCs,~,OCV_values,Qnom]= Get_OCV(cycle);
OCV.SOCs=OCV_SOCs;
OCV.values=OCV_values;
SOC_init=100;

%Load parameters
LookUpTables=load('LookUpTables.mat'); %R0, R, Tau parameters


%% Coulomb counting

 Qcumsum=(cumsum(current_meas.*[diff(time_meas);0]))/3600;
 SOC_cc=SOC_init+(Qcumsum./Qnom)*100; %SOC from 100% to O

%% EKF 

[SOC_EKF,~,~] = SOC_EKF(current_meas, V_meas, time_meas,Temp_meas, SOC_init, Qnom, LookUpTables, OCV);


f2=figure;
hold on
plot(meas.Time,SOC_cc,'DisplayName', 'SOC coulomb counting')
plot(meas.Time,SOC_EKF,'DisplayName', 'SOC EKF')
xlabel('time(s)')
ylabel('SOC(%)')
grid on

saveas(f2, [Rep_Meatech, filesep, 'EKF', filesep, 'SOCcompare.fig'])