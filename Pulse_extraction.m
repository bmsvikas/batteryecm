
clc
clear all
close all


%% Load 206 km pulse data, at BOL (cylce 0)

Working_rep='/Users/joel-louiskone/Documents/Meatech_code/';
Data_rep='/Users/joel-louiskone/Documents/Meatech case study/Data Folder/206k km case/Cycle 0';
mkdir([Working_rep,filesep,'pulses']) %create pulse repository
%% Fixed data

Temp= {'23'};
Crates_ch ={'1','2','4','8'}; % 4 charges C_rates
Crates_dis ={'1','5','10','15'}; % 4 discharges C_rates
Crates_tot={'1','2','4','5','8','10','15'};
load([Working_rep,filesep,'10-16-19_20.16 948_HPPC.mat']);
pulse_duration=10;
Qnom=2.455;% Capacity at cyle 0 (206km aging) 
SOC_init=100;

%% SOC calculation for SOC breakpoints

    Qcumsum=(cumsum(meas.Current.*[diff(meas.Time);0]))/3600;
    SOC=SOC_init+(Qcumsum./Qnom)*100; %SOC from 100% to O
    OCV_SOC=[100;95;91;83;74;66;58;50;42;34;31;28;26;25]; %Calculation during resting phase (at C/3) = begginning of each pulse 
%% Indices for pulses extraction in charge and discharge

meas.Current(abs(meas.Current)<1)=0; %clean current outside HPPC tests

% Initial and final indices for all charge pulses at all C-rates

ind_ch=find(meas.Current>0);
ind_init_ch=[ind_ch(1);ind_ch(find(diff(ind_ch)>50)+1)]; % initial indices of discharge pulses
ind_init_ch=ind_init_ch-1; %indice just before the pulse to extract the OCV
ind_final_ch=[ind_ch((diff(ind_ch)>50));ind_ch(end)];

%Initial and final indices for all discharge pulses at all C-rates
ind_dis=find(meas.Current<0);
ind_init_dis=[ind_dis(1);ind_dis(find(diff(ind_dis)>50)+1)]; % initial indices of discharge pulses
ind_init_dis=ind_init_dis-1; %indice just before the pulse to extract the OCV
ind_final_dis=[ind_dis((diff(ind_dis)>50));ind_dis(end)];


%% Creating and saving data : For loop for each T, total C_rates and SOC


for j= 1:length(Crates_tot) %Crate

    pulses.Temp=str2double(Temp); %Only 1 temperature so no need for  a 'for loop'
    pulses.Crate= str2double(Crates_tot(j));
    savename=strcat ('Meatech',Temp(1),'°', Crates_tot(j),'C', '.mat');

    % We isolate the signals for the selected C-rate
    if sum(contains(Crates_ch,Crates_tot(j)))
        select_Crate_ch=find(strcmp(Crates_ch,Crates_tot(j))); %indice of selected C-rate signal to isolate
        ind_init_ch_select=ind_init_ch(select_Crate_ch:length(Crates_ch):end);
        ind_final_ch_select=ind_final_ch(select_Crate_ch:length(Crates_ch):end);

    end

    %% Charge pulses

    for k= 1:length(OCV_SOC) %SOC measures (from 100% to 10%)

        pulses.('pulsesC').(strcat('measurement',num2str(k)))=[];

        if sum(strcmp(Crates_tot(j),Crates_ch))==1 %if chosen C-rate is not a discharge one, pulses.measurement=[]

            Time_ch=meas.Time(ind_init_ch_select(k):ind_final_ch_select(k));
            Time_ch=Time_ch-Time_ch(1);
            Voltage_ch=meas.Voltage(ind_init_ch_select(k):ind_final_ch_select(k));
            Current_ch=meas.Current(ind_init_ch_select(k):ind_final_ch_select(k));
            Ts_ch=meas.Battery_Temp_degC(ind_init_ch_select(k):ind_final_ch_select(k)); %surface T


            Time1s= find(Time_ch>1,1,'first') ; %1st indice after 1s of pulse (to avoid a potential small peak at 0)
            Time3s=find(Time_ch>3,1,'first');  %1st indice after 3s of pulse

            init_mean=mean(Current_ch(Time1s:Time3s));
            final_mean=mean(Current_ch(Time3s:find(Time_ch<=pulse_duration,1,'last')));

            if diff([abs(init_mean),abs(final_mean)])>-1e-3 % CC verification (eliminate CCCV diff negative)

                pulses.('pulsesC').(strcat('measurement',num2str(k)))=[Time_ch,Voltage_ch,Current_ch,Ts_ch];
            end
        end

    end



    %% Discharge pulses




    for k= 1:length(OCV_SOC) %SOC measures (from 100% to 10% at any C_rates)

        pulses.('pulsesD').(strcat('measurement',num2str(k)))=[];


        if sum(contains(Crates_dis,Crates_tot(j)))
            select_Crate_dis=find(strcmp(Crates_dis,Crates_tot(j))); %indice of selected C-rate signal to isolate
            ind_init_dis_select=ind_init_dis(select_Crate_dis:length(Crates_dis):end);
            ind_final_dis_select=ind_final_dis(select_Crate_dis:length(Crates_dis):end);

        end

        if sum(strcmp(Crates_tot(j),Crates_dis))==1 %if chosen C-rate is not a discharge one, pulses.measurement=[]


            Time_dis=meas.Time(ind_init_dis_select(k):ind_final_dis_select(k));
            Time_dis=Time_dis-Time_dis(1);
            Voltage_dis=meas.Voltage(ind_init_dis_select(k):ind_final_dis_select(k));
            Current_dis=meas.Current(ind_init_dis_select(k):ind_final_dis_select(k));
            Ts_dis=meas.Battery_Temp_degC(ind_init_dis_select(k):ind_final_dis_select(k)); %surface T


            Time1s= find(Time_dis>1,1,'first') ; %1st indice after 1s of pulse (to avoid a potential small peak at 0)
            Time3s=find(Time_dis>3,1,'first');  %1st indice after 3s of pulse

            init_mean=mean(Current_dis(Time1s:Time3s));
            final_mean=mean(Current_dis(Time3s:find(Time_dis<=pulse_duration,1,'last')));

            if diff([abs(init_mean),abs(final_mean)])>-1e-3 % CC verification (eliminate CCCV)

                pulses.('pulsesD').(strcat('measurement',num2str(k)))=[Time_dis,Voltage_dis,Current_dis,Ts_dis];
            end
        end

    end

    %% Saving for selected temperature and chosen Crate_tot

    save ([Working_rep  '/' 'pulses' '/' savename{:}],'pulses')
    clear pulses
end
