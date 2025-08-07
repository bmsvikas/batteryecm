%% This code runs the parameter estimation for the ECM model.

clc
clear all
close all

%% ECM functions repository
Rep_Meatech='/Users/joel-louiskone/Documents/Meatech_code';
addpath('/Users/joel-louiskone/Documents/Meatech_code')

Rep_pulses=[Rep_Meatech,filesep,'pulses']; %pulses repository
mkdir (Rep_pulses)
addpath(Rep_pulses)

Rep_results= [Rep_pulses,filesep,'results'];
mkdir(Rep_results)
addpath(Rep_results)

%%  User inputs
% User must select a number of options

EstimationOpts = struct;
EstimationOpts.nRCpairs           = 2;        % number of RC pairs to be fitted
EstimationOpts.pulse_duration     = 10;


%% OCV values at BOL and Cell capacity (206 km case, cycle 0)
cycle='0';
[OCV_SOCs,OCV_Temps,OCV_values,Cellcapa]= Get_OCV(cycle);

EstimationOpts.OCV_SOC0     = OCV_SOCs;
EstimationOpts.OCV_values0   =OCV_values;
EstimationOpts.Cellcapa     = Cellcapa;
%% Load pulses
cd (Rep_pulses)
files        = dir('*.mat');

%%
% 2nd part : Parameter estimation
% Loop through the mat files containing the pulses structures. For each
% matfile, a folder is created in which the fitting results and figures are
% saved.



for i=1:length(files)

    disp(['Running file',' ',num2str(i),'/',num2str(length(files))])

    for j=1:2%1:2 % Charge/Discharge

        %%

        CurrSign = ["Charge", "Discharge"];
        CurrSign = CurrSign{j};

        Testname = files(i).name;
        load(Testname)

        %%  Step 1 : ◊Run optimization process using HPPC at BOL
        %
        % OCV extraction for each pulses to perform: Voltage-pulse

        [OCV,pulses_signal,substracted] = OCV_pulses_Meatech(Testname,CurrSign,EstimationOpts);


        % ECM parameters optimization : calulate R0 and optimize RC pairs

        if sum(pulses_signal.Nmeas)~=0

            [Results]      = Parameter_Estimation_Multicore_Meatech(EstimationOpts,OCV,pulses_signal,substracted,CurrSign,Testname);


        end
    end


end


%% Step 2 : Generate lookuptable LUT
SOC_breakpoints=0:1:100;
Temp_breakpoints=23:1:26;
cycle='0';
%Generate Lookuptable LUT
[LookUpTables_HPPC]=LUT(Rep_results,SOC_breakpoints,Temp_breakpoints,cycle);

fields = fieldnames(LookUpTables_HPPC);  % 'fields' is now a cell array of strings

% Loop over each field name and assign it as a variable in the current workspace
for i = 1:numel(fields)
    fiedlName = fields{i};  % Use a clear name for readability
    eval([fiedlName ' = LookUpTables_HPPC.' fiedlName ';']);
end

% Save each of the unpacked fields (now variables) individually to the MAT-file
save([Rep_results,filesep, 'LookUpTables_HPPC.mat'], fields{:});
save([Rep_Meatech,filesep, 'LookUpTables_HPPC.mat'], fields{:});

%% Step 3 refine parameters using cycling profile (WLTP)

% The goal is to refine each parameters in SOC-range by filtering the
% signal around a certain temperature, C-rate and SOC tolerance

LookUpTables=load([Rep_Meatech, filesep,'LookUpTables_HPPC.mat']); 

load([Rep_Meatech, filesep,'11-04-19_18.28 1001_WLTP206c.mat']); %WLTP signal, cycle 8(300 km)


% Experimental data (Temp, C-rate, SOC)
SOC_init_WLTP=50;
Capa_300km=2.4309; %206 km, cycle 8(300 km)
current=meas.Current;
Crate=current/Capa_300km; %current vector is converted in C-rate
time=meas.Time;
voltage=meas.Voltage;
Temperature=meas.Battery_Temp_degC;
Qcumsum=(cumsum(current.*[diff(time);0]))/3600;
SOC=SOC_init_WLTP+(Qcumsum./Capa_300km)*100; % SOC in %

% Tolerance and median tempereture and C-rate for WLTP
tolerance_temp=1; %1°C tolerance on the temperature 
tolerance_Crate=0.5;% 50% tolerance on the C-rate
tolerance_soc=2; % SOC is in percentage, ±2 tolerance
Temp_median_WLTP =median(abs(Temperature)); % Assume WLTP test was done near 25°C
Crate_median_WLTP = median(abs(Crate)); 
SOC_breakpoints=LookUpTables.SOCs; %identification in SOC range
Temp_breakpoints=LookUpTables.Temps;
Crate_breakpoints=LookUpTables.Crates;

%Which SOC_breakpoints do we use for refining the parameters? Depend on the
%WLTP range SOC
validity_SOC_breakpoints = SOC_breakpoints(SOC_breakpoints >= min(SOC) & SOC_breakpoints <= max(SOC)); %filter using SOC breakpoints 


% Loop only over SOC breakpoints

% Preallocate result arrays (1D)
num_soc = length(validity_SOC_breakpoints);
Tau1_array = zeros(num_soc, 1);
Tau2_array = zeros(num_soc, 1);
R1_array = zeros(num_soc, 1);
R2_array = zeros(num_soc, 1);
temp_indices = zeros(num_soc, 1);
crate_indices = zeros(num_soc, 1);

    ODE.nRCpairs = 2;
    ODE.y0 = zeros(2,1);

parfor indice_soc = 1:num_soc
    % Filter indices
    indice_filter = find(abs(SOC - validity_SOC_breakpoints(indice_soc)) <= tolerance_soc & ...
                   abs(Temperature - Temp_median_WLTP) <= tolerance_temp & ...
                   abs(abs(Crate) - Crate_median_WLTP) <= tolerance_Crate);
    
  

    % Build signal
    Signal = zeros(numel(indice_filter), 3);
    Signal(:,1) = time(indice_filter);
    Signal(:,2) = voltage(indice_filter);
    Signal(:,3) = current(indice_filter);

    % Remove duplicates
    [~, ia] = unique(Signal(:,1), 'stable');
    Signal = Signal(ia,:);
    
    % Interpolation
    CurrentFit = fit(Signal(:,1), Signal(:,3), 'linearinterp');

    % Initial guess
    x0 = [2, 5, 5e-4, 5e-4];
    lb = [1 3 1e-4 1e-4];
    ub = [3 10 1e-3 1e-3];

    % ODE setup
   
    ODE_local=ODE;
    ODE_local.tspan = Signal(:,1);

    % Optimization
    fmincon_opts = optimoptions('fmincon', 'Display', 'off');
    [OPTOUT, ~, ~, ~] = fmincon(@(k)ECMCost_Meatech(k, Signal, CurrentFit, ODE_local), ...
        x0, [], [], [], [], lb, ub, [], fmincon_opts);

    % Save to local arrays
    Tau1_array(indice_soc) = OPTOUT(1);
    Tau2_array(indice_soc) = OPTOUT(2);
    R1_array(indice_soc) = OPTOUT(3);
    R2_array(indice_soc) = OPTOUT(4);

    temp_indices(indice_soc) = find(abs(Temp_breakpoints - Temp_median_WLTP) <= tolerance_temp, 1);
    crate_indices(indice_soc) = find(abs(Crate_breakpoints - Crate_median_WLTP) <= tolerance_Crate, 1);
end

% Assign to LookUpTables after the parfor
for i = 1:num_soc
    soc_i = i;
    temp_i = temp_indices(i);
    crate_i = crate_indices(i);
    
    LookUpTables.Tau1(soc_i, temp_i, crate_i) = Tau1_array(i);
    LookUpTables.Tau2(soc_i, temp_i, crate_i) = Tau2_array(i);
    LookUpTables.R1(soc_i, temp_i, crate_i) = R1_array(i);
    LookUpTables.R2(soc_i, temp_i, crate_i) = R2_array(i);
end

fields = fieldnames(LookUpTables);  % 'fields' is now a cell array of strings

% Loop over each field name and assign it as a variable in the current workspace
for i = 1:numel(fields)
    fiedlName = fields{i};  % Use a clear name for readability
    eval([fiedlName ' = LookUpTables.' fiedlName ';']);
end

% Save each of the unpacked fields (now variables) individually to the MAT-file
save([Rep_Meatech,filesep, 'LookUpTables.mat'], fields{:});
