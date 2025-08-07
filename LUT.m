%%Script calling Update_LookUpTable function to update the LookupTable with
%%all data in the 'results' folder
function [LookUpTables]=LUT(Rep_results,SOC_breakpoints,Temp_breakpoints,cycle)

%First step: fill LookUptables_init from 0:5:100 as beginning and end of SOC (thus SOC0) can be different
%for pulses at a given temperature and  C-rate

EstimationOpts.nRCpairs            = 2;        % number of RC pairs to be fitted.
EstimationOpts.xSOCs             = 0:5:100;
EstimationOpts.nSOCs              = length(EstimationOpts.xSOCs);


%% Files results
files_results        = dir(Rep_results);
files_results=files_results(contains({files_results.name},'°C'));

%% Initialisation of the LookupTable
LookUpTables_init=struct;
LookUpTables_init.nRCpairs = EstimationOpts.nRCpairs;
LookUpTables_init.SOCs    = EstimationOpts.xSOCs ;


for i= 1:length (files_results)
    load([Rep_results, filesep, files_results(i).name ,filesep, 'Results.mat']);
    LookUpTables_init.Temps(i)  = Results.Temp;
    LookUpTables_init.Crates(i) = Results.Crate;
end


% Remove repeated values and extend current to negative values

LookUpTables_init.Temps    =  unique(LookUpTables_init.Temps);
LookUpTables_init.Crates   =  unique(LookUpTables_init.Crates);
LookUpTables_init.Dims     =  [length(LookUpTables_init.SOCs), length(LookUpTables_init.Temps), length(LookUpTables_init.Crates)];

% Create tables for parameters
LookUpTables_init.R0       = ones(LookUpTables_init.Dims)*nan;

for i=1:EstimationOpts.nRCpairs

    LookUpTables_init.(strcat("Tau",num2str(i))) = ones(LookUpTables_init.Dims)*nan;
    LookUpTables_init.(strcat("R",num2str(i))) = ones(LookUpTables_init.Dims)*nan;

end



%%  Fill initial LUT (LookUpTables_init) with results save in the results structures
%each pulses at different C-rate can have different SOC0 vector, we will
%firstly interpolate on the SOC_breakpoints and then interpolate on
%temperature breakpoints

% SOC interpolation for a given temperature and C-rate

for j= 1:length(files_results) %files_Resuts : one folder by temperature and C-rate (charge or discharge) at different SOCs

    load([Rep_results,filesep, files_results(j).name,filesep, 'Results.mat']);

    %Set LookUpTablesinit for RRC parameters (SOC interpolation and fillmissing  )
    [LookUpTables_init] = Tables_Meatech(Results,LookUpTables_init);

end


%% Second step,SOC and temperature interpolation on breakpoints: Use LookUpTables_init to fill the final LUT (Lookuptables) with desired SOCs &t temperatures :
LookUpTables=struct;
LookUpTables.Temps=Temp_breakpoints; %Desired temperature interpolation
LookUpTables.Crates=LookUpTables_init.Crates;
LookUpTables.SOCs=SOC_breakpoints; % Desired SOC interpolation (performed in the first step)
LookUpTables.nRCpairs=LookUpTables_init.nRCpairs;

LookUpTables.Dims     =[length(LookUpTables.SOCs), length(LookUpTables.Temps), length(LookUpTables.Crates)];

% Initialisation ol LUT
LookUpTables.R0       = ones(LookUpTables.Dims)*nan;

for i=1:EstimationOpts.nRCpairs

    LookUpTables.(strcat("Tau",num2str(i))) = ones(LookUpTables.Dims)*nan;
    LookUpTables.(strcat("R",num2str(i))) = ones(LookUpTables.Dims)*nan;

end

%1D interpolation on SOC and temperature using fit function (in our case)
ind_interp_SOC=find(LookUpTables.SOCs>=LookUpTables_init.SOCs(1) & LookUpTables.SOCs<=LookUpTables_init.SOCs(end));


%2D interpolation on R0
for k=1:length(LookUpTables.Crates )

    %Create meshgrig for LookUpTables_init : fit function only accept
    %column vector

    SOC_vector=LookUpTables.SOCs(ind_interp_SOC); % SOC vector used for interpolation on LookUpTables_init
    LookUpTables.R0(ind_interp_SOC,1,k)=interp1(LookUpTables_init.SOCs,LookUpTables_init.R0(:,k),SOC_vector); %init is 1D vector

    %fill missing values on SOC and Temperature
    LookUpTables.R0(:,:,k)=fillmissing(LookUpTables.R0(:,:,k),'nearest'); %SOC Extrapolation
    LookUpTables.R0(:,:,k)=fillmissing(LookUpTables.R0(:,:,k),'nearest',2); %Temperature Extrapolation

    %2D interpolation for R1,Tau1,R2,Tau2...
    for nbr_RC=1:LookUpTables.nRCpairs

       
        %reshape and evaluate at ind_interp_SOC and ind_interp_T
        LookUpTables.(strcat("Tau",num2str(nbr_RC)))(ind_interp_SOC,1,k)=interp1(LookUpTables_init.SOCs,LookUpTables_init.(strcat("Tau",num2str(nbr_RC)))(:,k),SOC_vector);
        LookUpTables.(strcat("R",num2str(nbr_RC)))(ind_interp_SOC,1,k)=interp1(LookUpTables_init.SOCs,LookUpTables_init.(strcat("R",num2str(nbr_RC)))(:,k),SOC_vector);


        %Extrapolate
        LookUpTables.(strcat("Tau",num2str(nbr_RC)))(:,:,k)=fillmissing(LookUpTables.(strcat("Tau",num2str(nbr_RC)))(:,:,k),'nearest'); %SOC interp
        LookUpTables.(strcat("R",num2str(nbr_RC)))(:,:,k)=fillmissing(LookUpTables.(strcat("R",num2str(nbr_RC)))(:,:,k),'nearest'); %SOC interp
        
        LookUpTables.(strcat("Tau",num2str(nbr_RC)))(:,:,k)=fillmissing(LookUpTables.(strcat("Tau",num2str(nbr_RC)))(:,:,k),'nearest',2); %Temp interp
        LookUpTables.(strcat("R",num2str(nbr_RC)))(:,:,k)=fillmissing(LookUpTables.(strcat("R",num2str(nbr_RC)))(:,:,k),'nearest',2); %Temp interp
       
    end
end

%% OCV LUT

%Step 1 : OCV parameters
%initialisation
LookUpTables.OCV_SOCs=SOC_breakpoints;
LookUpTables.OCV_Temps=Temp_breakpoints;
LookUpTables.OCV=ones(length(LookUpTables.OCV_SOCs),length(LookUpTables.OCV_Temps))*NaN;

%Get OCV values for given cycle
[OCV_SOCs,~,OCV_values,~]= Get_OCV(cycle);

ind_SOC_interpOCV=find(LookUpTables.OCV_SOCs>=OCV_SOCs(1) & LookUpTables.OCV_SOCs<=OCV_SOCs(end));

%interpolate
OCV_Fit=fit(OCV_SOCs',OCV_values','pchip');
LookUpTables.OCV(ind_SOC_interpOCV,1)=OCV_Fit(LookUpTables.OCV_SOCs(ind_SOC_interpOCV));

%Fill missing values
LookUpTables.OCV(:,:)=fillmissing(LookUpTables.OCV(:,:),'nearest'); %SOC Extrapolation
LookUpTables.OCV(:,:)=fillmissing(LookUpTables.OCV(:,:),'nearest',2); %Temperature Extrapolation

end