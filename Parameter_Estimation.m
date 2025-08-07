function      [Results]      = Parameter_Estimation_Multicore_Meatech(EstimationOpts,OCV,pulses_signal,substracted,CurrSign,Testname)

% Fitting tool for the ECM characterization pulses.
% Pulse signals are solved using the ODE functions of the RC pairs.
% A local optimization function is used to find the local optimas for the
% fitting while a global optimization tries to find the global optimum.
% More in details, the GO tool passes the initial starting points to the
% local function and subesequently the later communicates the value of the
% objective function.

loaded=load(Testname);
pulses=loaded.pulses;


% Determine folder name in which results, such as fit figures, will be
% saved
if CurrSign == "Charge"
    Currname = "C";
    Pulsename  = "pulsesC";
    sign_CurrSign=1;
elseif CurrSign =="Discharge"
    Currname = "D";
    Pulsename  = "pulsesD";
    sign_CurrSign=-1;
end

% foldername format Crate - C/D - Temp -C = 6C25C
foldername =  strcat(num2str(pulses.Crate),Currname,num2str(pulses.Temp),"°C");

% Create folder
mkdir(strcat (pwd, filesep, 'results', filesep,foldername))
%1 - Get the data from the measurement
Nmeas             = pulses_signal.Nmeas; %Number of pulses measurements
Liste_meas        = pulses_signal.Liste_meas;
ind_meas          = pulses_signal.ind_meas;

% Remove structures to run parfor
Storage_R0        = [];

for i=1:Nmeas

    measname             = Liste_meas{ind_meas(i)};
    disp(['Pulse' ' ' num2str(i) '/' num2str(Nmeas) ' ' 'started'])

    %a) Get R0 (calculated from the original signal)
    Signal=pulses_signal.(Pulsename).(measname);  % Signal containing time, voltage and current.
    voltage=Signal(:,2);
    current=Signal(:,3);

    % Calculate the derivative of current to detect step
    dI = [0; diff(current)];

    % Find index of largest current step
    [~, stepIndex] = max(abs(dI));

    % Calculate delta V and delta I
    deltaV =abs(voltage(1) - voltage(stepIndex));
    deltaI = abs(current(1) - current(stepIndex));

    % Calculate R0
    R0 = deltaV / deltaI;

    %Used substracted signal to calculate pulse-RO*I
    Substracted          = substracted.(Pulsename).(measname);
    [~,ia]               = unique(Substracted(:,1),'stable'); %fit function needs unique components while ODE.tspan need to be striclty increasing or decreasing
    Substracted          = Substracted(ia,:);


    FitR0                = fit([0;OCV.Qsum(i)],[R0 ;R0] ,'linearinterp');
    vector_linear         =linspace(0,OCV.Qsum(i),length(Substracted(:,1)));

    %Substract R0I to Pulse-OCV (substracted voltage - RO*I)
    R0I                   =  FitR0(vector_linear).*Substracted(:,3);
    Substracted(:,2)      = Substracted(:,2) - R0I;


    %2 - Optimization
    %A) Estimation options

    %parameters for differential equation
    ODE=struct;
    ODE.tspan     = Substracted(:,1);
    ODE.y0        = zeros(EstimationOpts.nRCpairs,1); %
    ODE.nRCpairs  = EstimationOpts.nRCpairs;
    CurrentFit     = fit(Substracted(:,1), Substracted(:,3), 'linearinterp'); %to pass to the ODE solver


    %% B) Optimisation type

    %try PSO or fmincon

    ODE.Solver.name   = "ODE15s"; %otherwise too slow with ODE 45/113
    Anonym_function= @(k)ECMCost_Meatech(k,Substracted,CurrentFit,ODE);
    lb              = [1,3,1e-4,1e-4];
    ub              = [5,10,max(R0),max(R0)];
    %% Particle Swarm algorithm

    %PSO_opts          = optimoptions('particleswarm','HybridFcn',@fmincon,'Display','iter','PlotFcn','pswplotbestf');
    %[OPTOUT,fval,exitflag,output]   = particleswarm(Anonym_function,2*ODE.nRCpairs,lb,ub,PSO_opts);

    %% fmincon


    fmincon_opts=optimoptions('fmincon','PlotFcn','optimplotfval');
    x0 = [ 2, 7,max(R0),max(R0)];
    [OPTOUT,fval,exitflag,output]   = fmincon(Anonym_function,x0,[],[],[],[],lb,ub,[],fmincon_opts);
    %%
    Optimname        = strcat(measname,'Gopt');
    sse              = Anonym_function(OPTOUT);
    ODE.Solver.error  = sse; %if u want to save error later


    % C) Plot ODE
    pltODE_Meatech(OPTOUT,Substracted,CurrentFit,Optimname,foldername,ODE)

    % D) Save iteration data
    Storage_OptParams(i,:) = OPTOUT; %optimized parameters
    Storage_SOC0  (i,:)    = [OCV.SOCinit(i);OCV.SOCend(i)];
    Storage_R0   (i,:)     = [R0;R0]; %R0 vector with 2 elements: same beginning and after pulse
    Storage_lb(i,:)=lb;
    Storage_ub(i,:)=ub;

    disp(['Pulse' ' ' num2str(i) '/' num2str(Nmeas) ' ' 'finished'])

    % % End loop

end
% Save total
Results.EstimationOpts  = EstimationOpts;
Results.Storage_SOC0       = Storage_SOC0;
Results.Storage_R0         = Storage_R0;
Results.Storage_OptParams  = Storage_OptParams;
Results.Storage_lb         = Storage_lb;
Results.Storage_ub         = Storage_ub;
Results.Crate              = pulses.Crate*sign_CurrSign;
Results.Temp               = pulses.Temp;
Results.CurrSign           = CurrSign;

save(strcat(pwd,filesep,'results',filesep,foldername,filesep,'Results.mat'),"Results")
end


