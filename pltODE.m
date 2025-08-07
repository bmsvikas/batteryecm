function [] = pltODE_Meatech(k0,Substracted,CurrentFit,Optimname,foldername,ODE)

% Assess ODE
disp("Ploting results")
    [~, solved] = ode15s(@(t,y)ECMODE_Meatech(t,y,CurrentFit,k0,ODE), ODE.tspan, ODE.y0);

yexp = Substracted(:,2);
yode = sum(solved,2);
ysim = yode; %if yexp=Vpulse-OCV-R0*I;

% Plot results
h=figure("Visible","on");
subplot(2,1,1)
hold on
plot(ODE.tspan, ysim,'LineWidth',3)
plot(Substracted(:,1), yexp, 'LineWidth',3 )
legend("Simulated", "Measured")
xlabel("Time (s)")
ylabel("Voltage (V)")
box on
grid on

subplot (2,1,2)
plot(ODE.tspan, abs(ysim - yexp),'LineWidth',3,'Color','r')
xlabel("Time (s)")
ylabel("|Voltage error (V)|")
box on
grid on
savelocation = strcat(pwd,filesep,'results',filesep,foldername);
figname      = strcat(Optimname,'_',ODE.Solver.name,'test','.fig');
savefig(h,strcat(savelocation,filesep,figname))
end
