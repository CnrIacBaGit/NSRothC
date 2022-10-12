clear all; close all

% Read scenario file
disp('Select scenario Excel file')
filename = uigetfile(['*.xls; *.xlsx']);

tab = readtable(filename,'ReadVariableNames',0);
disp(sprintf('%s file read', filename))

% Create output folder
outname =tab.(2){1};
if not(exist(outname,'dir'))
    mkdir(outname)
end

% Global parameters
k       = tab.(3)(2:5);
alpha   = tab.(3)(6);
beta    = tab.(3)(7);
gamma   = tab.(3)(8);
eta     = tab.(3)(9);

% Input parameters
ar      = tab.(3)(10); 
mur1    = tab.(3)(11); 
mur2    = tab.(3)(12);
br1     = tab.(3)(13); 
br2     = tab.(3)(14);
sigmar1 = tab.(3)(15); 
sigmar2 = tab.(3)(16);
ag      = tab.(3)(17);  
mug     = tab.(3)(18); 
bg      = tab.(3)(19); 
sigmag  = tab.(3)(20);
af      = tab.(3)(21);  
muf     = tab.(3)(22); 
bf      = tab.(3)(23); 
sigmaf  = tab.(3)(24);

% Initial SOC content
C0      = tab.(3)(25:28);
IOM     = tab.(3)(29);

% Time parameters
t0      = tab.(3)(30);
T       = tab.(3)(31);
Nmon    = tab.(3)(32);
dt      = tab.(3)(33);


% Setting
tspan   = [t0,t0+Nmon];
period  = (t0:dt:t0+T)';
avecg   = [gamma 1-gamma 0 0];
avecf   = [eta eta 0 1-2*eta];

% rho and b
    G       = @(t,mu,sigma) 1/sigma/sqrt(2*pi)*exp(-(t-mu-floor((t-t0)/T)*T).^2/2/sigma^2);
    int1    = integral(@(t) G(t,mur1,sigmar1),t0,t0+T);
    int2    = integral(@(t) G(t,mur2,sigmar2),t0,t0+T);
    intg    = integral(@(t) G(t,mug,sigmag),t0,t0+T);
    intf    = integral(@(t) G(t,muf,sigmaf),t0,t0+T);
    rho     = ar+ G(period,mur1,sigmar1)/int1*br1 + G(period,mur2,sigmar2)/int2*br2;
    g       = ag+G(period,mug,sigmag)/intg*bg;
    f       = af +G(period,muf,sigmaf)/intf*bf;
    b       = g*avecg + f*avecf;
    
    out_table = table(period(1:end-1),rho(1:end-1),g(1:end-1),f(1:end-1));
    out_table.Properties.VariableNames =  {'time','rho','plant_input','FYM'};
    writetable(out_table,sprintf('%s\\input.xls',outname))        

    figure()
    plot([period(1:end-1);period(1:end-1)+T;period+2*T],[rho(1:end-1);rho(1:end-1);rho])
    title('rho')
    xlabel('Month')
    xlim([period(1),period(end)+2*T])
    xticks(period(1):period(end)+2*T)
    xticklabels({'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
        'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
        'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D','J'})
    savefig( sprintf('%s\\rho.fig',outname) )
    
    figure()
    plot([period(1:end-1);period(1:end-1)+T;period+2*T],[b(1:end-1,1);b(1:end-1,1);b(:,1)],'b',...
        [period(1:end-1);period(1:end-1)+T;period+2*T],[b(1:end-1,2);b(1:end-1,2);b(:,2)],'--r',...
        [period(1:end-1);period(1:end-1)+T;period+2*T],[b(1:end-1,3);b(1:end-1,3);b(:,3)],'k',...
        [period(1:end-1);period(1:end-1)+T;period+2*T],[b(1:end-1,4);b(1:end-1,4);b(:,4)],'g')
    legend('b_1','b_2','b_3','b_4')
    title('b')
    xlabel('Month')
    xlim([period(1),period(end)+2*T])
    xticks(period(1):period(end)+2*T)
    xticklabels({'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
        'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
        'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D','J'})
    savefig( sprintf('%s\\b.fig',outname) )

    b       = b';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Periodic equilibrium solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % RothC periodic solution
    
    CstarRC=RCps(period,dt,alpha,beta,k,rho,b);
    socRC=CstarRC([1:4:end,1])+CstarRC([2:4:end,2])+CstarRC([3:4:end,3])+CstarRC([4:4:end,4])+IOM;
    
    savefig( sprintf('%s\\periodicRC.fig',outname) )

    out_table = table(period,CstarRC([1:4:end,1]),CstarRC([2:4:end,2]),CstarRC([3:4:end,3]),CstarRC([4:4:end,4]),IOM*period,socRC);
    out_table.Properties.VariableNames =  {'time','DPM','RPM','BIO','HUM','IOM','SOC'};
    writetable(out_table,sprintf('%s\\periodic_solution_RC.xls',outname))        
      
   % ERE periodic solution
    
    CstarERE=EREps(period,dt,alpha,beta,k,rho,b);
    socERE=CstarERE([1:4:end,1])+CstarERE([2:4:end,2])+CstarERE([3:4:end,3])+CstarERE([4:4:end,4])+IOM;
    
    savefig( sprintf('%s\\periodicERE.fig',outname) )
    
    out_table = table(period,CstarERE([1:4:end,1]),CstarERE([2:4:end,2]),CstarERE([3:4:end,3]),CstarERE([4:4:end,4]),IOM*period,socERE);
    out_table.Properties.VariableNames =  {'time','DPM','RPM','BIO','HUM','IOM','SOC'};
    writetable(out_table,sprintf('%s\\periodic_solution_ERE.xls',outname))
    
    % NS periodic solution
   
    CstarNS=NSps(period,dt,alpha,beta,k,rho,b);
    socNS=CstarNS([1:4:end,1])+CstarNS([2:4:end,2])+CstarNS([3:4:end,3])+CstarNS([4:4:end,4])+IOM;
    
    savefig( sprintf('%s\\periodicNS.fig',outname) )
    
    out_table = table(period,CstarNS([1:4:end,1]),CstarNS([2:4:end,2]),CstarNS([3:4:end,3]),CstarNS([4:4:end,4]),IOM*period,socNS);
    out_table.Properties.VariableNames =  {'time','DPM','RPM','BIO','HUM','IOM','SOC'};
    writetable(out_table,sprintf('%s\\periodic_solution_NS.xls',outname))
    
       
figure()
plot(period,socRC,'b',period,socERE,'k',period,socNS,'--g','LineWidth', 2)
xlabel('Month')
xticks(period(1):period(end))
xlim([period(1),period(end)])
xticklabels({'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D','J'})
title('periodic SOC')
legend('RC','ERE','NS')
savefig( sprintf('%s\\SOCperiodic.fig',outname) )   

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

disp('Simulation started')


    % RC run
    [tout,CoutRC]=RC(tspan,dt,C0,alpha,beta,k,rho,b,1);
    socRC=sum(CoutRC,2);

    savefig( sprintf('%s\\runNS.fig',outname) )

    out_table = table(tout,CoutRC(:,1),CoutRC(:,2),CoutRC(:,3),CoutRC(:,4),IOM*ones(length(tout),1),socRC);
    out_table.Properties.VariableNames =  {'time','DPM','RPM','BIO','HUM','IOM','SOC'};
    writetable(out_table,sprintf('%s\\run_RC.xls',outname))


    % ERE run
    [tout,CoutERE]=ERE(tspan,dt,C0,alpha,beta,k,rho,b,1);
    socERE=sum(CoutERE,2);

    savefig( sprintf('%s\\runERE.fig',outname) )

    out_table = table(tout,CoutERE(:,1),CoutERE(:,2),CoutERE(:,3),CoutERE(:,4),IOM*ones(length(tout),1),socERE);
    out_table.Properties.VariableNames =  {'time','DPM','RPM','BIO','HUM','IOM','SOC'};
    writetable(out_table,sprintf('%s\\run_ERE.xls',outname))


    % NS run
    [tout,CoutNS]=NS(tspan,dt,C0,alpha,beta,k,rho,b,1);
    socNS=sum(CoutNS,2);

    savefig( sprintf('%s\\runNS.fig',outname) )

    out_table = table(tout,CoutNS(:,1),CoutNS(:,2),CoutNS(:,3),CoutNS(:,4),IOM*ones(length(tout),1),socNS);
    out_table.Properties.VariableNames =  {'time','DPM','RPM','BIO','HUM','IOM','SOC'};
    writetable(out_table,sprintf('%s\\run_NS.xls',outname))


N=length(period);
NT=length(tout);
np=(NT-1)/N;

figure()
plot(tout,socRC,'b',tout,socERE,'k',tout,socNS,'--g','LineWidth',2)
legend('RC','ERE','NS')
title('SOC')
xlim(tspan)
xticks(tout(1:N:end))
xticklabels(split(num2str(0:np)))
xlabel('Year')
savefig( sprintf('%s\\SOCrun.fig',outname) )


disp(sprintf('Simulation ended with success. Results have been saved in the folder %s', outname))
