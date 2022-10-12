clear all; close all

% Read scenario file
disp('Select scenario Excel file')
filename = uigetfile(['*.xls; *.xlsx']);

tab = readtable(filename,'ReadVariableNames',0);

% Create output folder
outname =tab.(2){1};
if not(exist(outname,'dir'))
    mkdir(outname)
end

% Weather file
weathername = tab.(2){2};

% Global parameters
k       = tab.(3)(3:6);
gamma   = tab.(3)(7);
eta     = tab.(3)(8);
clay    = tab.(3)(9);

% Initial SOC content
C0      = tab.(3)(10:13);
IOM     = tab.(3)(14);

% Time parameters
yr0     = tab.(3)(15);
T       = tab.(3)(16);

% Input files
YR          = tab.(3)(17:end);
inputname   = tab.(2)(17:end);

% Setting
% tspan   = [t0,t0+Nmon];
% period  = (t0:dt:t0+T)';
avecg   = [gamma 1-gamma 0 0];
avecf   = [eta eta 0 1-2*eta];

x = 1.67 * ( 1.85 + 1.6* exp(-0.0786*clay) );
biohumfac = 1/(x+1);
alpha= 0.46*biohumfac;
beta = 0.54*biohumfac;

weather = readtable(weathername,'ReadVariableNames',0);
temp    = str2double(weather.(2)(2:end));
rain    = str2double(weather.(3)(2:end));
evap    = str2double(weather.(4)(2:end));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

Nsim     = length(YR);
months   = {};
for it = yr0:YR(Nsim)
    months = [months; {'JAN'; 'FEB'; 'MAR'; 'APR'; 'MAY'; 'JUN'; 'JUL'; 'AUG'; 'SEP'; 'OCT'; 'NOV'; 'DEC'}] ;
end

t0=12*yr0;
dt=1;
C0RC  = C0;
C0ERE = C0;
C0NS  = C0;
Ttout  =t0;
TCoutRC  =C0';
TCoutERE =C0';
TCoutNS  =C0';

disp('Simulation started')


for it=1:Nsim
    
    yrend   = YR(it);
    input   = readtable(sprintf('%s',inputname{it}),'ReadVariableNames',0);
    g       = str2double(input.(2)(2:end));
    f       = str2double(input.(3)(2:end));
    cover   = str2double(input.(4)(2:end));
    
    b       = g*avecg + f*avecf;
    b       = b';
    rho     = rhofun(temp,rain,evap,cover,clay);

    tspan   = [t0,12*yrend+12];
    
    [tout,CoutRC]  = RC(tspan,dt,C0RC,alpha,beta,k,[rho;rho(1)],[b, b(:,1)],0);
    [~   ,CoutERE] = ERE(tspan,dt,C0ERE,alpha,beta,k,[rho;rho(1)],[b, b(:,1)],0);
    [~   ,CoutNS]  = NS(tspan,dt,C0NS,alpha,beta,k,[rho;rho(1)],[b, b(:,1)],0);  
    
    Ttout    = [Ttout;    tout(2:end)];
    TCoutRC  = [TCoutRC;  CoutRC(2:end,:)];
    TCoutERE = [TCoutERE; CoutERE(2:end,:)];
    TCoutNS  = [TCoutNS;  CoutNS(2:end,:)];
      
    t0 = tout(end);
    C0RC  = CoutRC(end,:)';
    C0ERE = CoutERE(end,:)';
    C0NS  = CoutNS(end,:)';

end

    
socRC  = sum(TCoutRC,2) +IOM;
socERE = sum(TCoutERE,2)+IOM;
socNS  = sum(TCoutNS,2) +IOM;
    
out_table = table(floor(Ttout(1:end-1)/12),months,TCoutRC(1:end-1,1),TCoutRC(1:end-1,2),TCoutRC(1:end-1,3),TCoutRC(1:end-1,4),...
                    IOM*ones(length(Ttout)-1,1),socRC(1:end-1));
out_table.Properties.VariableNames =  {'year','month','DPM','RPM','BIO','HUM','IOM','SOC'};
writetable(out_table,sprintf('%s\\outputRC.xls',outname))

out_table = table(floor(Ttout(1:end-1)/12),months,TCoutERE(1:end-1,1),TCoutERE(1:end-1,2),TCoutERE(1:end-1,3),TCoutERE(1:end-1,4),...
                    IOM*ones(length(Ttout)-1,1),socERE(1:end-1));
out_table.Properties.VariableNames =  {'year','month','DPM','RPM','BIO','HUM','IOM','SOC'};
writetable(out_table,sprintf('%s\\outputERE.xls',outname))

out_table = table(floor(Ttout(1:end-1)/12),months,TCoutNS(1:end-1,1),TCoutNS(1:end-1,2),TCoutNS(1:end-1,3),TCoutNS(1:end-1,4),...
                    IOM*ones(length(Ttout)-1,1),socNS(1:end-1));
out_table.Properties.VariableNames =  {'year','month','DPM','RPM','BIO','HUM','IOM','SOC'};
writetable(out_table,sprintf('%s\\outputNS.xls',outname))

   
%%% Plot SOC

figure()
plot(Ttout/12,socRC,'b',Ttout/12,socERE,'k',Ttout/12,socNS,'--g','LineWidth',2)
legend('RC','ERE','NS')
xlim([yr0,Ttout(end)/12])
xlabel('Year')
ylabel('SOC')
title(sprintf('%s',outname),'Interpreter','none')
savefig( sprintf('%s\\SOC.fig',outname) )


disp(sprintf('Simulation ended with success. Results have been saved in the folder %s', outname))

