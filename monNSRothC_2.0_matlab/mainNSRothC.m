clear all; close all

% Read scenario
disp('Select scenario Excel file')
filename = uigetfile(['*.xls; *.xlsx']);
tab = readtable(filename,'ReadVariableNames',0);

outname =tab.(2){1};
dataname = tab.(2){2};
parname = tab.(2){3};

% Create output folder
if not(exist(outname,'dir'))
    mkdir(sprintf('OUTPUT_%s',outname))
end

% Read files
datatab   = readtable(dataname,'ReadVariableNames',0);
partab  = readtable(parname,'ReadVariableNames',0);

% Global parameters
k       = partab.(3)(1:4);
r       = partab.(3)(5);
eta     = partab.(3)(6);
clay    = partab.(3)(7);
d       = partab.(3)(8);
meantemp= partab.(3)(9);

% Initial SOC content
C0      = partab.(3)(10:13);
IOM     = partab.(3)(14);

% Read data
years   = datatab.(1);
months  = datatab.(2);
temp    = datatab.(3);
rain    = datatab.(4);
pet     = datatab.(5);
g       = datatab.(6);
f       = datatab.(7);
cover   = datatab.(8);

% Settings
gamma   = r/(r+1);
avecg   = [gamma 1-gamma 0 0];
avecf   = [eta eta 0 1-2*eta];

b         = g*avecg + f*avecf;
b         = b';

x = 1.67 * ( 1.85 + 1.6* exp(-0.0786*clay) );
biohumfac = 1/(x+1);
alpha= 0.46*biohumfac;
beta = 0.54*biohumfac;

% Time parameters
t0=1;
T =length(years);
tspan   = [t0,T];
%dt=1; %month

% Rate modifying factors
[ka,kb,kc,acc] = rhofun(temp,rain,pet,cover,clay,d,meantemp);
rho=ka.*kb.*kc;

rho_table = table(years,months,ka,kb,kc,acc);
rho_table.Properties.VariableNames =  {'year','month','ka','kb','kc','acc TSMD'};
writetable(rho_table,sprintf('OUTPUT_%s\\RHO_%s.xls',outname,outname))

% Run
Lam=[0 0 0 0; 0 0 0 0; alpha alpha alpha alpha; beta beta beta beta];
A=-(eye(4)-Lam)*diag(k);
M=eye(4)-Lam;
Minv=inv(M);
At=-M*diag(k)*Minv;
Atinv=-M*diag(1./k)*Minv;

Cout=C0;
tout=t0:T;
for n=tout(1:end-1)
    if rho(n)==0
        C0=C0 + b(:,n);
    else
        F = Atinv  * ( expm( rho(n)*At) -eye(4) ) / rho(n);
        C0=C0 + F*(rho(n)*A*C0+ b(:,n));
    end
    Cout=[Cout, C0];
end
tout=tout';
Cout=Cout';

soc  = sum(Cout,2) +IOM;
  
% Save table
out_table = table(years,months,Cout(:,1),Cout(:,2),Cout(:,3),Cout(:,4),...
                    IOM*ones(T,1),soc);
out_table.Properties.VariableNames =  {'year','month','DPM','RPM','BIO','HUM','IOM','SOC'};
writetable(out_table,sprintf('OUTPUT_%s\\SOC_%s.xls',outname,outname))

   
%%% Plot SOC
figure()
plot(years,soc,'LineWidth',2)
xlabel('Year')
ylabel('SOC')
title(sprintf('%s',outname),'Interpreter','none')
savefig( sprintf('OUTPUT_%s\\SOC_%s.fig',outname,outname) )


disp(sprintf('Simulation ended with success. Results have been saved in the folder %s', outname))

