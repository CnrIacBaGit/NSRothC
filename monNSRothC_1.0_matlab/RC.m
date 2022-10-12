function [tout,Cout]=RC(tspan,dt,C0,alpha,beta,k,rho,b,fig)

Lam=[0 0 0 0; 0 0 0 0; alpha alpha alpha alpha; beta beta beta beta];
N=length(rho)-1;

F=zeros(4,4,N);
for n=1:N
    F(:,:,n)=Lam+(eye(4)-Lam)*expm(-dt*rho(n)*diag(k));
end

    t0=tspan(1); tend=tspan(2); 
    Cout=C0;
    tout=t0; 
    
    while t0 <= tend-N*dt 

       for n=1:N
           
            C0=F(:,:,n)*C0 + dt*b(:,n);
            t0=t0+dt;
            Cout=[Cout, C0];
            tout=[tout, t0];          
       
       end

    end
    
    tout=tout';
    Cout=Cout';
    
   NT=length(tout);
   np=(NT-1)/N;

if fig==1    
       
    figure()
    subplot(2,2,1)
    plot(tout,Cout(:,1),'LineWidth',2)
    title('DPM')
    xlim(tspan)
    xticks(tout(1:N:end))
    xticklabels(split(num2str(0:np)))
    xlabel('Year')
    subplot(2,2,2)
    plot(tout,Cout(:,2),'LineWidth',2)
    title('RPM')
    xlim(tspan)
    xticks(tout(1:N:end))
    xticklabels(split(num2str(0:np)))
    xlabel('Year')
    subplot(2,2,3)
    plot(tout,Cout(:,3),'LineWidth',2)
    title('BIO')
    xlim(tspan)
    xticks(tout(1:N:end))
    xticklabels(split(num2str(0:np)))
    xlabel('Year')
    subplot(2,2,4)
    plot(tout,Cout(:,4),'LineWidth',2)
    title('HUM')
    xlim(tspan)
    xticks(tout(1:N:end))
    xticklabels(split(num2str(0:np)))
    xlabel('Year')

end

end

