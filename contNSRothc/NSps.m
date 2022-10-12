function Cstar=NSps(period,dt,alpha,beta,k,rho,b)

Lam	= [0 0 0 0; 0 0 0 0; alpha alpha alpha alpha; beta beta beta beta];
A   = (Lam-eye(4))*diag(k);
invA    = inv(A);
M       = eye(4)-Lam;
Minv    = inv(M);
At      = -M*diag(k)*Minv;
invAt   = -M*diag(1./k)*Minv;
N   = length(period)-1;

    F = zeros(4*N); 
    tn= zeros(4*N,1);
    for i=0:N-2
        F(4*i+1:4*i+4,4*i+1:4*i+4)=Lam+(eye(4)-Lam)*diag(exp(-dt*rho(i+1)*k));
        F(4*i+1:4*i+4,4*i+5:4*i+8)=-eye(4);
        tn(4*i+1:4*i+4)=invAt*(expm(dt*rho(i+1)*At)-eye(4))*b(:,i+1)/rho(i+1);
    end
    F(4*N-3:4*N,4*N-3:4*N)=Lam+(eye(4)-Lam)*diag(exp(-dt*rho(N)*k));
    F(4*N-3:4*N,1:4)=-eye(4);
    tn(4*N-3:4*N)=invAt*(expm(dt*rho(N)*At)-eye(4))*b(:,N)/rho(N);
    Cstar=F\tn;
    Cstar=-Cstar;
   
    
    figure()
    subplot(2,2,1)
    plot(period,Cstar([1:4:end,1]),'g','LineWidth', 2)
    xlabel('Month')
    xticks(period(1):period(end))
    xlim([period(1),period(end)])
    xticklabels({'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D','J'})
    title('DPM')
    
    subplot(2,2,2)
    plot(period,Cstar([2:4:end,2]),'g','LineWidth', 2)
    xlabel('Month')
    xticks(period(1):period(end))
    xlim([period(1),period(end)])
    xticklabels({'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D','J'})
    title('RPM')
    
    subplot(2,2,3)
    plot(period,Cstar([3:4:end,3]),'g','LineWidth', 2)
    xlabel('Month')
    xticks(period(1):period(end))
    xlim([period(1),period(end)])
    xticklabels({'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D','J'})
    title('BIO')
    
    subplot(2,2,4)
    plot(period,Cstar([4:4:end,4]),'g','LineWidth', 2)
    xlabel('Month')
    xticks(period(1):period(end))
    xlim([period(1),period(end)])
    xticklabels({'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D','J'})
    title('HUM')

end