function [a,b,c,accTSMD]=rhofun(temp,rain,pet,cover,clay,d,meantemp)

n=length(temp);

% the rate modifying factor for temperature
a = 47.91./ ( 1+ exp( 106.06 ./ (temp-meantemp+106.06/log(46.91)) ) );

% the rate modifying factor for moisture
Mcover = -(20 + 1.3*clay-0.01*clay^2) * d/23;
Mbare= Mcover/1.8;
vec=rain-pet;
accTSMD=zeros(n,1);
i=1;
while vec(i)>0
    i=i+1;
end
for j=i:n
    if cover(j)==0
        accTSMD(j)= min( max( accTSMD(j-1)+vec(j) , Mcover ), 0);
    else
        accTSMD(j)= min( max( accTSMD(j-1)+vec(j) , Mbare ), 0);
    end
end

b=ones(n,1);
indc = cover==0 & accTSMD<=0.444*Mcover;
indb = cover==1 & accTSMD<=Mbare; %0.444 ???
b(indc) = 0.2+ (1-0.2)* (Mcover-accTSMD(indc))./(Mcover-0.444*Mcover);
b(indb) = 0.2+ (1-0.2)* (Mbare-accTSMD(indb))./(Mbare-0.444*Mbare);

c=ones(n,1);
c(cover==1)=0.6;

end
