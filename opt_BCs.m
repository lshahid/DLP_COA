function [PBF, RHC] = opt_BCs( RC, Qin, tend, dt, Rprox, Z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n=length(RC)/2;

RCR=zeros(1, 3*n);

options = optimoptions('fsolve','Display','off');
for i=1:n
    RCR(3*i-2)=Z;
    RCR(3*i-1)=RC(2*i);
    RCR(3*i)=RC(2*i-1)-Z;
    
    if RCR(3*i)<0
        RCR(3*i-2)=0.1*RC(2*i-1);
        RCR(3*i)=0.9*RC(2*i-1);
    end
        
end

count=0;

for t=0:dt:tend
    count=count+1;
    Q0=Qin(count);
  
    if count==1
        Qtm1=[0 0];
        Ptm1=[0 0 0];
    else
        Qtm1=Qt;
        Ptm1=Pt;
    end
    
fx=@(QR) solve_nonlinRCR( QR, Qtm1, Ptm1, RCR, Q0, dt, Rprox);

[QR, fval]=fsolve(fx, Qtm1(1), options);
error(count)=fval;
figure(1); semilogy(abs(error), 'k');  drawnow


Qt=[QR; Q0-QR];

dQ=(Qt-Qtm1)./dt;

for i=1:length(RCR)/3
    Z=RCR(3*i-2);
    R=RCR(3*i);
    C=RCR(3*i-1);
    RC=R*C;
    
    R1=Rprox(1,i);
    R2=Rprox(2,i);
    I=Rprox(3,i);
    
    Pt(i+1)= ( Qt(i)*(R+Z) + RC*Z*dQ(i) + RC/dt * Ptm1(i+1) ) /  (1+RC/dt) ;
end
Pt(1)= R1*Qt(i)*abs(Qt(i)) + R2*Qt(i) + I*dQ(i) + ...
        ( Qt(i)*(R+Z) + RC*Z*dQ(i) + RC/dt * Ptm1(i+1) ) /  (1+RC/dt) ;
  
Q(count, :)=Qt;
P(count, :)=Pt;
    
end

Plast=P(round(0.75*end):end,:);
Qlast=Q(round(0.75*end):end,:);
figure(2)
plot([Qin Q]);

figure(3); plot(P./1333.2)

PBF=sum(Qlast(:,1))/sum(Qlast, 'all');

Psys=max(Plast(:,1))/1333.2;
Pdia=min(Plast(:,1))/1333.2;
%Pavg=mean(Plast)/1333.2;
P_LPAs=max(Plast(:,3))/1333.2;
P_LPAd=min(Plast(:,3))/1333.2;

%RHC=[Psys Pdia P_LPAs P_LPAd];
RHC=[Psys Pdia P_LPAs];

output=[Psys Pdia P_LPAs P_LPAd PBF];

%display(num2str([output RCR]));

end

function [error] = solve_nonlinRCR( QR, Qtm1, Ptm1, RCR, Q0, dt, Rprox)
%UNTITL ED2 Summary of this function goes here
%   Detailed explanation goes here 

Pt1=zeros(1,length(RCR)/3);

Q=[QR; Q0-QR];

dQ=(Q-Qtm1)./dt;

for i=1:length(RCR)/3
    Z=RCR(3*i-2);
    R=RCR(3*i);
    C=RCR(3*i-1);
    RC=R*C;
    
    R1=Rprox(1,i);
    R2=Rprox(2,i);
    I=Rprox(3,i);
    
    Pt1(i)= R1*Q(i)*abs(Q(i)) + R2*Q(i) + I*dQ(i) + ...
        ( Q(i)*(R+Z) + RC*Z*dQ(i) + RC/dt * Ptm1(i+1) ) /  (1+RC/dt) ;
end

    error=(Pt1(1)-Pt1(2))/1333.2;
    
end