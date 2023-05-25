function [output] = opt_BCs_quadratic( RC, Qin, tend, dt, Rprox, ZR_ratio)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n=length(RC)/2;

RCR=zeros(1, 3*n);

options = optimoptions('fsolve','Display','off');
for i=1:n
    RCR(3*i-2)=ZR_ratio*RC(2*i-1);
    RCR(3*i-1)=RC(2*i);
    RCR(3*i)=(1-ZR_ratio)*RC(2*i-1);
end

count=0;

Qt0=Qin(1)/2;
for t=0:dt:tend
    count=count+1;
    QM=Qin(count);
  
    if count==1
        Qtm1=zeros(1,n);
        Ptm1=[0 0 0];
    end

[Ptm1, Qtm1] = solve_Eq( QM, Qtm1, Ptm1, RCR, dt, Rprox);

P(count,:)= Ptm1;
Q(count,:)= Qtm1;

end

Plast=P(round(0.75*end):end);
Qlast=Q(:,round(0.75*end):end);
figure(2)
plot(t, Q(:,1)); hold on; plot(t, Q(:,2));plot(t, Qin)

Qperc=sum(Qlast(1,:))/sum(Qlast, 'all');

Psys=max(Plast)/1333.2;
Pdia=min(Plast)/1333.2;
%Pavg=mean(Plast)/1333.2;
Pavg=(Psys+2*Pdia)/3;

output=[Psys Pdia Pavg Qperc];

display(num2str([output RCR]));

end

function [P, Q] = solve_Eq( QM, Qtm1, Ptm1, RCR, dt, Rprox)
%UNTITL ED2 Summary of this function goes here
%   Detailed explanation goes here 

%%
QR_tm1= Qtm1(1);
QL_tm1= Qtm1(2);

PR_tm1=Ptm1(2);
PL_tm1=Ptm1(3);

R2_R=Rprox(1,1);
R1_R=Rprox(2,1);
I_R=Rprox(2,1);

Z_R=RCR(1);
C_R=RCR(2);
R_R=RCR(3);
RC_R=R_R*C_R;

R2_L=Rprox(1,2);
R1_L=Rprox(2,2);
I_L=Rprox(2,2);

Z_L=RCR(4);
C_L=RCR(5);
R_L=RCR(6);
RC_L=R_L*C_L;

%%
AT=R2_R+R2_L;

BT=( R1_R + I_R/dt + 1/(1+RC_R/dt)*(R_R+Z_R) + RC_R*Z_R/dt ) + ...
    ( R1_L + I_L/dt + 1/(1+RC_L/dt)*(R_L+Z_L) + RC_L*Z_L/dt ) + ...
    2*R2_L*QM;

CT=( -(I_R + RC_R*Z_R/(1+RC_R/dt))*QR_tm1/dt + (RC_R/(dt+RC_R))*PR_tm1 ) + ...
   - ( -(I_L + RC_L*Z_L/(1+RC_L/dt))*QL_tm1/dt + (RC_L/(dt+RC_L))*PL_tm1 ) +...
   -R2_L*QM^2 - ( R1_L + I_L/dt + 1/(1+RC_L/dt)*(R_L+Z_L) + RC_L*Z_L/dt )*QM;

%%
QR=(-BT - sqrt(BT^2 - 4*AT*CT) )/ (2*AT);
QL=QM-QR;

PR= 1/(1+RC_R/dt) * ( QR*(R_R+Z_R) + RC_R*Z_R*(QR-QR_tm1)/dt + RC_R/dt * PR_tm1 );
PL= 1/(1+RC_L/dt) * ( QL*(R_L+Z_L) + RC_L*Z_L*(QL-QL_tm1)/dt + RC_L/dt * PL_tm1 );

PM1= PR + R2_R*QR^2 + R1_R*QR + I_R/dt * (QR - QR_tm1)/dt;
PM2= PL + R2_L*QL^2 + R1_L*QL + I_L/dt * (QL - QL_tm1)/dt;

PM=(PM1 + PM2)/2;

P=[PM PR PL];
Q=[QR QL];
    
end