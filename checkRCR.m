%% Check RCR
clc
t=VarTable.t;
Q=VarTable.Q3;
PDLP=VarTable.Pseg3;
Q(1)=0;
dt=t(2)-t(1);

BCID=3;

Z=BCs(BCID).RCR(1);
        C=BCs(BCID).RCR(2);
        R=BCs(BCID).RCR(3);
        RC=R*C;

P=0;
for i=2:length(t)
    dQ(i) = (Q(i)-Q(i-1))/dt;
    P(i)=( Q(i)*(R+Z) + RC*Z*dQ(i) + RC/dt*P(i-1) ) / (1+RC/dt);
    
end

figure; plot(Q); figure; plot(dQ); figure
plot(t, PDLP, 'linewidth', 2)
hold on
plot(t, P/1333)
