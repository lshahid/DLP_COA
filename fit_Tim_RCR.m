dt=0.01;
t=0:dt:4;
CO=4*1000/60;
Q=CO*(1-cos(2*pi*t));
dQdt=CO*2*pi*sin(2*pi*t);

plot(t,Q)
hold on
plot(t,dQdt)
%%
E=10*2.065e6;
h=3.2;
D=19.1;
A=pi*D^2/4;
rho=1;
PWV=sqrt(E*h/(D*rho));
Z=1.06*PWV/( A/100 );

%%
R=1405;
Z=0;
C=0.00016;
RC=R*C;


P=zeros(1,length(t));
for i=2:length(t)
    
    P(i)= 1/(1+RC/dt) * ( (R+Z)*Q(i) + (RC*Z)*dQdt(i) + RC * P(i-1)/dt) ;
    
end

plot(t,P./1333.2)

mP=mean(P(round(0.75*end):end)./1333.2)
PP=range(P(round(0.75*end):end)./1333.2)