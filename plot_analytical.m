%% Plot Telegrapher Eqn

L=60;
tend=1;

dx=0.01;
dt=0.01;

x=0:dx:L;
t=0:dt:tend;
D0=4;
R0=D0/2;
A0=pi*D0^2/4;
E=3600000;
h=0.1;

rho=1;
mu=0.009;

R = 8*mu*pi./A0^2;
I = rho/A0;
C = A0*D0/(E*h);

beta = I/L;
d=beta/2;
c=sqrt(1/(I*C));

P=zeros(size(t));

for n=1:100
   
    wn = sqrt( (n*pi/L)^2 * c^2 - beta^2/4);
    phin = atan(d/wn);
    An = (2/L*cos(phin)) * sin(n*pi/2);
    
    P = 0.5 * An * exp(-d.*t) .* ( sin( n*pi/2 - wn.*t + phin ) + ...
        sin( n*pi/2 + wn.*t - phin ) ) +P;
    
    
end

 plot(t,P); 
%% Plot lumped parameter

R = L*8*mu*pi./A0^2;
I = L*rho/A0;
C = L*A0*D0/(E*h);

beta = R/I;
d=beta/2;
c=sqrt(1/(I*C));

w=sqrt(beta^2-4*c^2);

P = 2*exp(-d.*t) .* ( cosh(0.5*t*w) + (1+beta)/w * sinh(0.5*t*w));

plot(t, P)
plot(t, c*t)

%% Plot lossless

L=30;
tend=0.5;

dx=0.01;
dt=0.01;

x=0:dx:L;
t=0:dt:tend;
D0=4;
R0=D0/2;
A0=pi*D0^2/4;
E=3600000;
h=0.1;

rho=1;
mu=0.009;

R = 8*mu*pi./A0^2;
I = rho/A0;
C = 0.5*A0*D0/(E*h);

beta = I/L;
d=beta/2;
c=sqrt(1/(I*C));

Z0=rho*c/A0;

Qin = zeros(size(t));
Qin(t<=0.1)=100;
plot(t,Qin);

Qout = zeros(size(t));
Qout(t>=0.1 & t<=0.2) = 200;
Qout(t>=0.3 & t<=0.4) = -200;


hold on
plot(t,Qout)

Qout_LC = 100*( 1-cos(t.*c/L) - heaviside(t-0.1).*( 1-cos((t-0.1).*c/L) ) );

plot(t,Qout_LC)
plot(1e-3:1e-3:0.5, varsave_w_P(:,3:4))

%%
Pin=zeros(size(t));

Pin(t<=0.1)=100*Z0;
Pin(t>=0.2 & t<=0.3) = -200*Z0;
Pin(t>=0.4 & t<=0.5) = 200*Z0;

figure
plot(t,Pin)
hold on
%

Pin_LC = 100*sqrt(I/C)*( sin(t.*c/L) - heaviside(t-0.1).*( sin((t-0.1).*c/L) ) );

%plot(t(2:end),L*I*diff(Qout_LC)./dt)
plot(t,Pin_LC)
%plot(1e-3:1e-2:0.5, varsave_w_P(:,7:8))

CP1=[Qin' Qout' Pin' zeros(size(Qin'))];
CP2=[Qin' Qout_LC' Pin_LC' zeros(size(Qin'))];
%% for RCR
clc

syms R C Z Cp L w

i=sym(sqrt(-1));

ZRCR = (R+Z+i*R*Z*C)/(1+i*w*R*C)

Z_cRCR = 1/(i*w*Cp/2 + 1/ZRCR)

Z_cRCR = simplify(Z_cRCR)

Z_tot = 1/(i*w*Cp/2 + 1/(i*w*L + Z_cRCR) );

Z_tot = simplify(Z_tot)

%% Plot Wave w RCR

L=13;
tend=1;

dx=0.01;
dt=0.01;

x=0:dx:L;
t=0:dt:tend;
D0=1.59;
R0=D0/2;
A0=pi*D0^2/4;
E=12720000;
h=0.03;

rho=1;
mu=0.009;

R = 8*mu*pi./A0^2;
I = rho/A0;
C = A0*D0/(E*h);

beta = I/L;
d=beta/2;
c=sqrt(1/(I*C))

Z0=rho*c/A0;

Q_in= zeros(size(t));
Q_wave_out = zeros(size(t));
P_wave_in = Q_wave_out;
P_wave_out = Q_wave_out;

Z=Z0;
C=0.00021;
Rd=2665-Z0;

for f=0:5
    
    w=2*pi*f;

    Z_RCR = (Rd+Z+sqrt(-1)*w*Rd*Z*C)/(1+sqrt(-1)*w*Rd*C) ;
    
    R = (Z_RCR-Z0)/(Z_RCR+Z0)
    R_alt = Rd/(Rd+2*(Z0 + sqrt(-1)*w*Rd*C*Z0) )
    
    scale = 1; if f>1; scale=1/f; end
    
    A = scale*50*exp(sqrt(-1)*w.*(t-0.2))./(exp(sqrt(-1)*w.*t) - R.*exp(sqrt(-1)*w*(t-2*L/c) ) );
    
    Q_in = scale*50*( exp(sqrt(-1)*w.*(t-0.2)) ) + Q_in;
    Q_wave_out = A.*( exp(sqrt(-1)*w.*(t-L/c)) - R*exp(sqrt(-1)*w.*(t+L/c-2*L/c) ) ) + Q_wave_out;
    
    P_wave_in = Z0*A.*( exp(sqrt(-1)*w.*(t-0/c)) + R*exp(sqrt(-1)*w.*(t+0/c-2*L/c) ) ) + P_wave_in;
    
    P_wave_out = Z0*A.*( exp(sqrt(-1)*w.*(t-L/c)) + R*exp(sqrt(-1)*w.*(t+L/c-2*L/c) ) ) + P_wave_out;
    
end

plot(t,Q_in); hold on; plot(t, Q_wave_out)

figure;
plot(t,P_wave_in./1333.2); hold on; plot(t, P_wave_out./1333.2)

mP = mean(real(P_wave_in)./1333.2)
PP = range(real(P_wave_in)./1333.2)

CP = real([Q_in' Q_wave_out' P_wave_in' P_wave_out']);