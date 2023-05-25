%% 
dx=1;
dt=1e-3;

L=30;
tend=0.5;

x=0:dx:L;
t=0:dt:tend;
D0=4;
R0=D0/2;
A0=pi*D0^2/4;
E=225000;
h=0.1;

rho=1;
mu=0.009;

Q=zeros(length(t), length(x));
A=zeros(length(t), length(x));
A(1,:) = A0*ones(1, length(x));

Qin = zeros(length(t),1);
Qin(1:100) = 50*(1-cos(2*pi*10*t(1:100)));
Q(1:end,1)=Qin;

idx=2:length(x)-1;

for n=2:501
    
    Qn=Q(n-1,:);
    An=A(n-1,:);
    
    Qt = mean(Q([n-1 n],:));
    At = An;
    
    Pn = 4/3 * (E*h/R0) .* (1- sqrt(A0./An));
    
    Am1 = [-Qn(3:end-1)./(An(3:end-1).*dx) 0];
    Ap1 = [0 Qn(2:end-2)./(An(2:end-2).*dx)];
    
    R = 8*mu*pi./An.^2;
    
    Ac = 2/dt + R(2:end);
    
    Amat=diag(Am1, -1) + diag(Ac, 0) + diag(Ap1, 1);
    
    B=zeros(size(Ac));
    B(2:end-1) = (-Qn(3:end-1).^2 ./ (2*dx)) .* ( 1./An(4:end) - 1./An(2:end-2)) + ...
        -An(3:end-1)./(2*rho*dx) .* (Pn(4:end) - Pn(2:end-2)) + ...
        Qn(3:end-1)./(2*dt);
    
    B(1) = (-Qn(2)^2 / (2*dx)) * ( 1./An(3) - 1./An(1)) + ...
        -An(2)./(2*rho*dx) .* (Pn(3) - Pn(1)) + Qt(1)*Qn(2)/(An(2)*dx) + ...
        Qn(2)./(2*dt);
    
    B(end) = (-Qn(end)^2 / (2*dx)) * ( 3./An(end) - 4./An(end-1) + 1./An(end-2) ) + ...
        -An(end)./(2*rho*dx) .* ( 3*Pn(end) - 4*Pn(end-1) + Pn(end-2) ) + ...
        Qn(end)./(2*dt);
    
    Qt(2:end) = B'\Amat;
    
    At(2:end-1) = -dt/(2*dx) * ( Qn(3:end)-Qn(1:end-2) ) + An(2:end-1);
    At(1) = -dt/(2*dx) * (-Qn(3) + 4*Qn(2) - 3*Qn(1) ) + An(1);
    
    Pt = 4/3 * (E*h/R0) .* (1- sqrt(A0./At));
    
    R = 8*mu*pi./At.^2;
    
    Q(n,2:end-1) = dt * ( -R(2:end-1).*Qt(2:end-1) + ...
        - Qt(2:end-1)./(At(2:end-1)*dx) .* ( Qt(3:end) - Qt(1:end-2) ) + ...
        - Qt(2:end-1).^2 .* ( 1./At(3:end) - 1./At(1:end-2) ) + ...
        -At(2:end-1)./(2*rho*dx) .* ( Pt(3:end) - Pt(1:end-2) ) )+ ...
        Qn(2:end-1);
    Q(n,end) = dt * ( -R(end).*Qt(end) + ...
        - Qt(end).^2 .* ( 3./At(end) - 4./At(end-1) + 1./At(end-2)  ) + ...
        -At(end)./(2*rho*dx) .* ( 3*Pt(end) - 4*Pt(end-1) + Pt(end-2)) )+ ...
        Qn(end);
    
    A(n,2:end-1) = -dt/(2*dx) * ( Qt(3:end)-Qt(1:end-2) ) + An(2:end-1);
    A(n,1) = -dt/(2*dx) * (-Qt(3) + 4*Qt(2) - 3*Qt(1) ) + An(1);
    A(n,end)=A0;
end