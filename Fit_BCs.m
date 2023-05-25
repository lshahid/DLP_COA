%% Fit BCs
RHC=[28 10 12 6];
goal_RHC = RHC(1:3)-RHC(4);
% 
% if RHC(3)-RHC(4)==0
%     goal_RHC(4)=goal_RHC(3)-1;
% end

goal_PBF = [0.8];
P_loc=91;
target.Qpercent = [0.89 0.11];
Q_loc=9;

Pin=table2array(VarTable(:,P_loc));
Qin=table2array(VarTable(:,Q_loc));
t=table2array(VarTable(:,1));
SV=trapz(t(round(.75*end):end), Qin(round(.75*end):end));
%%
BC_segs=[];
BC_def=[];
BCidout=[];
group=[];
for i=1:length(BCs)
    BC_segs(i)=BCs(i).segs;
    BC_def(i)=BCs(i).def;
end
BCidout=find(BC_def==1);

VarNames=VarTable.Properties.VariableNames;

group=zeros(size(BCidout));

for i=1:length(BCidout)
    group(i) =  BCs(BCidout(i)).group;
end

Qout=zeros(length(VarTable.t),max(group));
Pout=Qout;
for i=1:max(group)
    BCgroup=BCidout(group==i);
    idxQ=zeros(size(BCgroup));
    idxP=zeros(size(BCgroup));
    for j=1:sum(group==i)
    idxQ(j)=1+BCs(BCgroup(j)).segs;
    idxP(j)=1+BCgroup(j)+length(segment)+length(junctions);
    end
    Qout(:,i) = sum( table2array(VarTable(:,idxQ)), 2 );
    
    Pout(:,i)= mean( table2array(VarTable(:,idxP)), 2 );
end

dQout=[0 0; diff(Qout)]./dt;

plot([Pin Pout])
figure
plot([Qin Qout])
%%
Pgrad=repmat(Pin, [1 max(group)])-Pout;
Rprox=zeros(3,max(group));
for i=1:max(group)
    
    A=[abs(Qout(:,i)).*Qout(:,i) Qout(:,i) dQout(:,i)];
    b=Pgrad(:,i);
    x=lsqlin(A, b, [], []);
    Rprox(:,i)=x;
    
    subplot(1,max(group),i)
    plot(t, Pgrad(:,i)./1333.2, 'k.')
    hold on
   
    plot(t, (x(1).*Qout(:,i).^2+x(2).*Qout(:,i)+x(3).*dQout(:,i))./1333.2, 'r');
end
%%
[filename, pathname]=uigetfile('*.txt','Select RCR File');
RC = readtable ([pathname filename]);

%ZR_ratio=input('Type Z:R ratio: ');
InflowBCID= find(BC_def==0);
PWV=input('Type PWV (cm/s): ');

Z=1.06*PWV/( BCs(InflowBCID(1)).A /100 );


OutflowID=BC_segs(BC_def==1);
OutflowBCID= find(BC_def==1);
Atotal=zeros(length(RC.R),1);
for i=1:length(OutflowID)
    Atotal(BCs(OutflowBCID(i)).group) = Atotal(BCs(OutflowBCID(i)).group) + ...
     BCs(OutflowBCID(i)).A  ;   
end

RC0=zeros(1,2*length(RC.R));
RT=0;
for i=1:length(RC.R)
    RC0(2*i-1:2*i)=[RC.R(i) RC.C(i)]';
end

RT=1/(1/RC0(1)+1/RC0(3));
CT=RC0(2)+RC0(4);
%R2_R1=RC0(3)/RC0(1);
R1_R2=RC0(1)/RC0(3);
C1_C2=RC0(2)/RC0(4);
%% Fit to PBF
RC_PBF=RC0;

clc
dt=0.01;
tend=4*BCs(2).Flow.t(end);

for i=1:50

[PBF, RHC] = opt_BCs( RC_PBF, Qin, tend, dt, Rprox, Z);

display([PBF goal_PBF])

if abs( (PBF - goal_PBF ) / goal_PBF ) <0.1
    break
end

RC_PBF(1)= RC_PBF(1)/(1-(PBF - goal_PBF)/(goal_PBF) );
if RC_PBF(1)<100
    RC_PBF(1)=100;
elseif RC_PBF(1)>1e4
    RC_PBF(1)=1e4;
end
RC_PBF(2)=0.5/RC_PBF(1);
RC_PBF(3)= 1/(1/RT-1/RC_PBF(1));
if RC_PBF(3)<100
    RC_PBF(3)=100;
elseif RC_PBF(3)>1e4
    RC_PBF(3)=1e4;
end
RC_PBF(4)=0.5/RC_PBF(3);

end

R=[RC_PBF(1); RC_PBF(3)];
C=[RC_PBF(2); RC_PBF(4)];
RC_PBF=table(R, C);

BC_PBF=BCs;

 Atotal=zeros(length(RC.R),1);
for i=1:length(OutflowID)
   
    Atotal(BCs(OutflowBCID(i)).group) = Atotal(BCs(OutflowBCID(i)).group) + ...
     BCs(OutflowBCID(i)).A  ;
      
end

for i=1:length(OutflowID)
    Asplit = BCs(OutflowBCID(i)).A / Atotal(BCs(OutflowBCID(i)).group);
    Rtemp= RC_PBF.R(BCs(OutflowBCID(i)).group)/Asplit ;
    Ctemp= RC_PBF.C(BCs(OutflowBCID(i)).group)*Asplit ;
    
        Z>Rtemp
        
        BC_PBF(OutflowBCID(i)).RCR = [Z Ctemp Rtemp-Z]   ;
end

%% Fit to RHC
RC_cum=[];
error_cum=[];
clc
dt=0.01;
tend=4*BCs(2).Flow.t(end);
MPAs_goal=goal_RHC(1);
MPAd_goal=goal_RHC(2);
Pgrad_goal=goal_RHC(3);

RC_RHC=RC0;

for i=1:50

[PBF, RHC] = opt_BCs( RC_RHC, Qin, tend, dt, Rprox, Z);

error=(goal_RHC-RHC)./RHC;

display([goal_RHC-RHC sum(abs(goal_RHC-RHC))])

error_cum(i)=sum(abs(goal_RHC-RHC));
RC_cum(i,:)=RC_RHC;

if max(abs(error)) <0.1
    break
end

if max(abs(goal_RHC-RHC)) <1
    break
end

if i>20 && error_cum(i)-error_cum(i-1)<0.05
   break 
end

RC_RHC(1) = RC_RHC(1)*(1+( goal_RHC(2) - RHC(2) )/goal_RHC(2) ); 

if RC_RHC(1)<100
    RC_RHC(1)=100;
elseif RC_RHC(1)>1e4
    RC_RHC(1)=1e4;
end

RC_RHC(2) = RC_RHC(2)*(1-( (goal_RHC(1)-goal_RHC(2)) - (RHC(1)-RHC(2)) )/(goal_RHC(1)-goal_RHC(2))); 

RC_R=RC_RHC(1)*RC_RHC(2);

if RC_RHC(2)>0.1
    RC_RHC(2)=0.1;
end

RC_RHC(3) = RC_RHC(3)*(1+( goal_RHC(3) - RHC(3) )/goal_RHC(3) ); 

if RC_RHC(3)<100
    RC_RHC(3)=100;
    elseif RC_RHC(3)>1e4
    RC_RHC(3)=1e4;
end

RC_RHC(4)= 0.5/RC_RHC(3);
%RC_RHC(4) = RC_RHC(4)*(1-( (goal_RHC(3)-goal_RHC(4)) - (RHC(3)-RHC(4)) )/(goal_RHC(3)-goal_RHC(4))); 

% RC0
% RC_RHC
end

RC_RHC=RC_cum(error_cum==min(error_cum),:);

R=[RC_RHC(1); RC_RHC(3)];
C=[RC_RHC(2); RC_RHC(4)];
RC_RHC=table(R, C);

BC_RHC=BCs;

 Atotal=zeros(length(RC.R),1);
for i=1:length(OutflowID)
   
    Atotal(BCs(OutflowBCID(i)).group) = Atotal(BCs(OutflowBCID(i)).group) + ...
     BCs(OutflowBCID(i)).A  ;
      
end

for i=1:length(OutflowID)
    Asplit = BCs(OutflowBCID(i)).A / Atotal(BCs(OutflowBCID(i)).group);
    Rtemp= RC_RHC.R(BCs(OutflowBCID(i)).group)/Asplit ;
    Ctemp= RC_RHC.C(BCs(OutflowBCID(i)).group)*Asplit ;
    
    if Z>Rtemp
        1;
    end
        
        BC_RHC(OutflowBCID(i)).RCR = [Z Ctemp Rtemp-Z]   ;
end
