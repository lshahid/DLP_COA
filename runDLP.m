function [VarTable] = runDLP(dt, tend, segment, junctions, trim_centerline,mat_props, BCs)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(BCs)
    BC_segs(i)=BCs(i).segs;
    BC_def(i)=BCs(i).def;
end

BCid=find(BC_def==0);
InflowID=BC_segs(BC_def==0);
BCidout=find(BC_def==1);
OutflowID=BC_segs(BC_def==1);

Qt1=zeros(length(segment),1);
Qt=zeros(length(segment),1);

solveQID=1:length(segment);
for i=1:length(InflowID)
solveQID(solveQID==InflowID(i))=[];
end

for i=1:length(junctions)
    ID=junctions(i).segs(end);
    solveQID(solveQID==ID)=[]; 
end

solvePID=1:length(junctions)+length(BCs);
solvePID(length(junctions)+BCidout)=[];

NVar=length(segment);

options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
%options = optimoptions('fsolve','Algorithm','trust-region-dogleg');

error=zeros(NVar);
count=1;

flag=1;



%%
fx0 = @(varT1) solveDLP_steady( varT1, dt, segment, junctions, ...
    trim_centerline, mat_props, BCs );

var0=zeros(NVar,1);
[var0, fval, flag]=fsolve(fx0, var0, options);

Pt1=zeros(length(junctions)+length(BCs), 1);
Pt1(solvePID)=var0(length(solveQID)+1 : end);

Qt1=zeros(length(segment),1);
Qt1(solveQID)=var0(1:length(solveQID));

for i=1:length(InflowID)
    
    flowwave=BCs(BCid(i)).Flow;
    
    Period=flowwave.t(end);
    omega=2*pi/Period;
    Qtemp=interp1(flowwave.t, flowwave.Q, dt);
    
    if i==1
    Qt1(InflowID(i))=Qtemp(1);
    else
    Qt1(InflowID(i))= -Qtemp(1);
    end
end

for i=1:length(junctions)
    ID=junctions(i).segs(end);
    
    Qt1(ID)= -junctions(i).defs(end)* ( junctions(i).defs * Qt1(junctions(i).segs) );
end

for i=1:length(BCidout)
    Pidx=length(junctions)+BCidout(i);
    Qidx=OutflowID(i);
    
    Z=BCs(BCidout(i)).RCR(1);
        C=BCs(BCidout(i)).RCR(2);
        R=BCs(BCidout(i)).RCR(3);
        RC=R*C;

    Q=Qt1(Qidx);
        
    Pt1(Pidx)= Q*(R+Z) ;
end

varsave=[0; Qt1; Pt1; fval];
varT=zeros(size([Qt1; Pt1]));
var0=zeros(size(var0));
%%

for t=dt:dt:tend
count=count+1;
    
fx=@(varT1) solveDLP(varT1, varT, t, dt, segment, junctions, trim_centerline,...
    mat_props, BCs);

[varT1, fval, exitflag]=fsolve(fx, var0, options);

display(['t=' num2str(t) ' Max Error is ' num2str(max(abs(fval)))])

subplot 131; pie([(tend-t)/tend, t/tend]); title('% Completed')
subplot 132; semilogy(abs(error)'); title('Solution Residual (mmHg)')
subplot 133; plot(flag~=1)
drawnow

Pt= varT(length(segment)+1 : end);
Pt1=zeros(size(Pt));
Pt1(solvePID)=varT1(length(solveQID)+1 : end);

Qt1=zeros(length(segment),1);
Qt=varT(1:length(segment));
Qt1(solveQID)=varT1(1:length(solveQID));

for i=1:length(InflowID)
    
    flowwave=BCs(BCid(i)).Flow;
    
    Period=flowwave.t(end);
    omega=2*pi/Period;
    tlocal=rem([t], Period);
    Qtemp=interp1(flowwave.t, flowwave.Q, tlocal);
    
    if i==1
    Qt1(InflowID(i))=Qtemp(1);
    else
    Qt1(InflowID(i))= -Qtemp(1);
    end
end

for i=1:length(junctions)
    ID=junctions(i).segs(end);
    
    Qt1(ID)= -junctions(i).defs(end)* ( junctions(i).defs * Qt1(junctions(i).segs) );
end

for i=1:length(BCidout)
    Pidx=length(junctions)+BCidout(i);
    Qidx=OutflowID(i);
    
    Z=BCs(BCidout(i)).RCR(1);
        C=BCs(BCidout(i)).RCR(2);
        R=BCs(BCidout(i)).RCR(3);
        RC=R*C;
    i;
    Q=Qt1(Qidx);
        dQ = ( Qt1(Qidx) - Qt(Qidx) )/dt;
        Pold=Pt(Pidx);
        
    Pt1(Pidx)=( Q*(R+Z) + RC*Z*dQ + RC/dt*Pold ) / (1+RC/dt);
    dP=Pt1(Pidx)-Pold;
end

varT=[Qt1; Pt1];

varsave(:, count)=[t;varT; fval];
error(:, count)=fval;
flag(count)=exitflag;
end

%% Make Table of outputs
VarNames=cell(1, length(segment)+length(junctions)+sum(BC_def==1) +1 + NVar);
VarNames{1}='t';
for i=1:length(segment)
    VarNames{i+1}=['Q' num2str(i)];
end

n=i;
for i=1:length(junctions)
    VarNames{i+1+n}=['Pjcn' num2str(i)];
end
n=i+n;
for i=1:length(BCs)
    VarNames{i+1+n}=['Pseg' num2str(i)];
end
n=i+n;
for i=1:length(segment)
    VarNames{i+1+n}=['ErrorSeg' num2str(i)];
end

VarTable=array2table([varsave]', 'VariableNames', VarNames);

end

