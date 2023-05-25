function [VarTable] = runDLP_steady(dt, segment, junctions, trim_centerline,mat_props, BCs)
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

varT=zeros(length(segment)+length(junctions)+length(BCs),1);

%options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
options = optimoptions('fsolve','Algorithm','trust-region-dogleg');

varsave=[varT];
error=zeros(NVar);
count=1;

flag=1;

var0=zeros(NVar,1);
    
fx = @(varT1) solveDLP_steady( varT1, dt, segment, junctions, ...
    trim_centerline, mat_props, BCs );

[varT1, fval, exitflag]=fsolve(fx, var0, options);

display(['Max Error is ' num2str(max(abs(fval)))])

Pt1=zeros(length(junctions)+length(BCs), 1);
Pt1(solvePID)=varT1(length(solveQID)+1 : end);

Qt1=zeros(length(segment),1);
Qt1(solveQID)=varT1(1:length(solveQID));

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

varsave=[Qt1; Pt1; fval];

error=fval;
flag(count)=exitflag;


%% Make Table of outputs
VarNames=cell(1, length(segment)+length(junctions)+sum(BC_def==1) + NVar);

for i=1:length(segment)
    VarNames{i}=['Q' num2str(i)];
end

n=i;
for i=1:length(junctions)
    VarNames{i+n}=['Pjcn' num2str(i)];
end
n=i+n;
for i=1:length(BCs)
    VarNames{i+n}=['Pseg' num2str(i)];
end
n=i+n;
for i=1:length(segment)
    VarNames{i+n}=['ErrorSeg' num2str(i)];
end

VarTable=array2table(varsave', 'VariableNames', VarNames);

end

