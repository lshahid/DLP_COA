%% Distributed Lumped Parameter Hemodynamic Model

%% Read Centerline File

[filename, pathname]=uigetfile('*.txt', 'Select Mimics Centerline File');

[filename_STL, pathname_STL]=uigetfile('*.stl', 'Select Mimics Centerline File');

STL=stlread([pathname_STL filename_STL]);

[segment]=read_centerline(filename, pathname);

%%

trim_centerline=clean_centerline(segment, STL);

BCs=defineBCs(segment, STL, trim_centerline);

[junctions, segment] = define_jcts(segment,BCs, trim_centerline);

savename=input('Type file save name', 's');
save(savename);

%% 

tic
clc

dt=0.01;
tend=4*BCs(1).Flow.t(end);

mat_props.visc=0.035;
mat_props.density=1.06;

steady_flag=0; % Transient = 0, Steady = 1

if steady_flag==0
    
[VarTable] = runDLP(dt, tend, segment, junctions, trim_centerline,mat_props, BCs);

elseif steady_flag==1
[VarTable] = runDLP_steady(0.1, segment, junctions, trim_centerline,mat_props, BCs);

end
run_time=toc

save TJS_post_output.mat


%% skip this for coarctation
tic
clc

dt=0.01;
tend=4*BCs(2).Flow.t(end);

mat_props.visc=0.036;
mat_props.density=1.06;

steady_flag=0;

if steady_flag==0
    
[VarTable] = runDLP(dt, tend, segment_stented, junctions, trim_centerline,mat_props, BC_RHC);

elseif steady_flag==1
[VarTable] = runDLP_steady(0.1, segment, junctions, trim_centerline,mat_props, BCs);

end
run_time=toc

save PigW10_RHC_post.mat
%% skip this for coarc
period=0.5;
omega=2*pi/period;
Res=[];
for i=1:length(t)
    
    Qt1=table2array(VarTable(i, 2:length(segment2)+1))';
    [Rtemp] = calcRpost(Qt1, segment2, junctions, trim_centerline, mat_props, omega);
    Rtemp=Rtemp';
    
    Res(i,:)=Rtemp(:)';
end

%% Plot Results
BC_segs=[];
BC_def=[];
for i=1:length(BCs)
    BC_segs(i)=BCs(i).segs;
    BC_def(i)=BCs(i).def;
end
BCidout=find(BC_def==1);


group=zeros(size(BCidout));

for i=1:length(BCidout)
    group(i) =  BCs(BCidout(i)).group;
end

Qout=zeros(length(VarTable.t),max(group));
for i=1:max(group)
    BCgroup=BCidout(group==i);
    idx=zeros(size(BCgroup));
    for j=1:sum(group==i)
    idx(j)=1+BCs(BCgroup(j)).segs;
    end
    Qout(:,i) = sum( table2array(VarTable(:,idx)), 2 );
end

LPBFp=100*sum(Qout(round(end*.75):end,2))/( sum(Qout(round(end*.75):end,1)) + sum(Qout(round(end*.75):end,2)) ) 
%%

VarNames=VarTable.Properties.VariableNames;

% Flow
idx=2:1+length(segment);
%idx=[2 3 12]
figure
plot(VarTable.t, table2array(VarTable(:,idx)))
legend(VarNames(idx))
title('Flow Rate (mL/s)')
xlabel('Time (sec)')

% Pressure
idx=2+length(segment):1+length(segment)+length(junctions)+length(BCs) ;

%idx=[67 47 53]; %seg1 j2 j8
figure
plot(VarTable.t, abs(table2array(VarTable(:,idx))./1333.2))
legend(VarNames(idx))
title('Pressure (mmHg)')
xlabel('Time (sec)')

%%

idx=[87 68]; %seg2 j10
P=table2array(VarTable(round(end*.75):end,idx))./1333.2;
figure; plot(P);

MPAs=max(P(:,1));
MPAd=min(P(:,1));

LPAs=max(P(:,2));
LPAd=min(P(:,2));

Pgrad=MPAs-LPAs;

BC_segs=[];
BC_def=[];
for i=1:length(BCs)
    BC_segs(i)=BCs(i).segs;
    BC_def(i)=BCs(i).def;
end
BCidout=find(BC_def==1);


group=zeros(size(BCidout));

for i=1:length(BCidout)
    group(i) =  BCs(BCidout(i)).group;
end

Qout=zeros(length(VarTable.t),max(group));
for i=1:max(group)
    BCgroup=BCidout(group==i);
    idx=zeros(size(BCgroup));
    for j=1:sum(group==i)
    idx(j)=1+BCs(BCgroup(j)).segs;
    end
    Qout(:,i) = sum( table2array(VarTable(:,idx)), 2 );
end

figure;plot(Qout);
LPBFp=100*sum(Qout(round(end*.75):end,2))/( sum(Qout(round(end*.75):end,1)) + sum(Qout(round(end*.75):end,2)) ) ;


cp_data=[MPAs MPAd LPAs LPAd Pgrad LPBFp]