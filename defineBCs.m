function [ BCs ] = defineBCs( segment, STL, trim_centerline )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Define BC segs
BCs(1).segs=[];
BC_segs=[];
for i=1:length(segment)
    if length(segment(i).P)==1
BC_segs(end+1)=i;
BCs(length(BC_segs)).segs=i;
    end
end

%% Define BC inlet/outlet Area

for i=1:length(BCs)
    seg=BC_segs(i);
    if segment( seg ).P.verts(1) == segment( seg ).data(1,1)
        BCs(i).A=segment(seg).data(trim_centerline(seg,2), 14);  
    else
        BCs(i).A=segment(seg).data(trim_centerline(seg,1), 14);  
    end
end

%% Define BC inflow/outflow
n=0;
for highlight=BC_segs
    n=n+1;
figure(1)
hold off
trisurf(STL.faces, STL.vertices(:,1), STL.vertices(:,2),STL.vertices(:,3),...
    'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
hold on
daspect([1 1 1])
for i=1:length(segment)
    if i==highlight
        plot3(segment(i).data(:,1), segment(i).data(:,2),segment(i).data(:,3),'ro-', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
        plot3(segment(i).data(1,1), segment(i).data(1,2),segment(i).data(1,3),'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
    else
        figure(1)
    plot3(segment(i).data(:,1), segment(i).data(:,2),segment(i).data(:,3),'k-')
    end
end

   BCs(n).def=input('Type 0 for inlet, 1 for outlet, 2 for rest are outlets: ');
   BC_def(n)=BCs(n).def;
   if BC_def(n)==2
       for j=n:length(BCs)
       BCs(j).def=1;
       BC_def(j)=1;
       end
       break
   end 
end

%% Read in Inflow BCs
InflowID=BC_segs(BC_def==0);
InflowBCID= find(BC_def==0);

for i=1:length(InflowID)
    [filename, pathname]=uigetfile('*.txt','Select Flow File');
    Flow=readtable([pathname filename]);
    Flow.Q(Flow.Q<0)=0;
    tq=linspace(0, Flow.t(end), 101);
    Qq=pchip(Flow.t, Flow.Q, tq);
    
    BCs(InflowBCID(i)).Flow = table(tq', Qq', 'VariableNames', {'t', 'Q'});
end

%% Read in Outflow BCs

[filename, pathname]=uigetfile('*.txt','Select RCR File');
RC = readtable ([pathname filename]);


flag=input('Impediance from PWV (1/0): ');
if flag==1
    PWV=input('Type PWV (cm/s): ');
    Z=1.06*PWV/( BCs(InflowBCID(1)).A /100 );
else
    ZR_ratio=input('Type Z:R ratio: ');
end


%% Assign Outlets to groups
OutflowID=BC_segs(BC_def==1);
OutflowBCID= find(BC_def==1);

n=0;
for highlight=OutflowID
    n=n+1;
    figure(1)
hold off
trisurf(STL.faces, STL.vertices(:,1), STL.vertices(:,2),STL.vertices(:,3),...
    'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
hold on
daspect([1 1 1])
for i=1:length(segment)
    if i==highlight
        plot3(segment(i).data(:,1), segment(i).data(:,2),segment(i).data(:,3),'ro-', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
        plot3(segment(i).data(1,1), segment(i).data(1,2),segment(i).data(1,3),'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
    else
        figure(1)
    plot3(segment(i).data(:,1), segment(i).data(:,2),segment(i).data(:,3),'k-')
    end
end
    
    BCs(OutflowBCID(n)).group=input('Type RCR Group: ');
    group(n)=BCs(OutflowBCID(n)).group;
    
end

%% Distribute Resistance and compliance by Area
 Atotal=zeros(length(RC.R),1);
for i=1:length(OutflowID)
   
    Atotal(BCs(OutflowBCID(i)).group) = Atotal(BCs(OutflowBCID(i)).group) + ...
     BCs(OutflowBCID(i)).A  ;
      
end

for i=1:length(OutflowID)
    Asplit = BCs(OutflowBCID(i)).A / Atotal(BCs(OutflowBCID(i)).group);
    Rtemp= RC.R(BCs(OutflowBCID(i)).group)/Asplit ;
    Ctemp= RC.C(BCs(OutflowBCID(i)).group)*Asplit ;
    
    if flag==1
        if Z>Rtemp
            1
        end
        BCs(OutflowBCID(i)).RCR = [Z Ctemp Rtemp-Z]   ;
        
    else
        BCs(OutflowBCID(i)).RCR = [ZR_ratio*Rtemp Ctemp (1-ZR_ratio)*Rtemp]   ;
    end
end

end

