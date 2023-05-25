function [junctions, segment] = define_jcts(segment,BCs, trim_centerline)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Find unique junctions

junctions(1).segs=sort( [1 segment(1).P.connect] );
segment(1).P.junction=1;
n=1;
for it=2:length(segment)
    for j=1:size(segment(it).P, 2)
         temp=0; [it j]
        for k=1:length(junctions)
            if length([it segment(it).P(j).connect ]) == length(junctions(k).segs)
                
    if  sort([it segment(it).P(j).connect ]) == junctions(k).segs 
        temp=1;
        segment(it).P(j).junction=k;
        break
        
    end
            end

        end
        if temp==0
            n=n+1;
        junctions(n).segs=sort( [it segment(it).P(j).connect ]);
        segment(it).P(j).junction=n;
        end
    end
end
    
%% Label order of each segment from inlet

for i=1:length(BCs)
    BC_segs(i)=BCs(i).segs;
    BC_def(i)=BCs(i).def;
end

order=zeros(1,length(segment));

inlet=BC_segs(BC_def==0);
inlet=inlet(1);
order(inlet)=1;
order_old=ones(size(order));
while sum(order==0)>0
    maxorder=find(order==max(order));
    newmax=max(order)+1;
    for i=1:length(maxorder)
    
    for j=1:length(junctions)
    if ismember(maxorder(i), junctions(j).segs)
       
        jctID=j;

        temp=junctions(jctID).segs(junctions(jctID).segs~=inlet);
    
    for k=1:length(temp)
        if order(temp(k))==0
                order(temp(k))=newmax;
            if order_old==order
                display('Stuck')
            else
                display('Not Stuck')
            end

              
    order_old=order;
        end
        
    end
    
    end
     
    end

  
    
end

end


%% Define

for i=1:length(junctions)
    junctions(i).defs=- ones(size(junctions(i).segs));
    SupplyID=find(order(junctions(i).segs)==min(order(junctions(i).segs)));
    junctions(i).defs(SupplyID)=1;
end

%% Extract Tangent Vects

for i=1:length(junctions)
    
    ntemp=zeros(3, length(junctions(i).segs));
    Atemp=zeros(1, length(junctions(i).segs));
    
    for j=1:length(junctions(i).segs)
        
        seg=junctions(i).segs(j);
        
        
        if ismember (seg , BC_segs )
            
            if segment( seg ).P.verts(1) == segment( seg ).data(1,1)   
                ntemp(:,j)= segment(seg).data(trim_centerline(seg,1), 4:6)';  
                Atemp(j)= segment(seg).data(trim_centerline(seg,1), 14);  
            else
                ntemp(:,j)= -segment(seg).data(trim_centerline(seg,2), 4:6)';
                Atemp(j)= segment(seg).data(trim_centerline(seg,2), 14);
            end
        else
            if i==segment(seg).P(1).junction
                Atemp(j)= segment(seg).data(trim_centerline(seg,1), 14);  
            else
                ntemp(:,j)= segment(seg).data(trim_centerline(seg,2), 4:6)';
                Atemp(j)= segment(seg).data(trim_centerline(seg,2), 14);
            end
            
        end
        
          
    end
         junctions(i).n=ntemp;
         junctions(i).A=Atemp;
    end
    

end
