%% Plot Pressure

start=132;
fin=177;
it=2;
PCWP=16;
P_all=table2array(VarTable(:, 2+length(segment):2+length(segment)+length(junctions)+length(BCs)))./1333.2;

P_all=P_all+PCWP;

for frame = 146%start:it:fin

figure(1)
hold off
trisurf(STL.faces, STL.vertices(:,1), STL.vertices(:,2),STL.vertices(:,3),...
    'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
hold on
colormap jet
daspect([1 1 1])
for i=1:length(segment)
    
    x=segment(i).data(:,1);
    y=segment(i).data(:,2);
    z=segment(i).data(:,3);
    xx=[x x]; yy=[y y]; zz=[z z];
    
    if size(segment(i).P, 2)>1 % interior segment
        
        jcnID(1)=segment(i).P(1).junction;
        jcnID(2)=segment(i).P(2).junction;
        
        segID(1)=find(junctions(jcnID(1)).segs == i);
        segID(2)=find(junctions(jcnID(2)).segs == i);
        
        defs(1)= junctions(jcnID(1)).defs(segID(1)); % 1 supply, -1 collector
        defs(2)= junctions(jcnID(2)).defs(segID(2));
        
        % pressure gradient is defined as collector end - supply end
        
        P = linspace(P_all(frame, jcnID(1)), P_all(frame, jcnID(2)), length(x));
  
    else % segment w BC
        BCid=find( BC_segs == i);
        jcnID=segment(i).P(1).junction;
        
           if segment(i).P.verts==segment(i).data(1,1:3)
               P = linspace(P_all(frame, jcnID), P_all(frame, length(junctions)+BCid), length(x));
           else
               P = linspace(P_all(frame, length(junctions)+BCid),P_all(frame, jcnID), length(x));
           end
    end
    surf(xx,yy,zz,[P' P'],'EdgeColor','interp', 'LineWidth', 5, 'EdgeAlpha', 0.8)
   grid off
   xticklabels([])
   yticklabels([])
   zticklabels([])
   set(gcf, 'Color', 'w')
   set(gca, 'fontsize', 16)
   set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
    
end 
   %title([num2str(frame) ' / ' num2str(fin)])
    %caxis([min(P_all(:)), max(P_all(:))])
    caxis([16, 29.5])
    view([0 1 0])
  colorbar
    if frame==start
        gif('Pressure_K_post.gif', 'frame', gcf)
    else
        gif
    end

end