%% Plot Flow
start=132;
fin=177;
it=2;

Q_all=table2array(VarTable(:, 2:1+length(segment)));

for frame = 146%start:it:fin

figure(1)
hold off
trisurf(STL.faces, STL.vertices(:,1), STL.vertices(:,2),STL.vertices(:,3),...
    'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
hold on
colormap(flipud(bone))
daspect([1 1 1])
for i=1:length(segment)
    
    x=segment(i).data(:,1);
    y=segment(i).data(:,2);
    z=segment(i).data(:,3);
    xx=[x x]; yy=[y y]; zz=[z z];
        
        Q = abs(linspace(Q_all(frame, i), Q_all(frame, i), length(x)));
  
    surf(xx,yy,zz,[Q' Q'],'EdgeColor','interp', 'LineWidth', 5, 'EdgeAlpha', 0.8)
   grid off
   xticklabels([])
   yticklabels([])
   zticklabels([])
   set(gcf, 'Color', 'w')
   set(gca, 'fontsize', 16)
   set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
    
end 
   %title([num2str(frame) ' / ' num2str(fin)])
    caxis([0, max(Q_all(:))])
    %caxis([16, 29.5])
    view([0 1 0])
  colorbar
    if frame==start
        gif('Pressure_K_post_flow.gif', 'frame', gcf)
    else
        gif
    end

end