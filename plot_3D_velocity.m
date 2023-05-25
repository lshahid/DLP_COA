%% Plot Pressure

start=212;
fin=284;
it=2;
Q_all=table2array(VarTable(:, 2:1+length(segment)));
maxV=0;

for frame = start:it:fin
    
    for i=1:length(segment)
     A = segment(i).data(:,14);
    A(1:trim_centerline(i,1))=A(trim_centerline(i,1));
    A(trim_centerline(i,2):end)=A(trim_centerline(i,2));
     V=abs(Q_all(frame, i)./A);
    maxV=max([max(V) maxV]);
    end
end


for frame = start:it:fin

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
    
    A = segment(i).data(:,14);
    A(1:trim_centerline(i,1))=A(trim_centerline(i,1));
    A(trim_centerline(i,2):end)=A(trim_centerline(i,2));
    
    V=abs(Q_all(frame, i)./A);
    surf(xx,yy,zz,[V V],'EdgeColor','interp', 'LineWidth', 5, 'EdgeAlpha', 0.8)
end 
   %title([num2str(frame) ' / ' num2str(fin)])
    caxis([0, maxV])
    view([0 0 1])
  colorbar
    if frame==start
        gif('PigC_vel.gif', 'frame', gcf)
    else
        gif
    end

end