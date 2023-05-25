%% Place Stent
close all
stent_loc=11; % segment # that needs to be stented

segment_stented=segment;

idx=trim_centerline(stent_loc(1),1):trim_centerline(stent_loc(1),2);
A=segment_stented(stent_loc).data(:,14);

for highlight=stent_loc
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

end

figure(2)
plot(A, 'o-')

stent_end=input('Type Stent end: ');
stent_start=input('Type Stent start: ');

if stent_end>stent_start
    stent_idx=stent_start:stent_end;
else
    stent_idx=stent_end:stent_start;
end

segment_stented(stent_loc).data(stent_idx,14)=A(stent_start);

figure(2)
hold on
plot(segment_stented(stent_loc).data(:,14), 'o-')

segment_stented(stent_loc).data(stent_idx,13)=0;

segment_stented(stent_loc).data(stent_idx,1)=linspace( segment_stented(stent_loc).data(stent_idx(1),1), ...
    segment_stented(stent_loc).data(stent_idx(end),1), length(stent_idx));

segment_stented(stent_loc).data(stent_idx,2)=linspace( segment_stented(stent_loc).data(stent_idx(1),2), ...
    segment_stented(stent_loc).data(stent_idx(end),2), length(stent_idx));

segment_stented(stent_loc).data(stent_idx,3)=linspace( segment_stented(stent_loc).data(stent_idx(1),3), ...
    segment_stented(stent_loc).data(stent_idx(end),3), length(stent_idx));
