function [ trim_centerline ] = clean_centerline( segment, STL )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

trim_centerline=zeros(length(segment),2);


for highlight=1:length(segment)
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
        figure
        plot(segment(i).data(:,14), 'o-'); grid on
    else
        figure(1)
    plot3(segment(i).data(:,1), segment(i).data(:,2),segment(i).data(:,3),'k-')
    end
    end

   trim_centerline(highlight,:)=input('Type start and end points [x1 x2]: ')
end

close all
end



