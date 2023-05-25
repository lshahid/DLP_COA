function [ segment ] = read_centerline( filename, pathname )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
fid=fopen([pathname filename], 'r');

n = 1;
segment=[];
segment_lines=[];
data_lines=[];
tline = fgetl(fid);
while ischar(tline)
  tline = fgetl(fid);
  n = n+1;

if tline==-1
    break
end
  
  if length(tline)>=1
  if tline(1)=='['
      segment_lines(end+1)=n;
      
  elseif tline(1:2)=='Po'
       if length(segment)<length(segment_lines)
           temp=textscan(tline, '%s %s %s %s %s %s' );
           segment(length(segment_lines)).P(1).verts=str2num([temp{2}{1}(2:end-1) ' ' temp{3}{1}(1:end-1) ' ' temp{4}{1}(1:end-1)]);
         segment(length(segment_lines)).P(1).connect=[];
        % segment(length(segment_lines)).data=[];
       else
            temp=textscan(tline, '%s %s %s %s %s %s' );
           segment(length(segment_lines)).P(2).verts=str2num([temp{2}{1}(2:end-1) ' ' temp{3}{1}(1:end-1) ' ' temp{4}{1}(1:end-1)]);
            segment(length(segment_lines)).P(2).connect=[];
      end
      
    elseif tline(1:5)=='    B'
         temp=textscan(tline, '%s %s %n' );
         segment(length(segment_lines)).P(end).connect(end+1)=cell2mat(temp(3));
         
  elseif length(tline)>9
      if tline(1:10)=='        Px'
       data_lines(end+1)=n;
       
       n_seg=0;
       
       while length(tline>1)
           tline = fgetl(fid);
           if length(tline)==0
               break
           elseif tline==-1
               break
           end
           
           temp=textscan(tline, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'TreatAsEmpty', 'n/a');
           n_seg=n_seg+1;
           segment(length(segment_lines)).data(n_seg,:)=cell2mat(temp);
       end
       
      end
  end
  
  end
  
end

    
end

