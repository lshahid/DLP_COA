
options = optimoptions('fsolve','Algorithm','levenberg-marquardt');

%%

for i=1:length(BCs)
    BC_segs(i)=BCs(i).segs;
    BC_def(i)=BCs(i).def;
end

BCid=find(BC_def==0);
InflowID=BC_segs(BC_def==0);
BCidout=find(BC_def==1);
OutflowID=BC_segs(BC_def==1);

%%

varT1=zeros(length(segment)+length(BCidout),1);
varT=zeros(length(segment)+length(junctions)+2*length(BCs),1);

varsave=[];
error=[];
flag=[];
close all;
clc
dt=1e-2;
tend=2.44;
count=0;
for t=dt:dt:tend
count=count+1;
    var0=varT1;
fx=@(varT1) solveDLP_FSI(varT1, varT, t, dt, segment, junctions, trim_centerline,...
    mat_props, BCs, Pref);

[varT1, fval, exitflag]=fsolve(fx, var0, options);

display(['t=' num2str(t) ' Max Error is ' num2str(max(abs(fval)))])

subplot 131; pie([(tend-t)/tend, t/tend]); title('% Completed')
subplot 132; semilogy(abs(error)'); title('Solution Residual (mmHg)')
subplot 133; plot(flag~=1)
drawnow

Qt1=zeros(length(segment)+length(BCs),1);
Qt=varT(1:length(Qt1));

Vt= varT(length(Qt)+1 : end);
Vt1= zeros(length(Vt), 1);

solveQID=1:length(segment)+length(BCs);
for i=1:length(BCid)
solveQID(solveQID==BCid(i)+length(segment))=[];
end

Qt1(solveQID)=varT1(1:length(solveQID));

for i=1:length(InflowID)
    
    flowwave=BCs(BCid(i)).Flow;
    
    Period=flowwave.t(end);
    omega=2*pi/Period;
    tlocal=rem([t t-dt], Period);
    Qtemp=interp1(flowwave.t, flowwave.Q, tlocal);
    
    if i==1
    Qt1(BCid(i)+length(segment))=Qtemp(1);
    else
    Qt1(BCid(i)+length(segment))= -Qtemp(1);
    end
end
% Enforce Mass Conservation

for i=1:length(junctions)
    ID=junctions(i).segs(end);
    
    Vt1(i) = dt*( junctions(i).defs * Qt1(junctions(i).segs) ) + Vt(i) ;
end

for i=1:length(InflowID)
   
    Vt1(length(junctions)+BCid(i))=dt*(Qt1(BCid(i)+length(segment))-Qt1(InflowID(i)) ) + Vt(length(junctions)+BCid(i));

end

for i=1:length(BCidout)
    
    Vt1(length(junctions)+BCidout(i)) = dt*(Qt1(OutflowID(i))-Qt1(BCidout(i)+length(segment))) + Vt(length(junctions)+BCidout(i));

end

varT=[Qt1; Vt1];

varsave(:, count)=[t;varT; fval];
error(:, count)=fval;
flag(count)=exitflag;

end

%% Calculate Compliance

C=zeros(length(junctions)+length(BCs),1);
segment=segment
for i=1:length(junctions)
   for j=1:length(junctions(i).segs)
       seg=junctions(i).segs(j);
        
        if ismember (seg , BC_segs )
            
            if segment( seg ).P.verts(1) == segment( seg ).data(1,1)   
                idx = trim_centerline(seg,1):trim_centerline(seg,2);  
            else
               idx = trim_centerline(seg,2):-1:trim_centerline(seg,1);  
            end
        else
            if i==segment(seg).P(1).junction
                idx = trim_centerline(seg,1):trim_centerline(seg,2);  
            else
                idx = trim_centerline(seg,2):-1:trim_centerline(seg,1);  
            end 
        end
        
    A=segment(seg).data(idx,14)./100;
   
    R=sqrt(A./pi);
    
    N=linspace(1,0,length(idx));
    
%     minID=find(islocalmin(A));
%     
%     Rs=0;
%     
%     for j=1:length(minID)
%         if minID~=1 | minID~=length(R)
%             As=A(minID(j));
%             if length(minID)==1
%                 A0=mean([ max(A(1:minID(j))), max(A(minID(j):end)) ]);
%              elseif j==1
%                 A0=mean( [max(A(1:minID(j))), max(A(minID(j): minID(j+1)))] );
%             elseif j==length(minID)
%                 A0=mean([ max(A(minID(j-1):minID(j))), max(A(minID(j):end)) ]);
%             else
%                 A0=mean([ max(A(minID(j-1):minID(j))), max(A(minID(j):minID(j+1))) ]);
%             end
%             
%             Rs(j)= density*1.52 /(2* A0^2)*(A0/As - 1)^2 * abs(Qt1(i));
%              if As/A0 <0.5
%                  N(1:minID(j))=1;
%                  N(minID(J:end))=0;
%              end
%         end
%     end
%     
    xyz=segment(seg).data(idx,1:3)./10;
    dL=sqrt(sum(diff(xyz).^2,2));
    L=zeros(size(R));
    
    for j=2:length(R)
        L(j)=L(j-1)+dL(j-1);
    end
    
    Eh=segment(seg).E*segment(seg).h;
    
    C(i) = trapz(L, N'.*( 2.*A.*R./Eh) ) + C(i);
        
   end
end

for i=1:length(BCs)
           seg=BCs(i).segs;

           if length(segment)==1
                idx=trim_centerline(1):trim_centerline(2);
                if i==2; idx=fliplr(idx); end
           else
           
            if segment( seg ).P.verts(1) == segment( seg ).data(1,1)   
                idx = trim_centerline(seg,2):-1:trim_centerline(seg,1);  
            else
               idx = trim_centerline(seg,1):trim_centerline(seg,2);  
            end
           end
    A=segment(seg).data(idx,14)./100;
   
    R=sqrt(A./pi);
    
    N=linspace(1,0,length(idx));
    
%     minID=find(islocalmin(A));
%     
%     Rs=0;
%     
%     for j=1:length(minID)
%         if minID~=1 | minID~=length(R)
%             As=A(minID(j));
%             if length(minID)==1
%                 A0=mean([ max(A(1:minID(j))), max(A(minID(j):end)) ]);
%              elseif j==1
%                 A0=mean( [max(A(1:minID(j))), max(A(minID(j): minID(j+1)))] );
%             elseif j==length(minID)
%                 A0=mean([ max(A(minID(j-1):minID(j))), max(A(minID(j):end)) ]);
%             else
%                 A0=mean([ max(A(minID(j-1):minID(j))), max(A(minID(j):minID(j+1))) ]);
%             end
%             
%             Rs(j)= density*1.52 /(2* A0^2)*(A0/As - 1)^2 * abs(Qt1(i));
%              if As/A0 <0.5
%                  N(1:minID(j))=1;
%                  N(minID(J:end))=0;
%              end
%         end
%     end
%     
    xyz=segment(seg).data(idx,1:3)./10;
    dL=sqrt(sum(diff(xyz).^2,2));
    L=zeros(size(R));
    
    for j=2:length(R)
        L(j)=L(j-1)+dL(j-1);
    end
    
    Eh=segment(seg).E*segment(seg).h;
    
    C(i+length(junctions)) = trapz(L, N'.*( 2.*A.*R./Eh) );

end


%% Calc Pressure
V = varsave(length(Qt)+2: end-length(fval), : );
P=V./C+Pref;
plot(P'./1333.2)

mP=mean(P(5,round(0.75*end):end)./1333.2)
PP = range(P(5,round(0.75*end):end)./1333.2)
dP = min(P(5,round(0.75*end):end)./1333.2)
sP = max(P(5,round(0.75*end):end)./1333.2)
%%
varsave_w_P=[varsave(1:end-length(fval),:)' P' varsave(end+1-length(fval) : end,:)'];

%idx=find(varsave_w_P(:, 1)>=9);
%cp = varsave_w_P(idx(1:10:end), [3 4 7 8]);
%% Calc PWV

% Qin = varsave_w_P(1:70,3);
% Qout = varsave_w_P(1:200,4);
% 
% Qin_scale = (Qin-min(Qin))./range(Qin);
% Qout_scale = (Qout-min(Qout))./range(Qout);
% 
% plot(Qin_scale); hold on; plot(Qout_scale);
% 
% t=dt:dt:1;
% 
% idx_in=find((Qin_scale<0.8) & (Qin_scale>0.2));
% idx_in = idx_in(idx_in<50);
% pin=polyfit(t(idx_in), Qin_scale(idx_in), 1);
% tUin = -pin(2)/pin(1);
% tMin = (0.5-pin(2))/pin(1);
% 
% idx_out=find((Qout_scale<0.8) & (Qout_scale>0.2));
% idx_out = idx_out(idx_out<200);
% pout=polyfit(t(idx_out), Qout_scale(idx_out), 1);
% tUout = -pout(2)/pout(1);
% tMout = (0.5-pout(2))/pout(1);
% 
% PWVU=0.3/(tUout-tUin)
% PWVM=0.3/(tMout-tMin)