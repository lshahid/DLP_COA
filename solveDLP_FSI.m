function [ error ] = solveDLP( varT1, varT, t, dt, segment, junctions, trim_centerline, mat_props, BCs, Pref)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

error=zeros(length(varT1),1);

%% Define Material Properties
density=mat_props.density;
visc=mat_props.visc;

%% Define dirichlet inflow BCs

for i=1:length(BCs)
    BC_segs(i)=BCs(i).segs;
    BC_def(i)=BCs(i).def;
end

BCid=find(BC_def==0);
InflowID=BC_segs(BC_def==0);
BCidout=find(BC_def==1);
OutflowID=BC_segs(BC_def==1);

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
%% Enforce Mass Conservation

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

%% Calculate Compliance

C=zeros(length(junctions)+length(BCs),1);

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

Pt1=Vt1./C+Pref;
Pt=Vt./C+Pref;
%% Calculate deltaP for each segment

deltaP=zeros(length(segment),1);

for i=1:length(segment)
    
    if length(segment)==1
        deltaP(i)= Pt1(1) - Pt1(2); 
        
    elseif size(segment(i).P, 2)>1 % interior segment
        
        jcnID(1)=segment(i).P(1).junction;
        jcnID(2)=segment(i).P(2).junction;
        
        segID(1)=find(junctions(jcnID(1)).segs == i);
        segID(2)=find(junctions(jcnID(2)).segs == i);
        
        defs(1)= junctions(jcnID(1)).defs(segID(1)); % 1 supply, -1 collector
        defs(2)= junctions(jcnID(2)).defs(segID(2));
        
        % pressure gradient is defined as collector end - supply end
        
        deltaP(i) = -defs * [Pt1(jcnID(1)); Pt1(jcnID(2))];
  
    else % segment w BC
        BCid=find( BC_segs == i);
        jcnID=segment(i).P(1).junction;
        if i==InflowID(1) % if primary inflow
            deltaP(i)= Pt1(length(junctions)+BCid) - Pt1(jcnID) ;
        else % rest of BCs
            deltaP(i)= Pt1(jcnID)- Pt1(length(junctions)+BCid) ; 
        end
    end 
end

%% Calculate segment resistances
flag_PAS=0;
for i=1:length(segment)
    
    idx=trim_centerline(i,1):trim_centerline(i,2);
    A=segment(i).data(idx,14)./100;
   
    R=sqrt(A./pi);
    
    curvature=segment(i).data(idx,13).*10;
    
    xyz=segment(i).data(idx,1:3)./10;
    dL=sqrt(sum(diff(xyz).^2,2));
    L=zeros(size(R));
    for j=2:length(R)
        L(j)=L(j-1)+dL(j-1);
    end
    
    inertia(i)=(density/pi) * trapz(L, 1./R.^2);
    
    Re=((Qt1(i)./A)*density*2.*R)/visc;
    
    K=Re.*sqrt(R.*curvature);
    
    gamma=0.1033.*sqrt(K).*(sqrt(1+1.729./K)-1.315./sqrt(K)).^(-3);
    gamma(K<14)=1;

    alpha= R.*sqrt(density*omega/visc);
    Alpha=alpha.*(sqrt(-1)-1)/sqrt(2);
    zeta=abs( 1+(2/Alpha) * besselj(1,Alpha)./besselj(0,Alpha) );
    
    Rv= ( 8*visc/pi ) * trapz(L, max(gamma, zeta)./(R.^4)) ;
    %Rv= ( 8*visc/pi ) * trapz(L, 1./(R.^4)) ;
    minID=find(islocalmin(A));
    
    Rs=0;
    
    for j=1:length(minID)
        if minID~=1 | minID~=length(R)
            As=A(minID(j));
            if length(minID)==1
                A0=mean([ max(A(1:minID(j))), max(A(minID(j):end)) ]);
             elseif j==1
                A0=mean( [max(A(1:minID(j))), max(A(minID(j): minID(j+1)))] );
            elseif j==length(minID)
                A0=mean([ max(A(minID(j-1):minID(j))), max(A(minID(j):end)) ]);
            else
                A0=mean([ max(A(minID(j-1):minID(j))), max(A(minID(j):minID(j+1))) ]);
            end
            
            Rs(j)= density*1.52 /(2* A0^2)*(A0/As - 1)^2 * abs(Qt1(i));
%             if As/A0 <0.5
%                 flag_PAS(end+1)=i;
%             end
        end
    end

     Res(i,:)= [Rv sum(Rs) 0];
      
end

%% Calculate Junction Resistances

% for i=1:length(junctions)
%     
%    segIDs = junctions(i).segs ( Qt1(junctions(i).segs).*junctions(i).defs' <0 );
%    jcnsegIDs = find( Qt1(junctions(i).segs).*junctions(i).defs' <0 ) ;
%    
%    supplyID= junctions(i).segs ( Qt1(junctions(i).segs).*junctions(i).defs' > 0 );
%    supplyjcnID= find ( Qt1(junctions(i).segs).*junctions(i).defs' > 0 );
%   
%    if length(supplyID)==0
%        
%    else
%    if length(supplyID)>1
%         n_supply= junctions(i).n(:,supplyjcnID) * abs(Qt1(supplyID));
%         n_supply=n_supply./sqrt(n_supply'*n_supply);
%         
%         Q_supply=sum(abs(Qt1(supplyID)));
% 
%         U_supply= sum( Qt1(supplyID).^2 ./ (0.01*junctions(i).A(supplyjcnID) )') / Q_supply;
%         
%         A_supply= Q_supply/U_supply;
%    else
%        n_supply=junctions(i).n(:, supplyjcnID);
%        Q_supply=Qt1(supplyID);
%        A_supply=junctions(i).A(supplyjcnID)/100;
%    end
%    
%        KE_supply=0.5*density*Q_supply^2/A_supply^2;
%        
%  
%    for j=1:length(segIDs)
%        
%        nb= junctions(i).n(:, jcnsegIDs(j));
%        Ab= junctions(i).A(jcnsegIDs(j)) / 100;
%        Qb= abs(Qt1(segIDs(j)));
%        
%        theta= acos( nb' * n_supply );
%        
%        phi=3*(pi-theta)/4;
%        
%        flow_split = Qb/Q_supply;
%        area_split = A_supply/Ab;
%        
%        Res(segIDs(j),3)=KE_supply*(1+flow_split^2*area_split^2-2*flow_split*area_split*cos(phi) ) / Qb ;
%        
%    end
%    end
%     
% end
% 
%% Calculate Residual (Ax-b)
%  Res
%  Qt1
%  Qt
%  Pt1
%  Pt
%  deltaP
for i=1:length(error)
    
    if i<=length(segment)
   % if ismember(i, flag_PAS)
       Rseg=sum(Res(i,[1 2]));
%       Rseg=0;
%     else
%         Rseg=sum(Res(i,:));
%     end
        
        error(i) = ( inertia(i)*(Qt1(i)-Qt(i))/dt + Rseg*Qt1(i) - deltaP(i) )/1333.2;
        %error(i) = Rseg*Qt1(i) - deltaP(i);

    else
            j=i-length(segment);
            Pidx=length(junctions)+BCidout(j);
            Qidx=length(segment)+BCidout(j);
            Z=BCs(BCidout(j)).RCR(1);
            C_RCR=BCs(BCidout(j)).RCR(2);
            R=BCs(BCidout(j)).RCR(3);
            RC=R*C_RCR;
    
            Q=Qt1(Qidx);
            dQ = ( Qt1(Qidx) - Qt(Qidx) )/dt;
            Pold=Pt(Pidx);
        
            Ptemp=( Q*(R+Z) + RC*Z*dQ + RC/dt*Pold ) / (1+RC/dt);
            
            error(i)=Pt1(Pidx)-Ptemp;
        end
        
end


end

