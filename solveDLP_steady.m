function [ error ] = solveDLP_steady( varT1, dt, segment, junctions, trim_centerline, mat_props, BCs )
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

Qt1=zeros(length(segment),1);

solveQID=1:length(segment);
for i=1:length(InflowID)
solveQID(solveQID==InflowID(i))=[];
end

for i=1:length(junctions)
    ID=junctions(i).segs(end);
    solveQID(solveQID==ID)=[]; 
end

Qt1(solveQID)=varT1(1:length(solveQID));

for i=1:length(InflowID)
    
    flowwave=BCs(BCid(i)).Flow;
    
    Period=flowwave.t(end);
    omega=2*pi/Period;
    Qtemp=interp1(flowwave.t, flowwave.Q, dt);
    
    if i==1
    Qt1(InflowID(i))=Qtemp(1);
    else
    Qt1(InflowID(i))= -Qtemp(1);
    end
end

%% Enforce Mass Conservation

for i=1:length(junctions)
    ID=junctions(i).segs(end);
    
    Qt1(ID)= -junctions(i).defs(end)* ( junctions(i).defs * Qt1(junctions(i).segs) );
end

%% Define Outlet BCs
Pt1= zeros(length(junctions)+length(BCs), 1);

solvePID=1:length(junctions)+length(BCs);
solvePID(length(junctions)+BCidout)=[];
Pt1(solvePID)=varT1(length(solveQID)+1 : end);

for i=1:length(BCidout)
    Pidx=length(junctions)+BCidout(i);
    Qidx=OutflowID(i);
    
    Z=BCs(BCidout(i)).RCR(1);
        C=BCs(BCidout(i)).RCR(2);
        R=BCs(BCidout(i)).RCR(3);
        RC=R*C;
    
    Q=Qt1(Qidx);
        
    Pt1(Pidx)= Q*(R+Z) ;
end

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
        end
    end

     Res(i,:)= [Rv sum(Rs) 0];
      
end

%% Calculate Junction Resistances

for i=1:length(junctions)
    
   segIDs = junctions(i).segs ( Qt1(junctions(i).segs).*junctions(i).defs' <0 );
   jcnsegIDs = find( Qt1(junctions(i).segs).*junctions(i).defs' <0 ) ;
   
   supplyID= junctions(i).segs ( Qt1(junctions(i).segs).*junctions(i).defs' > 0 );
   supplyjcnID= find ( Qt1(junctions(i).segs).*junctions(i).defs' > 0 );
  
   if length(supplyID)==0
       
   else
   if length(supplyID)>1
        n_supply= junctions(i).n(:,supplyjcnID) * abs(Qt1(supplyID));
        n_supply=n_supply./sqrt(n_supply'*n_supply);
        
        Q_supply=sum(abs(Qt1(supplyID)));

        U_supply= sum( Qt1(supplyID).^2 ./ (0.01*junctions(i).A(supplyjcnID) )') / Q_supply;
        
        A_supply= Q_supply/U_supply;
   else
       n_supply=junctions(i).n(:, supplyjcnID);
       Q_supply=Qt1(supplyID);
       A_supply=junctions(i).A(supplyjcnID)/100;
   end
   
       KE_supply=0.5*density*Q_supply^2/A_supply^2;
       
 
   for j=1:length(segIDs)
       
       nb= junctions(i).n(:, jcnsegIDs(j));
       Ab= junctions(i).A(jcnsegIDs(j)) / 100;
       Qb= abs(Qt1(segIDs(j)));
       
       theta= acos( nb' * n_supply );
       
       phi=3*(pi-theta)/4;
       
       flow_split = Qb/Q_supply;
       area_split = A_supply/Ab;
       
       Res(segIDs(j),3)=KE_supply*(1+flow_split^2*area_split^2-2*flow_split*area_split*cos(phi) ) / Qb ;
       
   end
   end
    
end

%% Calculate Residual (Ax-b)
%  Res
%  Qt1
%  Pt1
%  deltaP
for i=1:length(error)
    
    if i<=length(segment) % Cons momentum each segment
        Rseg=sum(Res(i,:));
        error(i) = ( Rseg*Qt1(i) - deltaP(i) )/1333.2;
    end
end

end

