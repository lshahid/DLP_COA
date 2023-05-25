function [Res] = calcRpost(Qt1, segment, junctions, trim_centerline, mat_props, omega)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

density=mat_props.density;
visc=mat_props.visc;

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

Res(:,4)=sum(Res, 2);

end

