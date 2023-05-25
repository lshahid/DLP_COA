%% Calc RAC

A=[];
P=[];
L=[];
for seg=1:3
   
    idx = trim_centerline(seg,1):trim_centerline(seg,2);
    Atemp = segment(seg).data(idx,14)./100;
    
    A = [A; Atemp];
    
    if seg==1
        dP0 = range(varsave_w_P(round(0.75*end):end,25));
        dP1 = range(varsave_w_P(round(0.75*end):end,22));
    elseif seg==2
        dP0 = range(varsave_w_P(round(0.75*end):end,22));
        dP1 = range(varsave_w_P(round(0.75*end):end,23));
    else
        dP0 = range(varsave_w_P(round(0.75*end):end,23));
        dP1 = range(varsave_w_P(round(0.75*end):end,26));
    end
    
    Ptemp = linspace(dP0, dP1, length(Atemp) )';
    
    P = [P;Ptemp];
    
    xyz=segment(seg).data(idx,1:3)./10;
    dL=sqrt(sum(diff(xyz).^2,2));
    Ltemp=zeros(size(Atemp));
    
    for j=2:length(Atemp)
        Ltemp(j)=Ltemp(j-1)+dL(j-1);
    end
    
    if seg>1
    Ltemp=Ltemp+max(L)*ones(size(Ltemp));
    end
    L=[L;Ltemp];
    
end

figure; plot(L,A)
figure
plot(L,P)

D = sqrt(4.*A./pi);

dA = A.*D.*P./(E*h);

figure
plot(L,dA)

RAC = 100*dA./A;

figure; plot(10*L,flipud(RAC))
axis([0 120 0 70])