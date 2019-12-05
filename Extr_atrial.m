function [xp] = Extr_atrial(x, GSig,X,L,w_nwind,w_noverlap,loc_v,lambda,c)
[m,n] = size(x);
FreqResol = w_nwind;
numberOfFrame = floor(n/w_nwind*2)-1;
Xp = zeros(size(X));
GSig1=GSig;
[U,Gamma] = eig(L);
for l = 1:FreqResol
    for k= 1:numberOfFrame
         if find(loc_v==k)~=0
            X1 = U'*GSig(:,k,l);                    
            Xest=((1-lambda*c) * eye(m) + lambda*Gamma) \ X1;
            Aest = X1-Xest;
            aa_est = U*Aest;
            GSig1(:,k,l) =aa_est;
         end
        tmp1 = permute(GSig1(:,k,l),[2,3,1]);
        Xp(l,k,:) = tmp1;
    end
end
xp = [];
for i = 1:m
    Xp1 = Xp(:,:,i);
    output = OverlapAdd2(abs(Xp1),angle(Xp1),w_nwind,w_noverlap);
    ol = length(output);
    if ol<n
       output = [output;x(i,ol+1:n)'];
    end
    xp =[xp,output];
end
xp = xp';
