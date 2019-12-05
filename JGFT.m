function [GSig,GFTSig,X,TV,TVn] = JGFT(x,L,w_nwind, anWin, w_noverlap)
[V,D] = eig(L);
[n,m] = size(x);
FreqResol = w_nwind;
numberOfFrame = floor(n/w_nwind*2)-1;
X = zeros(FreqResol,numberOfFrame,m);
for i = 1:m
    xx = enframe(x(:,i), anWin, w_noverlap)';
    for k = 1:numberOfFrame
        X(:,k,i) = fft(xx(:,k),w_nwind);
    end
end
GSig = zeros(m,numberOfFrame,FreqResol);
GFTSig = zeros(m,numberOfFrame,FreqResol);
for l = 1:FreqResol
    for k= 1:numberOfFrame
        tmp = permute(X(l,k,:),[3,1,2]);
        GSig(:,k,l) = tmp;
        GFTSig(:,k,l) = V'*tmp;
    end
end
TV = zeros(numberOfFrame,FreqResol);
TVn = zeros(numberOfFrame,FreqResol);
for l = 1:FreqResol
    for k= 1:numberOfFrame
           TVn(k,l) = GFTSig(:,k,l)'* L *GFTSig(:,k,l)/(GFTSig(:,k,l)'*GFTSig(:,k,l));
           TV(k,l) = abs(GSig(:,k,l)')*L*abs(GSig(:,k,l));
    end
end