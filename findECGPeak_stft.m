function [L, LC,locs] = findECGPeak_stft1(ECG, w_nwind,type)
[n,m] = size(ECG);
numberOfFrame = floor(n/w_nwind*2)-1;
ecg = abs(ECG).^2;
if type == 'SR'
    [locs, v] = findpeaks(ecg, [], 600);
else
    [locs, v] = findpeaks(ecg, [], 300);
end
LC = ceil((locs/(w_nwind/2)));
numP=size(LC,1);
numInter = numP*2;
L = zeros(numP+numInter,1);
for i = 1:numP
    L((i-1)*3+1,1) = LC(i)-1;
    L((i-1)*3+2,1) = LC(i);
    L((i-1)*3+3,1) = LC(i)+1;
end
if L(1) == 0 
    L(1) = [];
end
if L(end) == numberOfFrame+1
    L(end) = [];
end
end