%% Example of graph-based atrial activity extraction

clear;
clc;
%% -----------------construct graph---------------------%
load 'coord.mat';
W = mygraph(coord);
G= gsp_graph(W);
G.coords = coord;
param = struct('show_edges',1);
param.vertex_highlight = 10;
gsp_plot_graph(G,param);
title('Graph Structure')
set(gca, 'fontsize', 16)

D = diag(sum(W-diag(diag(W))));
L = D - W;

%% -----------------load data---------------------%
load egm1.mat
x_AF = egm1;
load pureaa1.mat;
pureAA_AF = pureaa1;
load pureva1.mat;
pureVA_AF = pureva1;

[m,n] = size(x_AF);
fs = 1000;
t = (0:n-1)/1000;


%% -----------------window-----------------
w_nwind = fix(0.1 * fs);  % 0.1s long
w_noverlap = fix(0.5 * w_nwind);  % 50% overlapping
anWin = hanning(w_nwind); % hanning window
FreqResol = w_nwind;
numberOfFrame = floor(n/w_nwind*2)-1;
ff = (0:FreqResol-1) * fs / FreqResol;

%% ---------------------Find peak of ventricular signal AF signals------------------%
[L_AF,LC_AF,locs_AF] = findECGPeak_stft(pureVA_AF(1,:), w_nwind,'AF');
z_AF = zeros(size(L_AF,1),1);
zc_AF = zeros(size(LC_AF,1),1);
y_AF = zeros(size(locs_AF,1),1)-1.2;

%% ---------------------JGFT and Total variance---------------------
[V,Gamma] = eig(L);
x = x_AF';
[GSig_AF,GFTSig_AF,X_AF,TV_AF,TVn_AF] = JGFT(x,L,w_nwind, anWin, w_noverlap);

%% ------------------detect the presence of ventricular signal---------------------
GFrm_AF = [];
for i = 1:numberOfFrame
    if abs(GFTSig_AF(1,i,3))>30
        GFrm_AF = [GFrm_AF, i];  
    end
end

%% ------------------Extract the atrial signal---------------------
c = 0.14;
lambda = 0.8;
[xp_AF] = Extr_atrial(x_AF,GSig_AF,X_AF,L,w_nwind,w_noverlap,GFrm_AF,lambda,c);


%% ---------------plot pure, mixed and extracted signals----------------
figure;
subplot 411;
plot(locs_AF/1000,y_AF,'or','markersize',2.5) 
hold on;
plot(t,pureAA_AF(1,:),'LineWidth',1.5);
title('Atrial activity');
set(gca, 'fontsize', 16)
axis([0,max(t),-1.2,1.2])
xlabel('time (s)')
subplot 412;
plot(locs_AF/1000,y_AF,'or','markersize',2.5) 
hold on;
plot(t,pureVA_AF(1,:),'LineWidth',1.5);
title('Ventricular activity')
set(gca, 'fontsize', 16)
axis([0,max(t),-1.2,1.2])
xlabel('time (s)')
subplot 413;
plot(locs_AF/1000,y_AF,'or','markersize',2.5) 
hold on;
plot(t,x_AF(1,:),'LineWidth',1.5);
title('Mixed EGM')
set(gca, 'fontsize', 16)
axis([0,max(t),-1.2,1.2])
xlabel('time (s)')
subplot 414;
plot(locs_AF/1000,y_AF,'or','markersize',2.5) 
hold on;
plot(t,xp_AF(1,:),'LineWidth',1.5);
axis([0,max(t),-1.2,1.2])
title('Extracted atrial activity')
xlabel('time (s)')
set(gca, 'fontsize', 16)



