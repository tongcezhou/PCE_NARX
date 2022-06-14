%% Distribution fitting
%==========================================================================
clear,clc,pack
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\';
cd(local_dir)

load('ElCentro.mat');
EQ = ElCentro(:,2); 
N = length(EQ);
time = 1:N;
save('validEQ.txt','EQ','-ascii')

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

load('PEER_Nleq10000_ZETAinterp.mat','RSSref','zeta_ref');
minRSS = min(RSSref,[],2);
RSSindx_IRF = find(minRSS > 0 & minRSS < 1.93);
minRSS = minRSS(minRSS > 0 & minRSS < 1.93);
load('PEER_ALL_MODopt_Arias.mat','RSS','Ia','D545s','D595s');
RSSindx_MOD = find(RSS < 1.93);
EQindx = intersect(RSSindx_IRF,RSSindx_MOD);
load('PEER_ALL_IRFinterp.mat');
wmid = theta(EQindx,1);
wdot = theta(EQindx,2);
zeta = zeta_ref(EQindx);
Ia = sqrt(Ia(EQindx))';
D545s = D545s(EQindx)';
D595s = D595s(EQindx)';


% omega_mid
%--------------------------------------------------------------------------
xindx = 0:0.001:20;
[phat] = mle(wmid,'distribution','loglogistic');
PD_wmid = ProbDistUnivParam('loglogistic',phat);

figure(1)
set(gca,'Fontname','TimesNewRoman')
subplot(9,3,[1 4 7]),set(gca,'Fontname','TimesNewRoman','Fontsize',8)
histfit(wmid,15,'loglogistic')
h = get(gca,'Children');
set(h(2),'FaceColor',[.6 .6 0.8])
xlabel('$\omega_{\mbox{mid}}$ (Hz)','Fontname','TimesNewRoman','Interpreter','Latex','Fontsize',11)
ylabel('Frequency','Fontname','TimesNewRoman','Fontsize',11)
text(10,370,'Log-Logistic(1.64,0.238)','Fontname','TimesNewRoman','Fontsize',11,'HorizontalAlignment','Center')
xlim([0 20])


% omega'
%--------------------------------------------------------------------------
xindx = -0.4:0.001:0.4;
[phat] = mle(wdot,'distribution','logistic');
PD_wdot = ProbDistUnivParam('logistic',phat);

subplot(9,3,[2 5 8]),set(gca,'Fontname','TimesNewRoman','Fontsize',8)
histfit(wdot,15,'logistic')
h = get(gca,'Children');
set(h(2),'FaceColor',[.6 .6 0.8])
xlabel('$\omega''$ (Hz/s)','Fontname','TimesNewRoman','Interpreter','Latex','Fontsize',11)
ylabel('Frequency','Fontname','TimesNewRoman','Fontsize',11)
text(0,840,'Logistic(-0.0516,0.0336)','Fontname','TimesNewRoman','Fontsize',11,'HorizontalAlignment','Center') 
axis([-0.5 0.5 0 900])

% zeta
%--------------------------------------------------------------------------
xindx = 0:0.01:1;
[phat] = mle(zeta,'distribution','beta');
PD_zeta = ProbDistUnivParam('beta',phat);

subplot(9,3,[3 6 9]),set(gca,'Fontname','TimesNewRoman','Fontsize',8)
histfit(zeta,15,'beta')
h = get(gca,'Children');
set(h(2),'FaceColor',[.6 .6 0.8])
xlabel('$\zeta_f$','Fontname','TimesNewRoman','Interpreter','Latex','Fontsize',11)
ylabel('Frequency','Fontname','TimesNewRoman','Fontsize',11)
text(0.5,235,'Beta(1.38,3.70)','Fontname','TimesNewRoman','Fontsize',11,'HorizontalAlignment','Center') 
xlim([0 1])


% Ia
%--------------------------------------------------------------------------
xindx = [0:1e-4:0.05]';
[phat] = mle(Ia,'distribution','exponential');
PD_Ia = ProbDistUnivParam('exponential',phat);

subplot(9,3,[13 16 19]),set(gca,'Fontname','TimesNewRoman','Fontsize',8)
histfit(Ia,15,'exponential')
h = get(gca,'Children');
set(h(2),'FaceColor',[.6 .6 0.8])
xlabel('$\sqrt{I_a}$','Interpreter','Latex','Fontsize',11)
ylabel('Frequency','Fontname','TimesNewRoman','Fontsize',11)
text(0.025,750,'Exp(0.00551)','Fontname','TimesNewRoman','Fontsize',11,'HorizontalAlignment','Center') 
xlim([0 0.05])


% D5-45
%--------------------------------------------------------------------------
Nbins = 1000;
xindx = linspace(0,0.5,Nbins)';
options = statset('Display','final','Maxiter',1000);
[phat] = mle(D545s,'distribution','beta');
PD_D545s = ProbDistUnivParam('beta',phat);

subplot(9,3,[14 17 20]),set(gca,'Fontname','TimesNewRoman','Fontsize',8)
histfit(D545s,15,'beta')
h = get(gca,'Children');
set(h(2),'FaceColor',[.6 .6 0.8])
xlabel('$t_{\mbox{mid}}$','Interpreter','Latex','Fontsize',11)
ylabel('Frequency','Fontname','TimesNewRoman','Fontsize',11)
text(0.25,150,'Beta(3.79,24.23)','Fontname','TimesNewRoman','Fontsize',11,'HorizontalAlignment','Center') 
axis([min(xindx) max(xindx) 0 160])


% D5-95
%--------------------------------------------------------------------------
xindx = linspace(0,1,Nbins)';
options = statset('Display','final','Maxiter',1000);
phat = mle(D595s,'distribution','beta');
PD_D595s = ProbDistUnivParam('beta',phat);

subplot(9,3,[15 18 21]),set(gca,'Fontname','TimesNewRoman','Fontsize',8)
histfit(D595s,15,'beta')
h = get(gca,'Children');
set(h(2),'FaceColor',[.6 .6 0.8])
xlabel('$D_{5-95}$','Interpreter','Latex','Fontsize',11)
ylabel('Frequency','Fontname','TimesNewRoman','Fontsize',11)
text(0.5,160,'Beta(4.15,5.15)','Fontname','TimesNewRoman','Fontsize',11,'HorizontalAlignment','Center') 
axis([min(xindx) max(xindx) 0 170])
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'histfit']),close;
end

clear theta *RSS* extflg zeta_ref EQindx

figure(3),
subplot(221),set(gca,'Fontname','TimesNewRoman')
hold on
plot(D595s,D545s,'+')
coef = D545s/D595s;
plot([0 1],[0 coef],'-r','Linewidth',1.5)
plot([0 1],[0 min(D545s./D595s)],'--r')
plot([0 1],[0 max(D545s./D595s)],'--r')
xlabel('$D_{5-95}$','Interpreter','Latex')
ylabel('$t_{\mbox{mid}}$','Interpreter','Latex')
legend({'Samples','Least Square fit','Min and Max ratio'},'Fontname','TimesNewRoman',...
    'Location','NorthWest')
grid on 
box on
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'D545_D595']),close;
end
save('PDFs.mat','PD*')


%% Space frame modal analysis (linear material properties) 
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\PreliminaryAnalysis\')
% Number of experiments
N = 50; 
% Number of variables (Density, Ex = Ey = Ez, Prxy = Prxz = Prxyz)
M = 2; 

% Latin Hypercube Sampling
Q = lhsdesign(N,M,'criterion','maximin','iterations',10);
Q = 2*(Q'-0.5);
Qtrans = [2e11 + 1e10*Q(1,:);1 + 0.2*Q(2,:)]';
Qtrans = Qtrans(:);
save('BeamProps.txt','Qtrans','-ascii')

% Run Ansys
delete('file.*')
dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "C:\Users\minas\Dropbox\MATLAB\ICOSSAR2013\PreliminaryAnalysis\SpaceFrame2DModal.txt" -o "output.txt"');



%% Space frame -- bilinear isotropic material (increasing force) 
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\PreliminaryAnalysis\')

NoS = 2000;
NoExp = 1;
% Generation of the input force
Finput = linspace(1e2,1e4,NoS)';                                        % Linear
save('InputForce.txt','Finput','-ascii')

% Run Ansys
delete('file.*')
dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\PreliminaryAnalysis\SpaceFrame2DForce.txt" -o "output.txt"');


%% Space frame -- bilinear isotropic material (increasing force) -- Figure 
% =========================================================================
load('SF_ForceResponse.txt');
time = SF_ForceResponse(:,1);
Fi = SF_ForceResponse(:,2)*36;
sumreact = SF_ForceResponse(:,3);
ux = SF_ForceResponse(:,4:8);
Sx = SF_ForceResponse(:,9:13);

figure,
subplot(4,1,1),plot(time,sumreact/1e3)
xlabel('Time (s)')
ylabel('F_{X_{reac}} (kN)')
subplot(4,2,3:2:7),plot(time,ux)
legend({'Floor 1','Floor 2','Floor 3','Floor 4','Floor 5'},'Location','NorthWest')
xlabel('Time (s)')
ylabel('v_X (m)')
axis tight
subplot(4,2,4:2:8),plot(time,Sx/1e6)
legend({'Floor 1','Floor 2','Floor 3','Floor 4','Floor 5'},'Location','NorthWest')
xlabel('Time (s)')
ylabel('\sigma_X (MPa)')
axis tight

figure,
plot(sumreact,ux(:,5))
xlabel('F_{X_{reac}} (kN)')
legend({'Floor 5'},'Location','NorthWest')
ylabel('v_X (m)')

%% Transient nonlinear analysis (various acceleration levels) 
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\')

load('SimulatedEQs.mat','Q*','EQsim','K','N','M','Ts')
% K : Number of experiments
% N: Number of samples
% M: Number of variables (beta, gamma, n)
% Convert synthetic earthquakes in m/s^2
EQsim = EQsim*9.81;

% Latin Hypercube Sampling
Eunc = lhsdesign(K,1,'criterion','maximin','iterations',10);
Eunc = 2*(Eunc'-0.5);
Eunctrans = 1 + 0.1*Eunc(1,:)';
BeamProps = Eunctrans(:);
save('BeamProps.txt','BeamProps','-ascii')

Exc = EQsim(:);
save('AccEQ.txt','Exc','-ascii')

Q = [Q Eunc']; 
UNCpars = [Qtrans Eunctrans];
save('InputVars.mat','Q','Exc','UNCpars')

% Run Ansys
delete('file.*')
dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\SpaceFrameEQ2D.txt" -o "output.txt"');


%% Transient nonlinear analysis (various acceleration levels) -- Figure
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\')
Dt = 0.025;
Fs = 1/Dt;

figure(1),
clr = colormap(lines(3));
subplot(311),hold on
subplot(312),hold on
subplot(313),hold on

K = 200;

for k = 1:K
    SFRes = load(['Responses\SF2D_EQ_Res',num2str(k),'.txt']);
    time = SFRes(:,1);
    acci = SFRes(:,2);
    sumreact = SFRes(:,3);
    ux = SFRes(:,4:8);
    acc =  SFRes(:,9:13);
    [Py(:,k),F] = pwelch(acc(:,5),512,480,512,Fs);
    [Txy(:,k)] = tfestimate(acci,acc(:,5),512,480,512,Fs);
    [Cxy(:,k),F] = mscohere(acci,acc(:,5),512,480,512,Fs);
    if k<4
        figure(1),
        subplot(311),plot(time,acci,'Color',clr(k,:))
        subplot(312),plot(time,sumreact/1000,'Color',clr(k,:))
        subplot(313),plot(time,ux(:,5),'Color',clr(k,:))
    end
end

figure(1),
subplot(311),set(gca,'Xticklabel',[])
ylabel('Input acc. (m/s^2)')
axis([time(1) time(end) -12 12])
subplot(312),ylabel('\Sigma(Reaction forces) (kN)')
set(gca,'Xticklabel',[])
axis([time(1) time(end) -3000 3000])
subplot(313),ylabel('Top floor displacement (m)')
xlabel('Time (s)')
axis([time(1) time(end) -0.8 0.8])


figure(2),
plot(ux(:,5),sumreact/1000,'-o')
ylabel('Total force (kN)')
xlabel('Top floor displacement (m)')
axis tight

figure(3)
surf(1:K,F,20*log10(abs(Py)))
shading interp
zlabel('Magnitude (dB)')
ylabel('Frequency (Hz)')
xlabel('Experiment No.')
title('PSD')
axis tight

figure(4)
surf(1:K,F,20*log10(abs(Txy)))
shading interp
zlabel('Magnitude (dB)')
ylabel('Frequency (Hz)')
xlabel('Experiment No.')
title('FRF')
axis tight

figure(5)
surf(1:K,F,Cxy)
shading interp
zlabel('Coherence function')
ylabel('Frequency (Hz)')
xlabel('Experiment No.')
axis tight


%% Conventional ARX - NARX identification 
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\';
cd(local_dir);

options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';

P = 2;
q = 1;
maxterms = 2;
K = 200;
N = 1000;

disp('Linear    Nonlinear')
disp('===================')
for ii = 1:K
    out = load(['Responses\SF2D_EQ_Res',num2str(ii),'.txt']);
    Y = out(1:N,13);
    X = out(1:N,2);
    na = 10;
    % ARX estimation
    [theta,res,criteria] = narx(Y,X,[na 0 na],[],1,1,[],[],options);
    rss_sss(ii) = criteria.rss_sss;
    % NARX estimation
    [theta,NLres,NLcriteria,regstr] = narx(Y,X,[na 0 na],[],P,q,maxterms,[],options);
    NLrss_sss(ii) = NLcriteria.rss_sss;
    s = sprintf('%2.5f \t %2.5f',rss_sss(ii),NLrss_sss(ii));
    disp(s)
    
end

figure(1)
plot(rss_sss,'-o')
hold on
plot(NLrss_sss,'r-d')
legend({'ARX(10,10)','NARX(10,10)_{(p = 2,maxterms = 2)}'})
xlabel('Simulation experiment')
ylabel('RSS/SSS (%)')


%% FOLS-NARX (nonlinear regressors selection)
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\';
cd(local_dir);

options.maxsize = 1e7;
options.focus = 'prediction';

P = 3;
q = 1;
maxterms = 2;
N = 1000;
K = 200;
disp('Linear    Nonlinear')
disp('===================')
for ii = 1:K
    cd(local_dir);
    out = load(['Responses\SF2D_EQ_Res',num2str(ii),'.txt']);
    Y = out(:,13);
    X = out(:,2);
    na = 10;
    % ARX estimation
    [theta,res,criteria] = narx(Y,X,[na 0 na],[],1,1,[],[],options);
    rss_sss(ii) = criteria.rss_sss;
    % FOLS NARX model structure selection
    options.criterion = 'rss';
    [Nr,RegSTR,regindx{ii},NLres,NLcriteria] = folsNARX(Y,X,[na 0 na],[],P,q,maxterms,40,[],options);
    NLrss_sss(ii) = NLcriteria.rss_sss;
    s = sprintf('%2.5f \t %2.5f',rss_sss(ii),NLrss_sss(ii));
    disp(s)
end

[~,~,~,InitRegstr] = narx(Y,X,[na 0 na],[],P,q,maxterms,[],options);
TotIndx = [];
for ii=1:K
    TotIndx = [TotIndx  regindx{ii}];
end
N = hist(TotIndx,length(InitRegstr));
FinalRegstr  = InitRegstr(N >= 30);

save('FOLSNARXresults.mat','FinalRegstr','InitRegstr','K','N','rss_sss','NLrss_sss','Nr','P','q','maxterms')


%% GA PC functional subspace selection
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\';
cd(local_dir)

% Number of experiments
K = 200; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, synthesized earthquakes parameters)
M = 5; 

[Y,X] = deal(zeros(NoS,K));
for i = 1:K
    out = load([local_dir,'Responses\SF2D_EQ_Res',num2str(i),'.txt'],'-ascii');
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('Responses\InputVars.mat','Q');

options.maxsize = 1e8;
options.basis = 'legen';
options.GA.PopulationSize = 75;
options.GA.Display = 'iter';
options.GA.PlotFcns = {@gaplotbestf,@gaplotbestindiv};

FinalRegstr{1} = '(Y(t-1,:).^1)'; FinalRegstr{2} = '(Y(t-2,:).^1)';
FinalRegstr{3} = '(Y(t-3,:).^1)'; FinalRegstr{4} = '(Y(t-4,:).^1)';
FinalRegstr{5} = '(Y(t-5,:).^1)'; FinalRegstr{6} = '(Y(t-6,:).^1)';
FinalRegstr{7} = '(Y(t-7,:).^1)'; FinalRegstr{8} = '(Y(t-8,:).^1)';
FinalRegstr{9} = '(Y(t-9,:).^1)'; FinalRegstr{10} = '(Y(t-10,:).^1)';
FinalRegstr{11} = '(X(t-0,:).^1)';FinalRegstr{12} = '(X(t-1,:).^1)';
FinalRegstr{13} = '(X(t-2,:).^1)';FinalRegstr{14} = '(X(t-3,:).^1)';
FinalRegstr{15} = '(X(t-4,:).^1)';FinalRegstr{16} = '(X(t-5,:).^1)';
FinalRegstr{17} = '(X(t-6,:).^1)';FinalRegstr{18} = '(X(t-7,:).^1)';
FinalRegstr{19} = '(X(t-8,:).^1)';FinalRegstr{20} = '(X(t-9,:).^1)';
FinalRegstr{21} = '(X(t-10,:).^1)';FinalRegstr{22} = '(Y(t-1,:).^3)';
FinalRegstr{23} = '(Y(t-2,:).^3)';FinalRegstr{24} = '(Y(t-3,:).^3)';
FinalRegstr{25} = '(Y(t-1,:).^1).*(Y(t-2,:).^2)';FinalRegstr{26} = '(Y(t-1,:).^2).*(Y(t-2,:).^1)';
FinalRegstr{27} = '(Y(t-1,:).^1).*(Y(t-3,:).^2)';FinalRegstr{28} = '(Y(t-1,:).^2).*(Y(t-3,:).^1)';
FinalRegstr{29} = '(Y(t-1,:).^1).*(Y(t-4,:).^2)';FinalRegstr{30} = '(Y(t-1,:).^2).*(Y(t-4,:).^1)';
FinalRegstr{31} = '(Y(t-1,:).^1).*(Y(t-5,:).^2)';FinalRegstr{32} = '(Y(t-1,:).^2).*(Y(t-5,:).^1)';
FinalRegstr{33} = '(Y(t-1,:).^1).*(Y(t-6,:).^2)';FinalRegstr{34} = '(Y(t-1,:).^2).*(Y(t-6,:).^1)';
FinalRegstr{35} = '(Y(t-1,:).^1).*(Y(t-7,:).^2)';FinalRegstr{36} = '(Y(t-1,:).^2).*(Y(t-7,:).^1)';
FinalRegstr{37} = '(Y(t-1,:).^1).*(Y(t-8,:).^2)';FinalRegstr{38} = '(Y(t-1,:).^2).*(Y(t-8,:).^1)';
FinalRegstr{39} = '(Y(t-1,:).^1).*(Y(t-9,:).^2)';FinalRegstr{40} = '(Y(t-1,:).^2).*(Y(t-9,:).^1)';
FinalRegstr{41} = '(Y(t-1,:).^1).*(Y(t-10,:).^2)';FinalRegstr{42} = '(Y(t-1,:).^2).*(Y(t-10,:).^1)';

na = 10;
p = 3;
q = 1;
rng(2^10);
[INDX,thetaij,res,criteria,GAoutput] = gaPCNARX(Y,X,Q',[na 0 na],FinalRegstr,p,q,options);
save('GAPCNARXresults.mat','INDX','FinalRegstr','GAoutput');


%% PC-NARX estimation
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\';
cd(local_dir)

% Number of experiments
K = 200; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, synthesized earthquakes parameters)
M = 5; 

[Y,X] = deal(zeros(NoS,K));
for i = 1:K
    out = load([local_dir,'Responses/SF2D_EQ_Res',num2str(i),'.txt'],'-ascii');
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('Responses\InputVars.mat','Q');
load('GAPCNARXresults.mat','INDX','FinalRegstr')

na = 10;
% PE estimation
options.maxsize = 1e8;
options.basis = 'legen';
options.focus = 'prediction';
[THij,res,criteria] = pcnarx(Y,X,Q',[na 0 na],INDX,FinalRegstr,[],options);

% SE estimation
options.focus = 'simulation';
[THij_sim,res_sim,criteria_sim] = pcnarx(Y,X,Q',[na 0 na],INDX,FinalRegstr,THij,options);

save('PCNARXmodel.mat','TH*','criter*','res*',...
    'INDX','FinalRegstr','X','Y','K','M','NoS','Q');


%% PC-NARX estimation results (prediction - simulation errors)
% =========================================================================
clear,clc,close all
local_dir = '/run/media/minas/8A20-0822/Documents/PublishedPapers/ICOSSAR2013/m-files/';
cd(local_dir)
load('PCNARXmodel.mat','FinalRegstr','INDX','THij*','Q')


write_dir = [local_dir,'Figures/'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
K = 200; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, synthesized earthquakes parameters)
M = 5; 
options.basis = 'legen';
na = 10;

[Y,X,RES,RESsim] = deal(zeros(NoS,K));
for i = 1:K
    disp(i)
    out = load([local_dir,'Responses/SF2D_EQ_Res',num2str(i),'.txt'],'-ascii');
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
    [~,theta_pred] = PCparam(Q(i,:)',3,INDX{1},THij,options);
    [~,theta_sim] = PCparam(Q(i,:)',3,INDX{1},THij_sim,options);
    Ysim = narxsim(X(:,i),na,FinalRegstr,theta_sim);
    RESsim(:,i) = Y(:,i) - Ysim;
    [Ypred,E] = narxpred(Y(:,i),X(:,i),na,FinalRegstr,theta_pred);
    RES(:,i) = Y(:,i) - Ypred;
    rss_sss(i) = 100*norm(RES(:,i))^2/norm(Y(:,i))^2; 
    sss_sss(i) = 100*norm(RESsim(:,i))^2/norm(Y(:,i))^2;
end

%%
figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(6,1,3:6),plot(1:K,rss_sss,'-bo','Linewidth',1,'Markersize',4)
hold on,plot(1:K,sss_sss,'r-x','Linewidth',1,'Markersize',4)
xlabel('Experiment #','Fontsize',12,'FontName','TimesNewRoman')
ylabel('Normalized Sum of Squared Errors (%)','Fontsize',12,'FontName','TimesNewRoman')
grid on
axis([0 201 0 120])
set(gca,'Fontsize',10)
legend({'PC-NARX prediction','PC-NARX simulation'},'Location','NorthEast','Fontsize',12)
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'K_preds_sims']),close;
     eps2xxx([write_dir,'K_preds_sims.eps'],{'pdf'})
end

%% El Centro earthquake (ANSYS run)
%==========================================================================
clear,clc,pack
local_dir = '/run/media/minas/8A20-0822/Documents/PublishedPapers/ICOSSAR2013/m-files/';
cd(local_dir)

load('ElCentro.mat');
EQ = ElCentro(:,2); 
N = length(EQ);
time = 1:N;
save('validEQ.txt','EQ','-ascii')

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

load([local_dir,'PDFs.mat'])
% Sampling perios
Ts = 0.025;
% El Centro characteristics
[Ia,t5,t45,t95,Iat,D5_45,D5_95] = EQarias(EQ/9.81,Ts);
options.points = 9;
omegas = spectralINTERP(EQ/9.81,Ts,options);

Qtrans = [omegas(1) omegas(2) sqrt(Ia) D5_95 200e9];
Q(:,1) = cdf(PD_wmid,Qtrans(:,1));
Q(:,2) = cdf(PD_wdot,Qtrans(:,2));
Q(:,3) = cdf(PD_Ia,Qtrans(:,3));
Q(:,4) = cdf(PD_D595s,Qtrans(:,4));
Q(:,5) = 0.5;
Q = -1+2*Q;

% PC-NARX parameters
yr = load('Responses/SF2D_EQ_Res_valid.txt');
load('PCNARXmodel.mat','INDX','THij*','FinalRegstr','criteria*')

% PC-NARX simulation
options.basis = 'legen';
[~,theta0_sim] = PCparam(Q',3,INDX{1},THij_sim,options);
tic;
Ysim = narxsim(EQ,10,FinalRegstr,theta0_sim);
cputime = toc;

% PC-NARX prediction
[~,theta0] = PCparam(Q',3,INDX{1},THij,options);
Ypred = narxpred(yr(:,13),EQ,10,FinalRegstr,theta0);

% RSS/SSS
rss_sss = 100*norm(yr(:,13)-Ypred)^2/norm(yr(:,13))^2;
% SSS/SSS
sss_sss = 100*norm(yr(:,13)-Ysim)^2/norm(yr(:,13))^2;


figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(6,1,3:6),plot(time,yr(:,13),'-','Linewidth',2,'Color',[0.1 0 0])
hold on,plot(time,Ypred,'-','Linewidth',1,'Color',[0 0.6 0.0])
plot(time,Ysim,'-','Linewidth',1,'Color',[0 0 1])
xlabel('Discrete time (samples)','Fontsize',14,'FontName','TimesNewRoman')
ylabel('$\alpha_{12}\ (m/s^2)$','Fontsize',14,'Interpreter','Latex')
grid on
axis([1 600 -2.5 2])
set(gca,'Fontsize',10)
legend({'FE response','PC-NARX prediction','PC-NARX simulation'},'Location','SouthEast','Fontsize',12)
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'preds_sims']),close;
end


%% Input-Output simulation data -- Figure
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\';
write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

Dt = 0.025;
Fs = 1/Dt;

figure(1),
clr = [0.2 0.0 0.0; 0.0 0.6 0.0; 0.0 0.0 1];
style = ['-';'-';'-'];
subplot(15,4,[1:3 5:7 9:11]),hold on
subplot(15,4,[17:19 21:23 25:27]),hold on
subplot(15,4,[33:35 37:39 41:43]),hold on
subplot(15,4,[49:51 53:55 57:59]),hold on

K = 3;

for k = 1:K
    SFRes = load(['Responses\SF2D_EQ_Res',num2str(k),'.txt']);
    time = SFRes(:,1);
    acci = SFRes(:,2);
    sumreact = SFRes(:,3);
    ux = SFRes(:,4:8);
    acc =  SFRes(:,9:13);
    [Py(:,k),F] = pwelch(acc(:,5),512,480,512,Fs);
    [Txy(:,k)] = tfestimate(acci,acc(:,5),512,480,512,Fs);
    [Cxy(:,k),F] = mscohere(acci,acc(:,5),512,480,512,Fs);
    if k<4
        figure(1),
        subplot(15,4,[1:3 5:7 9:11]),plot(time,acci,'Color',clr(k,:),'Linestyle',style(k,:),'Linewidth',1)
        subplot(15,4,[17:19 21:23 25:27]),plot(time,sumreact/1000,'Color',clr(k,:),'Linestyle',style(k,:),'Linewidth',1)
        subplot(15,4,[33:35 37:39 41:43]),plot(time,ux(:,5),'Color',clr(k,:),'Linestyle',style(k,:),'Linewidth',1)
        subplot(15,4,[49:51 53:55 57:59]),plot(ux(:,5),sumreact/1000,'Color',clr(k,:),'Linestyle',style(k,:),'Linewidth',1)
    end
end

figure(1),
subplot(15,4,[1:3 5:7 9:11]),set(gca,'Fontname','TimesNewRoman','Fontsize',9)
xlabel('Time (s)','Fontname','TimesNewRoman','Fontsize',11)
text(-2.5,0,sprintf('%s\n%s','Input acc.','(m/s^2)'),'Fontname','TimesNewRoman','Rotation',90,'Horizontalalignment','Center','Fontsize',11)
axis([time(1) time(end) -5.2 5.2])
grid on;box on ;
h = legend({'1^{st} Simulation','2^{nd} Simulation','3^{rd} Simulation'},'Fontname','TimesNewRoman');
set(h,'Orientation','Horizontal','Position',[0.115 0.95 0.6 0.05],'Fontsize',10);
subplot(15,4,[17:19 21:23 25:27]),set(gca,'Fontname','TimesNewRoman','Fontsize',9)
text(-2.5,0,sprintf('%s\n%s','\Sigma(Reaction forces)','(kN)'),'Fontname','TimesNewRoman','Rotation',90,'Horizontalalignment','Center','Fontsize',11)
xlabel('Time (s)','Fontname','TimesNewRoman','Fontsize',11)
axis([time(1) time(end) -400 400])
grid on;box on ;
subplot(15,4,[33:35 37:39 41:43]),set(gca,'Fontname','TimesNewRoman','Fontsize',9)
text(-2.5,-0.025,sprintf('%s\n%s','Top floor displ.','(m)'),'Fontname','TimesNewRoman','Rotation',90,'Horizontalalignment','Center','Fontsize',11)
xlabel('Time (s)','Fontname','TimesNewRoman','Fontsize',11)
axis([time(1) time(end) -0.15 0.1])
grid on;box on ;
subplot(15,4,[49:51 53:55 57:59]),set(gca,'Fontname','TimesNewRoman','Fontsize',9)
xlabel('Top floor displacement (m)','Fontname','TimesNewRoman','Fontsize',11)
text(-0.176,0,sprintf('%s\n%s','\Sigma(Reaction forces)','(kN)'),'Fontname','TimesNewRoman','Rotation',90,'Horizontalalignment','Center','Fontsize',11)
axis([-0.15 0.1 -400 400])
grid on;box on ;
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'SF2D_simulations']),close;
end

%% Input random parameters -- Figure
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\';
cd(local_dir)

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';
% Number of experiments
N = 200; 

load('Responses\InputVars.mat','UNCpars');
% Material 1: Vertical beams 
EX  = 200E9*UNCpars(:,5);        % N/m^2	

figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(10,8,[1:4 9:12]),plot(1:N,EX/1e9)
ylabel('$E_{x}$ (GPa)','Fontsize',11,'Interpreter','Latex')
set(gca,'Xticklabel',[],'Fontsize',7,'xlim',[1,N],'ylim',[180 220])
subplot(10,8,[5 6 13 14]),
[Nhist,Xhist] = hist(EX,15);
barh(Xhist,Nhist)
set(gca,'yticklabel',[],'Xticklabel',[],'Fontsize',7,'xlim',[1,70],'ylim',[180e9 220e9])
subplot(10,8,[17:20 25:28]),plot(1:N,UNCpars(:,1)/2/pi)
ylabel('$\omega_{\mbox{mid}}$ (Hz)','Fontsize',11,'Interpreter','Latex')
set(gca,'Xticklabel',[],'Fontsize',7,'xlim',[1,N],'ylim',[0 22])
subplot(10,8,[21 22 29 30]),
[Nhist,Xhist] = hist(UNCpars(:,1)/2/pi,15);
barh(Xhist,Nhist)
set(gca,'yticklabel',[],'Xticklabel',[],'Fontsize',7,'xlim',[1,70],'ylim',[0 15])

subplot(10,8,[33:36 41:44]),plot(1:N,UNCpars(:,2)/2/pi)
set(gca,'Xticklabel',[],'Fontsize',7,'xlim',[1,N],'ylim',[-0.25 0.25])
ylabel('$\omega''$ (Hz/s)','Fontsize',11,'Interpreter','Latex')
subplot(10,8,[37 38 45 46]),
[Nhist,Xhist] = hist(UNCpars(:,2)/2/pi,15);
barh(Xhist,Nhist)
set(gca,'yticklabel',[],'Xticklabel',[],'Fontsize',7,'xlim',[1,70],'ylim',[-0.25 0.25])

subplot(10,8,[49:52 57:60]),plot(1:N,UNCpars(:,3))
set(gca,'Xticklabel',[],'Fontsize',7,'xlim',[1,N],'ylim',[0 0.035])
ylabel('$\sqrt{I_a}$','Fontsize',11,'Interpreter','Latex')
subplot(10,8,[53 54 61 62]),
[Nhist,Xhist] = hist(UNCpars(:,3),15);
barh(Xhist,Nhist)
set(gca,'Xticklabel',[],'yticklabel',[],'Fontsize',7,'xlim',[1,70],'ylim',[0 0.035])

subplot(10,8,[65:68 73:76]),plot(1:N,UNCpars(:,4))
set(gca,'Fontsize',7,'xlim',[1,N],'ylim',[0 0.9])
xlabel('Experiment number','Fontsize',11,'Fontname','TimesNewRoman')
ylabel('$D_{5-95}$','Fontsize',11,'Interpreter','Latex')
subplot(10,8,[69 70 77 78]),
[Nhist,Xhist] = hist(UNCpars(:,4),15);
barh(Xhist,Nhist)
set(gca,'yticklabel',[],'Fontsize',7,'xlim',[1,70],'ylim',[0 0.9])
xlabel('Frequency','Fontsize',11,'Fontname','TimesNewRoman')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'input_pars']),close;
end


%% Uncrosscorrelatedeness (estimated covariance matrix) --- OK
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\';
cd(local_dir)

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r150';

load('PCNARXmodel.mat')
RES = reshape(res,200,NoS-10)';

figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
surf(cov(RES))
clr = colormap(gray);
colormap(flipud(clr));
shading interp
set(gca,'Fontsize',9,'projection','Perspective')
view([60 75])
xlim([1 200])
ylim([1 200])
zlim([-0.005 0.025])
text(100,-40,0,'row index $k_1$','Interpreter','Latex','Fontsize',14,'Rotation',-50,'Horizontalalignment','center')
text(240,80,0,'column index $k_2$','Interpreter','Latex','Fontsize',14,'Rotation',22,'Horizontalalignment','center')
zlabel('cov($\hat{e}_{k_1}\; , \hat{e}_{k_2} \;$)','Interpreter','Latex','Fontsize',14)
grid on
if print_to_eps=='y';
     print(pictype,resolution,'-zbuffer',[write_dir,'covariance']),close;
end


%% Parameter surfaces
print_to_eps = 'y';
clear Qsurf
NoI1 = 200;
NoI2 = 200;
aux = linspace(-1,1,NoI1);
Qsurf(1,:) = kron(ones(1,NoI2),aux);
Qsurf(3,:) = zeros(1,NoI1*NoI2);
Qsurf(2,:) = zeros(1,NoI1*NoI2);
Qsurf(4,:) = zeros(1,NoI1*NoI2);
Qsurf(5,:) = kron(aux,ones(1,NoI2));

cd('C:\Users\sminas\Dropbox\MATLAB\PCE\m-files\');
[~,an_interp] = PCparam(Qsurf,3,INDX{1},aij,options);
    
figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
surf(linspace(0.9,1.1,NoI2),linspace(0.9,1.1,NoI1),reshape(an_interp(:,1),NoI1,NoI2))
shading interp
set(gca,'Fontsize',6)
xlim([0.9 1.1])
ylim([0.9 1.1])
xlabel('$h$ (\%)','Interpreter','Latex','Fontsize',7)
ylabel('$E_x$ (\%)','Interpreter','Latex','Fontsize',7)
zlabel(['$a_',num2str(1),'(\xi)$'],'Interpreter','Latex','Fontsize',7)
grid on
box on
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'ARXsurf']),close;
     result = eps2xxx([write_dir,'ARXsurf.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% Parameter curves
print_to_eps = 'y';
clear Qsurf NoI*
NoI1 = 1000;
Qsurf(5,:) = zeros(1,NoI1);
Qsurf(2,:) = zeros(1,NoI1);
Qsurf(3,:) = zeros(1,NoI1);
Qsurf(4,:) = zeros(1,NoI1);
Qsurf(1,:) = linspace(-1,1,NoI1);

[~,an_interp] = PCparam(Qsurf,3,INDX{1},aij,options);

figure(2)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
plot(linspace(190,210,NoI1),an_interp(:,4),'-')
xlim([190 210])
set(gca,'fontsize',6)
xlabel('$h$ (\%)','interpreter','latex','fontsize',12)
ylabel(['$\theta_',num2str(1),'(\xi)$'],'interpreter','latex','fontsize',12)
grid on

if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'ARXcurve']),close;
     result = eps2xxx([write_dir,'ARXcurve.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% Steel Properties
% Properties                    Carbon     Alloy        Stainless     Tool
% ==========================================================================
% Density (1000 kg/m3)        7.85          7.85        7.75-8.1    7.72-8.0
% Elastic Modulus (GPa)       190-210       190-210     190-210     190-210
% Poisson's Ratio             0.27-0.3      0.27-0.3    0.27-0.3    0.27-0.3
% Tensile Strength (MPa)      276-1882      758-1882    515-827     640-2000
% Yield Strength (MPa)        186-758       366-1793    207-552     380-440
% ==========================================================================

%% ANSYS SpaceFrameRes*.tst file output
% 1     Time
% 2     Input acceleration 
% 3     Sum of reaction Forces 
% 4     Displacement of 1st floor (node 4)
% 5     Displacement of 2nd floor (node 6)
% 6     Displacement of 3rd floor (node 8)
% 7     Displacement of 4th floor (node 10)
% 8     Displacement of 5th floor (node 12)
% 9     Acceleration response of 1st floor (node 4)
% 10    Acceleration response of 2nd floor (node 6)
% 11    Acceleration response of 3rd floor (node 8)
% 12    Acceleration response of 4th floor (node 10)
% 13    Acceleration response of 5th floor (node 12)
% 13    Acceleration response of 5th floor (node 54)