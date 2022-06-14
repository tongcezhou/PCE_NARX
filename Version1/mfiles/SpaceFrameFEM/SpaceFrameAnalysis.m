%% Space frame modal analysis (linear material properties)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\')
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
dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\SpaceFrameModal.f" -o "output.txt"');


%% Space frame modal analysis (linear material properties) -- Figure
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\')
load('Wn.txt')
figure
subplot(211),plot(2*pi*Wn(1:2:end,:),'-o')
xlabel('Experiment number')
ylabel('Natural Frequencies (Hz)')
grid on
subplot(212),plot(100*Wn(2:2:end,:),'-o')
xlabel('Experiment number')
ylabel('Damping ratios (%)')
grid on


%% Space frame -- bilinear isotropic material (increasing force)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\')

NoS = 2000;
NoExp = 1;
% Generation of the input force
Finput = linspace(1000,2.5e4,NoS)';                                        % Linear
save('InputForce.txt','Finput','-ascii')

% Run Ansys
delete('file.*')
dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\SpaceFrameForce.f" -o "output.txt"');


%% Space frame -- bilinear isotropic material (increasing force) -- Figure
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\')
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
load('SF_ForceResponse.txt');
time = SF_ForceResponse(:,1);
Fi = SF_ForceResponse(:,2);         % *135;
sumreact = SF_ForceResponse(:,3);
ux = SF_ForceResponse(:,4:8);
Sx = SF_ForceResponse(:,9:13);

figure,
subplot(8,1,1:2),plot(time,sumreact/1e3,'Linewidth',2)
set(gca,'xticklabel',[],'Fontname','TimesNewRoman')
ylabel('$\sum F_{X_\mathrm{reac}}$ (kN)','Interpreter','Latex')
grid on
subplot(8,1,3:5),
plot(time,ux(:,1),'-',time,ux(:,2),':',time,ux(:,3),'-.',time,ux(:,4),'--',time,ux(:,5),'-','Linewidth',1)
legend({'1^{st} Storey','2^{nd} Storey','3^{rd} Storey','4^{th} Storey','5^{th} Storey'},'Location','NorthWest')
set(gca,'xticklabel',[],'Fontname','TimesNewRoman')
ylabel('$v_X$ (m)','Interpreter','Latex')
axis tight
grid on
subplot(8,1,6:8),
plot(time,Sx(:,1)/1e6,'-',time,Sx(:,2)/1e6,':',time,Sx(:,3)/1e6,'-.',time,Sx(:,4)/1e6,'--',time,Sx(:,5)/1e6,'-','Linewidth',1)
legend({'1^{st} Storey','2^{nd} Storey','3^{rd} Storey','4^{th} Storey','5^{th} Storey'},'Location','NorthWest')
xlabel('Time (s)','Fontname','TimesNewRoman')
set(gca,'Fontname','TimesNewRoman')
ylabel('$\sigma_X$ (MPa)','Interpreter','Latex')
axis tight
grid on
print('-depsc2','-r300',[write_dir,'force.eps'])
eps2xxx([write_dir,'force.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe')

%% Transient analysis through ANSYS (increasing acceleration) 
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\')

NoS = 2000;
NoExp = 1;

% Generation of the random signal excitation
InputAcc = (linspace(0,10,NoS).*sin(linspace(-10*pi,10*pi,NoS)))';
% InputAcc = 10.*sin(linspace(-10*pi,10*pi,NoS-200))';
save('InputAcc.txt','InputAcc','-ascii')

% Run Ansys
delete('file.*')
dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\SpaceFrame.f" -o "output.txt"');


%% Transient analysis through ANSYS (increasing acceleration) -- Figure
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\')

load('SF_SinusForceResponse.txt')
time = SF_SinusForceResponse(:,1);
acci = SF_SinusForceResponse(:,2);
sumreact = SF_SinusForceResponse(:,3);
ux = SF_SinusForceResponse(:,4:8);
acc =  SF_SinusForceResponse(:,9:13);

figure,
subplot(221),plot(time,acci)
xlabel('Time (s)')
ylabel('Input acc. (x-axis; m/s^2)')
subplot(222),plot(time,sumreact/1e3)
xlabel('Time (s)')
ylabel('Sum of reaction forces (x-axis; kN)')
subplot(223),plot(time,ux)
legend({'Floor 1','Floor 2','Floor 3','Floor 4','Floor 5'})
xlabel('Time (s)')
ylabel('Displacement (ES node; m)')
subplot(224),plot(sumreact(1:1000),ux(1:1000,5))
xlabel('Sum of reaction forces (x-axis; N)')
ylabel('Top floor displacement (ES node; m)')
axis tight
print('-dpng','-r300','sinus_increasing_force.png')



%% Transient nonlinear analysis (various acceleration levels)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\')

% Number of experiments
N = 50; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 

% Latin Hypercube Sampling
Q = lhsdesign(N,M,'criterion','maximin','iterations',10);
Q = 2*(Q'-0.5);
Qtrans = [1 + 0.05*Q(1,:); 1 + 0.2*Q(2,:); 1 + 0.2*Q(3,:)]';
BeamProps = Qtrans(:);
save('BeamProps.txt','BeamProps','-ascii')

PDFsigma = ProbDistUnivParam('lognormal',[0.7 0.5]);
xindx = 0:0.001:20;     Fx = pdf(PDFsigma,xindx);
sigma = icdf('lognormal',0.5*Q(4,:)+0.5,0.7,0.5);

% Generation of the random signal excitation
for k=1:N
    rng(k)
    Exc((k-1)*NoS+1:k*NoS,1) = sigma(k)*randn(NoS,1);
end
save('AccRandom.txt','Exc','-ascii')

UNCpars = [Qtrans sigma'];
Exc = reshape(Exc,NoS,N);
save('InputVars.mat','Q','Exc','UNCpars')

% Run Ansys
delete('file.*')
dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "D:\SpaceFrame\SpaceFrame.txt" -o "output.txt"');



%% Transient nonlinear analysis (various acceleration levels) -- Figures
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\Responses50')
Dt = 0.025;
Fs = 1/Dt;

figure(1),
clr = colormap(lines(3));
subplot(311),hold on
subplot(312),hold on
subplot(313),hold on

for k = 1:6
    SFRes = load(['SFres',num2str(k),'.txt']);
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
surf(1:200,F,20*log10(abs(Py)))
shading interp
zlabel('Magnitude (dB)')
ylabel('Frequency (Hz)')
xlabel('Experiment No.')
title('PSD')
axis tight

figure(4)
surf(1:200,F,20*log10(abs(Txy)))
shading interp
zlabel('Magnitude (dB)')
ylabel('Frequency (Hz)')
xlabel('Experiment No.')
title('FRF')
axis tight

figure(5)
surf(1:200,F,Cxy)
shading interp
zlabel('Coherence function')
ylabel('Frequency (Hz)')
xlabel('Experiment No.')
axis tight


%% Material properties plot
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
N = 50; 
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 

load('Responses50/InputVars.mat','UNCpars');
% Material 1: Vertical beams 
EX  = 200E9*UNCpars(:,1);        % N/m^2	
% Yielding stress
sy  = 250e6*UNCpars(:,2);        % N/m^2
p = 1;
for h =  [0.30:-0.05:0.20 0.20 0.20]
    height(:,p) = h*UNCpars(:,3);
    p = p+1;
end
height(:,[6 7]) = 0.10*ones(N,2);
sigma = UNCpars(:,4);


figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(9,1,1:3),plot(1:N,height(:,1),'-',1:N,height(:,2),'--',1:N,height(:,3),'-.','Linewidth',1.2)
set(gca,'Fontsize',7,'xlim',[1,N],'ylim',[0.15 0.36])
ylabel('Cross-section height (m)','Fontsize',9)
h = legend({'1^{st}';'2^{nd}';'3^{rd} - 4^{th} - 5^{th}'},'Location','NorthOutside','Orientation','Horizontal','Fontsize',9);
set(h,'position',[0.25 0.95 0.5 0.05])
set(gca,'xticklabel',[],'Fontname','TimesNewRoman')

subplot(9,1,4:5),plot(1:N,EX/1e9,'Linewidth',1.2)
ylabel('$E_{x}$ (GPa)','Fontsize',9,'Fontsize',9,'Interpreter','Latex')
set(gca,'Fontsize',7,'xlim',[1,N],'ylim',[190 210])
set(gca,'xticklabel',[],'Fontname','TimesNewRoman')

subplot(9,1,6:7),plot(1:N,sy/1e6,'Linewidth',1.2)
set(gca,'Fontsize',7,'xlim',[1,N],'ylim',[200 300])
ylabel('$\sigma_Y$ (MPa)','Fontsize',9,'Interpreter','Latex')
set(gca,'xticklabel',[],'Fontname','TimesNewRoman')

subplot(9,1,8:9),plot(1:N,sigma,'Linewidth',1.2),hold on
[maxs,Jmax] = max(sigma);
[mins,Jmin] = min(sigma);
plot(Jmin,mins,'ro',Jmax,maxs,'ro','Linewidth',1.2)
set(gca,'Fontsize',7,'xlim',[1,N],'ylim',[0 8])
ylabel('$\sigma_{\alpha}$','Fontsize',9,'Fontsize',9,'Interpreter','Latex')
xlabel('Experiment number','Fontsize',9)

if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'uncertain_props']),close;
     result = eps2xxx([write_dir,'uncertain_props.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end




%% Time histories plot (Node 54 acceleration response) - [min(sigma_a) and max(sigma_a) cases]
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
N = 3; 
% Number of samples
NoS = 2500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

load('Responses50\InputVars.mat','UNCpars');
sigma = UNCpars(:,4);
[maxs,Jmax] = max(sigma);
[mins,Jmin] = min(sigma);

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
clr = colormap(lines(3));
subplot(5,2,1:2:3),hold on
set(gca,'xticklabel',[],'Fontsize',10,'Fontname','TimesNewRoman')
grid on
subplot(5,2,5:2:7),hold on
xlabel('Time (s)','Fontname','TimesNewRoman')
set(gca,'Fontsize',10,'Fontname','TimesNewRoman')
grid on

linestyle = {'- ','--'};
p = 0;
for i = [Jmin Jmax]
    p = p+1;
    % Load data
    out = load(['Responses50\SFres',num2str(i),'.txt']);
    time = out(:,1);
    acci = out(:,2);
    acc =  out(:,13);
    
    figure(1)
    subplot(5,2,[1 3] + (p-1)*4),plot(time,acc,'Color',clr(p,:))
    box on
    P(:,p) = pwelch(acc,512,480,512,1/Ts);
    Txy(:,p) = tfestimate(acci,acc,512,480,512,1/Ts);
    [Cxy(:,p),F] = mscohere(acci,acc,512,480,512,1/Ts);
end

% Plot signals and PSDs
figure(1)
subplot(5,2,1:2:3),box on,% axis([time(1) time(end) -12 12])
ylabel('$y[t]\ \mathrm{ (m/s^2)}$','Interpreter','Latex','Fontsize',12)
legend({['Exp. # ',num2str(Jmin)]})
subplot(5,2,5:2:7),box on,% axis([time(1) time(end) -0.025 0.025])
ylabel('$y[t]\ \mathrm{ (m/s^2)}$','Interpreter','Latex','Fontsize',12)
legend({['Exp. # ',num2str(Jmax)]})
subplot(5,2,2:2:4),plot(F,20*log10(abs(Txy)),'linewidth',0.5)
axis([0 1/(2*Ts) -60 20])
grid on
ylabel('FRF magnitude (dB)','Fontname','TimesNewRoman')
set(gca,'Fontsize',10,'Fontname','TimesNewRoman')
set(gca,'xticklabel',[],'Fontsize',10,'Fontname','TimesNewRoman')
legend({['Exp. # ',num2str(Jmin)],['Exp. # ',num2str(Jmax)]})

subplot(5,2,6:2:8),plot(F,Cxy,'linewidth',0.5)
axis([0 1/(2*Ts) 0 1])
grid on
ylabel('Coherence','Fontname','TimesNewRoman')
xlabel('Frequency (Hz)','Fontname','TimesNewRoman')
set(gca,'Fontsize',10,'Fontname','TimesNewRoman')
legend({['Exp. # ',num2str(Jmin)],['Exp. # ',num2str(Jmax)]},'Location','SouthWest')
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'signals']),close;
    result = eps2xxx([write_dir,'signals.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end



%% PC-NARX model structure selection (FOLS)
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

% Number of experiments
N = 50; 
% Number of samples
NoS = 500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

% Input-Output data
[Y,X] = deal(zeros(NoS,N));
for i = 1:N
    out = load(['Responses50/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('Responses50/InputVars.mat','Q');

% Nonlinear polynomial regressors
P = 3;      q = 1;      maxterms = 1;
% Model orders
na = 10;     nb = 10;     nd = 0;
% Nonlinear regressors (complete vector)
regressors = polyregNARX([na nd nb],P,q,maxterms);
% FOLS structure selection
options.criterion = 'rss';
PCorder = 3;
indx{1} = combinations(M,PCorder,1);
[Nr,RegSTR,regindx,basisindx,NLres,NLcriteria] = folsPCNARX(Y,X,Q,[na nd nb],indx,[],[],[],regressors,[],[],options);

% Run candidate models for 1-50 most important regressors
options.maxsize = 1e9;
options.wls_iter = 5;
options.method = 'wls';
options.focus = 'prediction';
options.PEinit = 'n';
for j = 1:length(regindx)
    disp(j)
    PCBASES = {indx{1}(basisindx(1:j),:)};
    NLREGRESSORS = regressors(regindx(1:j));
    options.focus = 'prediction';
    [thetaPE,resPE{j},PEcriteria{j}] = pcnarxuncommon(Y,X,Q,[na nd nb],NLREGRESSORS,PCBASES,[],options);
    rssPE(j) = PEcriteria{j}.rss_sss;
    bicPE(j) = PEcriteria{j}.bic;
%     options.focus = 'simulation';
%     if j == 1
%         TH0 = thetaPE.wls;
%     else
%         TH0 = thetaPE.wls;
%         TH0 = [thetaSIM{j-1}.sim;0];
%     end
%     [thetaSIM{j},resSIM{j},SIMcriteria{j}] = pcnarxuncommon(Y,X,Q,[na nd nb],NLREGRESSORS,PCBASES,TH0,options);
    [resSIM{j},SIMcriteria{j}] = pcnarxuncommonSIM(thetaPE.wls,Y,X,Q,[na nd nb],NLREGRESSORS,PCBASES,options);
    rssSIM(j) = SIMcriteria{j}.rss_sss;
    bicSIM(j) = SIMcriteria{j}.bic;
end

figure
plot(1:100,rssPE(1:100),'-bo',1:100,rssSIM(1:100),'-rd')
set(gca,'yscale','log')
%%
save('SpaceFrame_RandomExc_LHS_FOLS.mat','P','q','maxterms','na','nb','nd','N','regressors','PCorder',...
      'indx','regindx','basisindx','options','res*','rss*','bic*','theta*','*criteria')


%% PC-NARX model structure selection (FOLS)
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

% Number of experiments
N = 50; 
% Number of samples
NoS = 500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

% Input-Output data
[Y,X] = deal(zeros(NoS,N));
for i = 1:N
    out = load(['Responses50/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('Responses50/InputVars.mat','Q');

% Nonlinear polynomial regressors
P = 3;      q = 1;      maxterms = 1;
% Model orders
na = 40;     nb = 40;     nd = 0;
% Nonlinear regressors (complete vector)
regressors = polyregNARX([na nd nb],P,q,maxterms);

% FOLS structure selection
options.criterion = 'rss';
PCorder = 3;
indx{1} = combinations(M,PCorder,1);
[Nr,RegSTR,regindx,basisindx,NLres,NLcriteria] = folsPCNARX(Y,X,Q(:,1:N),[na nd nb],indx,P,q,maxterms,[],200,[],options);

%% Run candidate models for 1-50 most important regressors
options.maxsize = 1e9;
options.wls_iter = 5;
options.method = 'wls';
options.focus = 'prediction';
for j = 1:length(regindx)
    disp(j)
    PCBASES = {indx{1}(basisindx(1:j),:)};
    NLREGRESSORS = regressors(regindx(1:j));
    [thetaNL,NLresPE{j},PEcriteria] = pcnarxuncommon(Y,X,Q,[na nd nb],NLREGRESSORS,PCBASES,[],options);
    rssPE(j) = PEcriteria.rss_sss;
    bicPE(j) = PEcriteria.bic;
    [NLresSIM{j},SIMcriteria] = pcnarxuncommonSIM(thetaNL.wls,Y,X,Q,[na nd nb],NLREGRESSORS,PCBASES,options);
    rssSIM(j) = SIMcriteria.rss_sss;
    bicSIM(j) = SIMcriteria.bic;
end

figure
plot(1:length(regindx),rssPE,'-o',1:length(regindx),rssSIM,'-rd')
set(gca,'yscale','log')

%% PC-NARX estimation (full space)
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

% Number of experiments
N = 50; 
% Number of samples
NoS = 500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

% Input-Output data
[Y,X] = deal(zeros(NoS,N));
for i = 1:N
    out = load(['Responses50/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('Responses50/InputVars.mat','Q');

% Nonlinear polynomial regressors
P = 3;      q = 1;      maxterms = 1;
% Model orders
na = 60;     nb = 60;     nd = 0;
options.method = 'ols';
options.focus = 'prediction';
options.maxsize = 1e9;  
options.wls_iter = 5;

for ord = 40
    for pow = 3
        for pc = 3
            disp([ord,pow,pc])
            % Model orders
            na = ord;     nb = ord;     nd = 0;
            P = pow;      q = 1;      maxterms = 1;
            regressors = polyregNARX([na nd nb],P,q,maxterms);
            options.criterion = 'rss';
            indx{1} = combinations(M,pc,1);
            indx{2} = indx{1};
            [TH,res,criteriaPE,output] = pcnarx(Y,X,Q,[na 0 na],indx,regressors,[],options);
            rssPE(ord,pow,pc) = criteriaPE.rss_sss;
            [resSIM,criteriaSIM] = pcnarxSIM(TH.ols,Y,X,Q,[na 0 na],indx,regressors,options);
            rssSIM(ord,pow,pc) = criteriaSIM.rss_sss;
        end
    end
end


%%
figure
subplot(211),plot(2:20,rssPE(2:20,:,1),'-bo',2:20,rssPE(2:20,:,2),'-rd',2:20,rssPE(2:20,:,3),'-k+')
subplot(212),plot(2:20,rssSIM(2:20,:,1),'-bo',2:20,rssSIM(2:20,:,2),'-rd',2:20,rssSIM(2:20,:,3),'-k+')
ylim([0 100])
% set(gca,'yscale','log')


%% PC-NARX estimation (na = 20)
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

% Number of experiments
N = 200; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

% Input-Output data
[Y,X] = deal(zeros(NoS,N));
for i = 1:N
    out = load(['Responses/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('InputVars.mat','Q');

% Nonlinear polynomial regressors
P = 1;      q = 1;      maxterms = 1;
% Model orders
na = 100;     nb = 100;     nd = 0;
% PC degree
pc = 3;
options.method = 'ols';
options.focus = 'prediction';
options.maxsize = 1e7;  
options.wls_iter = 5;
regressors = polyregNARX([na nd nb],P,q,maxterms);
indx{1} = combinations(M,pc,1);
indx{2} = indx{1};
[TH,res,criteriaPE,output] = pcnarx(Y,X,Q,[na 0 na],indx,regressors,[],options);

options.focus = 'simulation';
[TH,res,criteriaPE,output] = pcnarx(Y,X,Q,[na 0 na],indx,regressors,TH.ols,options);

[resSIM,criteriaSIM] = pcnarxSIM(TH.ols,Y,X,Q,[na 0 na],indx,regressors,options);



%% PC-ARX Identification
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

write_dir = [local_dir,'Figures\'];
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
N = 200; 
% Number of samples
NoS = 500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

[Y,X] = deal(zeros(NoS,N));
for i = 1:N
    out = load(['Responses/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('InputVars.mat','Q');

n_max = 100;
q = 1;
options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';

% ARX modelling
for na = 1:n_max
    for p = 1:4
        disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
        indx{1} = combinations(M,p,q);
        indx{2} = indx{1};
        B = length(indx{1});
        [aij,bij,res,criteria] = pcarx(Y,X,Q,[na 0 na],indx,[],options);
        rss_sss(na,p) = criteria.rss_sss;
        connum(na,p) = criteria.connum;
        spp(na,p) = criteria.spp;
        sigma = var(res);
        bic(na,p) = log(sigma) + ((na+na+1)*B)*log(N*NoS)/(N*NoS);
    end
    save('PCARXresults','rss_sss','bic','spp','connum');
end

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
clr = colormap(lines(4));
subplot(7,1,1:3),plot(2:n_max,rss_sss(2:n_max,1),'-o','Color',clr(1,:)),hold on
plot(2:n_max,rss_sss(2:n_max,2),'-d','Color',clr(2,:))
plot(2:n_max,rss_sss(2:n_max,3),'-s','Color',clr(3,:))
plot(2:n_max,rss_sss(2:n_max,4),'-*','Color',clr(4,:))
text(11,6.5,'PC-ARX($n,n$) models','Interpreter','Latex','Fontsize',12,'HorizontalAlignment','center')
ylabel('RSS/SSS (\%)','Interpreter','Latex','Fontsize',12)
set(gca,'xticklabel',[],'Fontsize',10,'Fontname','Timesnewroman','xlim',[1.5 50.5])
grid on
legend({'$P=1$','$P=2$','$P=3$','$P=4$'},'Interpreter','Latex','Fontsize',10)
subplot(7,1,4:6),plot(2:n_max,connum(2:n_max,1),'-o','Color',clr(1,:)),hold on
plot(2:n_max,connum(2:n_max,2),'-d','Color',clr(2,:))
plot(2:n_max,connum(2:n_max,3),'-s','Color',clr(3,:))
plot(2:n_max,connum(2:n_max,4),'-*','Color',clr(4,:))
ylabel('$\mathbf \Phi(\mathbf \Xi)$ condition number','Interpreter','Latex','Fontsize',12)
set(gca,'yscale','log','Fontsize',10,'Fontname','Timesnewroman','xlim',[1.5 50.5])
grid on
legend({'$P=1$','$P=2$','$P=3$','$P=4$'},'Interpreter','Latex','Fontsize',10)
xlabel('AR/X models order ($n$)','Interpreter','Latex','Fontsize',12)
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'MSScriteria']),close;
     result = eps2xxx([write_dir,'MSScriteria.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end





%% FOLS-NARX (prediction simulation criteria)
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

write_dir = [local_dir,'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\'];
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
K = 200; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

[Y,X] = deal(zeros(NoS,K));
for i = 1:K
    out = load(['Responses/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end

for k = 1:K
    [Cxy,F] = mscohere(X(:,k),Y(:,k),128,120,256,1/Ts);
    NLmeas(k) = 1 - mean(Cxy(F<40));
end
plot(1:K,NLmeas,'-o');
ii = input('Given run with stronger nonlinearities: ');
X = X(1:NoS,ii);
Y = Y(1:NoS,ii); 

options.maxsize = 1e7;
options.focus = 'prediction';
%%
for ii = 10:100
    for jj = 1:3
        % Nonlinear polynomial regressors
        P = jj;      q = 1;      maxterms = 1;
        % Model orders
        na = ii;     nb = ii;     nd = 0;
        regressors = polyregNARX([na nd nb],P,q,maxterms);
        [thetaNL,NLresPE,PEcriteria] = narx(Y,X,[na nd nb],regressors,[],[],[],[],options);
        rss_sss(ii,jj) = PEcriteria.rss_sss;
        % acf(NLresPE,100,0.8)
    end
end

%%
for jj = 3
    % Nonlinear polynomial regressors
    P = 3;      q = 1;      maxterms = jj;
    % Model orders
    na = 41;     nb = 41;     nd = 0;
    regressors = polyregNARX([na nd nb],P,q,maxterms);
    [thetaNL,NLresPE,PEcriteria] = narx(Y,X,[na nd nb],regressors,[],[],[],[],options);
    acf(NLresPE,100,0.8)
end


%% FOLS-NARX (prediction simulation criteria)
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

write_dir = [local_dir,'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\'];
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
K = 200; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

[Y,X] = deal(zeros(NoS,K));
for i = 1:K
    out = load(['Responses/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end

for k = 1:K
    [Cxy,F] = mscohere(X(:,k),Y(:,k),128,120,256,1/Ts);
    NLmeas(k) = 1 - mean(Cxy(F<40));
end
plot(1:K,NLmeas,'-o');
ii = input('Given run with stronger nonlinearities: ');
X = X(1:NoS,ii);
Y = Y(1:NoS,ii); 

options.maxsize = 1e7;
options.focus = 'prediction';
% Nonlinear polynomial regressors
P = 1;      q = 1;      maxterms = 1;
% Model orders
na = 60;     nb = 60;     nd = 0;
regressors = polyregNARX([na nd nb],P,q,maxterms);

options.criterion = 'rss';
[Nr,RegSTR,regindx,NLres,NLcriteria] = folsNARX(Y,X,[na nd nb],[],P,q,maxterms,100,[],options);

for j = 1:length(regindx)
    disp(j)
    options.focus = 'prediction';
    [thetaNL,NLresPE{j},PEcriteria] = narx(Y,X,[na nd nb],regressors(regindx(1:j)),[],[],[],[],options);
    rssPE(j) = PEcriteria.rss_sss;
    bicPE(j) = PEcriteria.bic;
    options.focus = 'simulation';
    options.PEinit = 'y';
    [thetaNL,NLresSIM{j},SIMcriteria] = narx(Y,X,[na nd nb],regressors(regindx(1:j)),[],[],[],[],options);
    rssSIM(j) = SIMcriteria.rss_sss;
    bicSIM(j) = SIMcriteria.bic;
end

save('SpaceFrame_RandomExc_LHS_FOLS_results_linear.mat','regressors',...
     'NLres*','regindx','rss*','bic*','na','nb','nd','P','q','maxterms','NoS','Ts','K','ii');

 
%% FOLS-NARX (prediction simulation criteria) -- Figures
% =========================================================================
load('SpaceFrame_RandomExc_LHS_FOLS_results.mat');
n = 100;
figure(1)
subplot(211),plot(1:n,rssPE(1:n),'-o',1:n,rssSIM(1:n),'-rd')
set(gca,'yscale','log')
xlim([1 n])
subplot(212),plot(1:n,bicPE(1:n),'-o',1:n,bicSIM(1:n),'-rd')
xlim([1 n])

figure(2)
subplot(211),plot(2:n,abs(diff(rssPE))./(1+abs(rssPE(1:n-1))),'-o',2:n,abs(diff(rssSIM))./(1+abs(rssSIM(1:n-1))),'-rd');
set(gca,'yscale','log')
xlim([1 n])
subplot(212),plot(2:n,abs(diff(bicPE))./abs(bicPE(1:n-1)),'-o',2:n,abs(diff(bicSIM))./abs(bicSIM(1:n-1)),'-rd');
 

%% GA PC-ARX Identification
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

write_dir = [local_dir,'Figures\'];
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
N = 200; 
% Number of samples
NoS = 500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

[Y,X] = deal(zeros(NoS,N));
for i = 1:N
    out = load(['SpaceFrameRes',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('InputVars.mat','Q');

options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';
options.criterion = 'RSS';
options.common = 'y';
options.GA.PopulationSize = 100 ;
options.GA.Display = 'iter';
options.GA.PlotFcns = {@gaplotbestf,@gaplotbestindiv};

n_max = 100;
p = 4;
q = 1;

cd('C:\Users\sminas\Dropbox\MATLAB\PCE\m-files\');
% ARX modelling
for na = 58:n_max
    disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
    indx{1} = combinations(M,p,q);
    indx{2} = indx{1};
    B = length(indx{1});
    [INDX,aij,bij,res,criteria,GAoutput] = gaPCARX(Y,X,Q,[na 0 na],p,1,options);
    rss_sss(na) = criteria.rss_sss;
    connum(na) = criteria.connum;
    spp(na) = criteria.spp;
    sigma = var(res);
    bic(na) = log(sigma) + ((na+na+1)*B)*log(N*NoS)/(N*NoS);
    basisindx{na} = INDX{1};
    save('GAPCARXresults','rss_sss','bic','spp','connum','basisindx');
end



%% PC-ARX(43) Identification
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
N = 200; 
% Number of samples
NoS = 500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

[Y,X] = deal(zeros(NoS,N));
for i = 1:N
    out = load(['SpaceFrameRes',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('InputVars.mat','Q');


cd('C:\Users\sminas\Dropbox\MATLAB\PCE\m-files\');
load('GAPCARXresults','rss_sss','bic','spp','connum','basisindx');
q = 1;
na = 43;
indx{1} = basisindx{43};
indx{2} = indx{1};
B = length(indx{1});


options.maxsize = 1e9;
options.basis = 'legen';
options.focus = 'prediction';
[aij0,bij0,res0,criteria0] = pcarx(Y,X,Q,[na 0 na],indx,[],options);
options.focus = 'simulation';
[aij,bij,res,criteria] = pcarx(Y,X,Q,[na 0 na],indx,[aij0;bij0],options);


%% PC-ARX(50) Identification
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM';
cd(local_dir)

write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
N = 200; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

[Y,X] = deal(zeros(NoS,N));
for i = 1:N
    out = load(['SpaceFrameRes',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('InputVars.mat','Q');


cd('C:\Users\sminas\Dropbox\MATLAB\PCE\m-files\');
load('PCARXresults','rss_sss','bic','spp','connum');
na = 50;
p = 2;
q = 1;
indx{1} = combinations(M,p,q);
indx{2} = indx{1};
B = length(indx{1});

options.maxsize = 1e9;
options.basis = 'legen';
options.focus = 'prediction';
[aij0,bij0,res0,criteria0] = pcarx(Y,X,Q,[na 0 na],indx,[],options);
options.focus = 'simulation';
[aij,bij,res,criteria] = pcarx(Y,X,Q,[na 0 na],indx,[aij0;bij0],options);


%%
res0 = reshape(res0,500-53,200);
res = reshape(res,500,200);

figure(1)
ccf(X(:,1).^2,Y(:,1)-mean(Y(:,1)),100,0.8);
figure(2)
ccf(Y(:,1)-mean(Y(:,1)),(Y(:,1)-mean(Y(:,1))).^2,100,0.8);
% Cross-corellation(z,z^2)
for ii=1:20
    figure(ii)
    ccf(res(:,ii)-mean(res(:,ii)),(res(:,ii)-mean(res(:,ii))).^2,100,0.8);
end
figure(4)
ccf(NLres-mean(NLres),(NLres-mean(NLres)).^2,100,0.8);

%% Parameter surfaces
print_to_eps = 'y';
clear Qsurf
NoI1 = 100;
NoI2 = 100;
aux = linspace(-1,1,NoI1);
Qsurf(1,:) = kron(ones(1,NoI2),aux);
Qsurf(3,:) = kron(aux,ones(1,NoI2));
Qsurf(2,:) = zeros(1,NoI1*NoI2);
Qsurf(4,:) = zeros(1,NoI1*NoI2);

[~,an_interp] = PCparam(Qsurf,4,indx{1},aij,options);
[~,bn_interp] = PCparam(Qsurf,4,indx{1},bij,options);

figure(2)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
for i = 1:4
    subplot(2,4,i),plot3(Qsurf(1,:),Qsurf(2,:),an_interp(:,i),'.')
    xlabel('E_x')
    ylabel('Yield stress','Fontsize',14)
    zlabel(['$a_',num2str(i),'(\xi)$'],'Interpreter','Latex','Fontsize',14)
    grid on
end

       
figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
for i = 1:4
    subplot(2,4,i),surf(linspace(0.8,1.2,NoI2),linspace(0.8,1.2,NoI1),reshape(an_interp(:,i),NoI1,NoI2))
    shading interp
    set(gca,'Fontsize',6)
    xlim([0.8 1.2])
    ylim([0.8 1.2])
    xlabel('$h$ (\%)','Interpreter','Latex','Fontsize',7)
    ylabel('$E_x$ (\%)','Interpreter','Latex','Fontsize',7)
    zlabel(['$a_',num2str(i),'(\xi)$'],'Interpreter','Latex','Fontsize',7)
    grid on
    box on
end
% for i = 1:4
%     subplot(2,4,i+4),surf(linspace(0.8,1.2,NoI2),linspace(0.8,1.2,NoI1),reshape(bn_interp(:,i),NoI1,NoI2))
%     shading interp
%     set(gca,'Fontsize',6)
%     xlim([0.8 1.2])
%     ylim([0.8 1.2])
%     xlabel('$h$ (\%)','Interpreter','Latex','Fontsize',7)
%     ylabel('$E_x$ (\%)','Interpreter','Latex','Fontsize',7)
%     zlabel(['$b_',num2str(i),'(\xi)$'],'Interpreter','Latex','Fontsize',7)
%     grid on
%     box on
% end
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'ARXsurf']),close;
     result = eps2xxx([write_dir,'ARXsurf.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% Parameter curves
print_to_eps = 'y';
clear Qsurf NoI*
NoI1 = 1000;
Qsurf(3,:) = linspace(-1,1,NoI1);
Qsurf(2,:) = zeros(1,NoI1);
Qsurf(1,:) = zeros(1,NoI1);
Qsurf(4,:) = zeros(1,NoI1);

[~,an_interp] = PCparam(Qsurf,4,indx{1},aij,options);
[~,bn_interp] = PCparam(Qsurf,4,indx{1},bij,options);

figure(2)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
for i = 1:4
    subplot(2,4,i),plot(linspace(0.8,1.2,NoI1),an_interp(:,i),'-')
    xlim([0.8 1.2])
    set(gca,'Fontsize',6)
    xlabel('$h$ (\%)','Interpreter','Latex','Fontsize',12)
    ylabel(['$a_',num2str(i),'(\xi)$'],'Interpreter','Latex','Fontsize',12)
    grid on
end
% for i = 1:4
%     subplot(2,4,i+4),plot(linspace(0.8,1.2,NoI1),bn_interp(:,i),'-')
%     xlim([0.8 1.2])
%     set(gca,'Fontsize',6)
%     xlabel('$h$ (\%)','Interpreter','Latex')
%     ylabel(['$b_',num2str(i),'(\xi)$'],'Interpreter','Latex','Fontsize',14)
%     grid on
% end
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'ARXcurve']),close;
     result = eps2xxx([write_dir,'ARXcurve.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end




%% Predictions
print_to_eps = 'y';

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
    subplot(5,1,1:3),plot(time(1:1000),Y(:,2),'-b')
    hold on,plot(time(44:1000),Y(44:1000,2)-RES(:,2),'--r')
    
    %xlim([0.8 1.2])
    set(gca,'Fontsize',9)
    xlabel('Time (s)','Fontsize',12)
    ylabel('\alpha_{54} (m/s^2)','Fontsize',12)
    grid on
legend({'Simulation run #2','PC-ARX(43,43) predictions'})
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'preds']),close;
     result = eps2xxx([write_dir,'preds.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end






%% GA - NARX
clear,clc,close all
local_dir = 'J:\Documents\MATLABroutines\PCNARXjournal\SpaceFrameFEM\';
cd(local_dir)

% Number of experiments
N = 200; 
% Number of samples
NoS = 500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

out = load(['SpaceFrameRes169.txt']);
Y = out(1:NoS,13);
X = out(1:NoS,2);

cd('J:\Documents\MATLABroutines\PCE\m-files')
clc, clear *rss*
options.focus = 'prediction';
options.GA.PopulationSize = 100 ;
options.GA.PlotFcns = {@gaplotbestf,@gaplotbestindiv};
options.GA.Display = 'iter';
P = 4;
q = 0.4;
maxterms = 4;

disp('Linear    Nonlinear')
disp('===================')
for ii = 20
    na = ii;
    [theta,res,criteria,regstr] = narx(Y,X,[na 0 na],[],1,1,[],[],options);
    rss_sss(ii) = criteria.rss_sss;
    sigma(ii) = var(res);
    bic(ii) = log(sigma(ii)) + (2*na+1)*log(N)/N;
    [regstr,regindx{ii},theta,NLres,NLcriteria] = gaNARX(Y,X,[na 0 na],[],P,q,maxterms,options);    
    NLrss_sss(ii) = NLcriteria.rss_sss;
    NLsigma(ii) = var(NLres);
    NLbic(ii) = log(NLsigma(ii)) + (2*na+1)*log(N)/N;
    s = sprintf('%2.5f \t %2.5f',rss_sss(ii),NLrss_sss(ii));
    disp(s)
end

[~,~,~,InitRegstr] = narx(Y,X,[na 0 na],[],P,q,maxterms,[],options);
TotIndx = [];
for ii=1:10
    TotIndx = [TotIndx  regindx{ii}];
end
hist(TotIndx,length(InitRegstr))

% Cross-corellation(x^2,y)
figure(1)
ccf(X.^2,Y-mean(Y),100,0.8);
figure(2)
ccf(Y-mean(Y),(Y-mean(Y)).^2,100,0.8);
% Cross-corellation(z,z^2)
figure(3)
ccf(res-mean(res),(res-mean(res)).^2,100,0.8);
figure(4)
ccf(NLres-mean(NLres),(NLres-mean(NLres)).^2,100,0.8);


%% Comparison with conventional AR models 
clear,clc,close all
local_dir = 'J:\Documents\MATLABroutines\PCNARXjournal\SpaceFrameFEM\';
cd(local_dir)

% Number of experiments
N = 200; 
% Number of samples
NoS = 2500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 4; 
% Sampling period
Ts =  0.025;

out = load('SpaceFrameRes169.txt');
Y = out(1:NoS,9);
X = out(1:NoS,2);

% Conventional AR identification
for na = 10% :150
    IdData = iddata(Y,X,Ts);
    model = arx(IdData,na);
    AR(i,1:na)= model.a(2:end); 
end

Qtrans = load([local_dir,'ANSYS\ShearFrameNL\MProps.txt'],'-ascii');
Q = [(Qtrans(:,2)-2e11)/1e10 ,(Qtrans(:,4)-0.065)/0.025 ,(Qtrans(:,6)-400e6)/200e6]';
y = zeros(NoS,N);

for i = 1:N
    out = load([local_dir,'ANSYS\ShearFrameNL\SFres',num2str(i),'.txt']);
    y(:,i) = out(1:NoS,9);
end

options.maxsize = 1e7;
options.basis = 'legen';
na = 8;
p = 5;
% AR modelling
disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
indx = combinations(M,p,q);
[A,res,criteria] = pcear(y,Q,na,indx,options);
% PCE-AR parameters
[phi,an] = PCEARparam(Q,p,indx,A,options);

% Conventional AR identification
for i = 1:N
    disp(['Exp. #: ',num2str(i)]);
    IdData = iddata(y(:,i),[],1/128);
    model = arx(IdData,na);
    AR(i,1:na)= model.a(2:end); 
end

close all
figure(1),
plot(1:N,AR,'-o','Markersize',6,'linewidth',2)
grid on
ylabel('AR parameters')
xlabel('Experiment number')
legend({'\alpha_1','\alpha_2','\alpha_3','\alpha_4','\alpha_5','\alpha_6'})
hold on
plot(1:N,an,'--x','Markersize',6,'linewidth',2)



figure(2)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(5,1,1:4),hold on
symb = ['o','x','d','s','>','^','*','h'];
clr = colormap(lines(8));
for i=1:na
    plot(1:N,an(:,i),'-','Markersize',5,'linewidth',1,'Marker',symb(i),'Color',clr(i,:))
end
legend({'$a_{1,n}$','$a_{2,n}$','$a_{3,n}$','$a_{4,n}$','$a_{5,n}$','$a_{6,n}$','$a_{7,n}$','$a_{8,n}$'},...
    'Location','NorthOutside','Orientation','Horizontal','Interpreter','Latex')
xlim([0.5 50.5])
grid on
box on
set(gca,'Fontname','TimesNewRoman','Fontsize',9)
ylabel('AR parameters','Fontname','TimesNewRoman','Fontsize',12)
xlabel('Experiment number ($n$)','Fontname','TimesNewRoman','Fontsize',12)
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'ARparams']),close;
     result = eps2xxx([write_dir,'ARparams.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin32c.exe');
end


%% Parameter surfaces
clear,clc,close all
local_dir = 'C:\Users\sminas\Desktop\PCE\';
cd('C:\Users\sminas\Desktop\PCE\m-files')

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

N = 50;
M = 3;
NoS = 640;
q = 1;

Qtrans = load([local_dir,'ANSYS\ShearFrameNL\MProps.txt'],'-ascii');
Q = [(Qtrans(:,2)-2e11)/1e10 ,(Qtrans(:,4)-0.065)/0.025 ,(Qtrans(:,6)-400e6)/200e6]';
y = zeros(NoS,N);

for i = 1:N
    out = load([local_dir,'ANSYS\ShearFrameNL\SFres',num2str(i),'.txt']);
    y(:,i) = out(1:NoS,9);
end

options.maxsize = 1e7;
options.basis = 'legen';
na = 8;
p = 5;
% AR modelling
disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
indx = combinations(M,p,q);
[A,res,criteria] = pcear(y,Q,na,indx,options);
[phi,an] = PCEARparam(Q,p,indx,A,options);

clear Qsurf
NoI1 = 100;
NoI2 = 166;
aux = linspace(-1,1,NoI1);
Qsurf(1,:) = kron(ones(1,NoI2),aux);
Qsurf(3,:) = kron(aux,ones(1,NoI2));
Qsurf(2,:) = zeros(1,NoI1*NoI2);
[phi,an_interp] = PCEARparam(Qsurf,p,indx,A,options);

figure(2)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
for i = 1:na
    subplot(2,4,i),plot3(Qsurf(1,:),Qsurf(2,:),an_interp(:,i),'.')
    xlabel('E_x')
    ylabel('Yield stress','Fontsize',14)
    zlabel(['$a_',num2str(i),'(\xi)$'],'Interpreter','Latex','Fontsize',14)
    grid on
end

       
figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
for i = 1:na
    subplot(2,4,i),surf(linspace(-1,1,NoI2),linspace(-1,1,NoI1),reshape(an_interp(:,i),NoI1,NoI2))
    shading interp
    set(gca,'Fontsize',6)
    xlabel('$\sigma_y$ (Pa)','Interpreter','Latex','Fontsize',6)
    ylabel('$E_x$ (Pa)','Interpreter','Latex','Fontsize',6)
    zlabel(['$a_',num2str(i),'(\xi)$'],'Interpreter','Latex','Fontsize',6)
    grid on
    box on
end
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'ARsurf']),close;
     result = eps2xxx([write_dir,'ARsurf.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin32c.exe');
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
% 4     Displacement of 1st floor (node 18)
% 5     Displacement of 2nd floor (node 27)
% 6     Displacement of 3rd floor (node 36)
% 7     Displacement of 4th floor (node 45)
% 8     Displacement of 5th floor (node 54)
% 9     Acceleration response of 1st floor (node 18)
% 10    Acceleration response of 2nd floor (node 27)
% 11    Acceleration response of 3rd floor (node 36)
% 12    Acceleration response of 4th floor (node 45)
% 13    Acceleration response of 5th floor (node 54)