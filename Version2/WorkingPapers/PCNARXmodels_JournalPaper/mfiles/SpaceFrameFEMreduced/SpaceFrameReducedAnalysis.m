%% Space frame modal analysis (linear material properties)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced\')
% Number of experiments
N = 50; 
% Number of variables (Density, Ex = Ey = Ez, Prxy = Prxz = Prxyz)
M = 1; 

% Latin Hypercube Sampling
Q = lhsdesign(N,M,'criterion','maximin','iterations',10);
Q = 2*(Q'-0.5);
Qtrans = [2e11 + 1e10*Q(1,:)]';
save('BeamProps.txt','Qtrans','-ascii')

% Run Ansys
delete('file.*')
% dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\SpaceFrameModal.f" -o "output.txt"');


%% Space frame modal analysis (linear material properties) -- Figure
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced\')
load('Wn.txt')
figure
subplot(211),plot(Wn(:,1:5),'-o')
xlabel('Experiment number')
ylabel('Natural Frequencies (Hz)')
grid on
subplot(212),plot(100*Wn(:,10:18),'-o')
xlabel('Experiment number')
ylabel('Damping ratios (%)')
grid on


%% Space frame -- bilinear isotropic material (increasing force)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced\')

NoS = 2000;
NoExp = 1;
% Generation of the input force
Finput = linspace(5e3,5e4,NoS)';                                        % Linear
save('InputForce.txt','Finput','-ascii')

% Run Ansys
delete('file.*')
% dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEM\SpaceFrameForce.f" -o "output.txt"');


%% Space frame -- bilinear isotropic material (increasing force) -- Figure
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced\')
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



%% Transient nonlinear analysis (various acceleration levels)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced\')

% Number of experiments (estimation set)
K1 = 40; 
% Number of experiments (validation set)
K2 = 40; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, std_a)
M = 2; 

% sigma_a distribution
PDFsigma = ProbDistUnivParam('lognormal',[2 0.5]);
xindx = 0:0.001:20;     Fx = pdf(PDFsigma,xindx);

% Latin Hypercube Sampling (estimation set)
rng(100)
Q_est = lhsdesign(K1,M,'criterion','maximin','iterations',10);
Q_est = 2*(Q_est'-0.5);
Qtrans_est = 1 + 0.05*Q_est(1,:)';
sigma_est = icdf('lognormal',0.5*Q_est(2,:)+0.5,2,0.5);

% Latin Hypercube Sampling (validation set)
rng(200)
Q_val = lhsdesign(K2,M,'criterion','maximin','iterations',10);
Q_val = 2*(Q_val'-0.5);
Qtrans_val = 1 + 0.05*Q_val(1,:)';
sigma_val = icdf('lognormal',0.5*Q_val(2,:)+0.5,2,0.5);

% Beam properties file
BeamProps = [Qtrans_est(:) ; Qtrans_val(:)];
save('BeamProps.txt','BeamProps','-ascii')

% Generation of the random signal excitation
for k=1:K1+K2 
    rng(k)
    if k<=K1
        Exc((k-1)*NoS+1:k*NoS,1) = sigma_est(k)*randn(NoS,1);
    else
        Exc((k-1)*NoS+1:k*NoS,1) = sigma_val(k-K1)*randn(NoS,1);
    end
end
save('AccRandom.txt','Exc','-ascii')

UNCpars_est = [Qtrans_est sigma_est'];
UNCpars_val = [Qtrans_val sigma_val'];
Exc = reshape(Exc,NoS,K1+K2);
save('InputVars.mat','Q_est','Q_val','Exc','UNCpars*')

% Run Ansys
delete('file.*')
% dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "D:\SpaceFrame\SpaceFrameReduced.f" -o "output.txt"');


%% Transient nonlinear analysis (various acceleration levels - normal pdf)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced\')

% Number of experiments (estimation set)
K1 = 40; 
% Number of experiments (validation set)
K2 = 40; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, std_a)
M = 2; 

% Latin Hypercube Sampling (estimation set)
rng(100)
Q_est = lhsnorm(zeros(M,1),eye(M),K1);
Qtrans_est = 1 + 0.02*Q_est(:,1);
sigma_est = exp(2 + 0.5*Q_est(:,2));
% Latin Hypercube Sampling (validation set)
rng(200)
Q_val = lhsnorm(zeros(M,1),eye(M),K2);
Qtrans_val = 1 + 0.02*Q_val(:,1);
sigma_val = exp(2+ 0.5*Q_val(:,2));
% Beam properties file
BeamProps = [Qtrans_est(:) ;Qtrans_val(:)];
save('BeamProps.txt','BeamProps','-ascii')

% Generation of the random signal excitation
for k=1:K1+K2 
    rng(k)
    if k<=K1
        Exc((k-1)*NoS+1:k*NoS,1) = sigma_est(k)*randn(NoS,1);
    else
        Exc((k-1)*NoS+1:k*NoS,1) = sigma_val(k-K1)*randn(NoS,1);
    end
end
save('AccRandom.txt','Exc','-ascii')

UNCpars_est = [Qtrans_est sigma_est];
UNCpars_val = [Qtrans_val sigma_val];
Exc = reshape(Exc,NoS,K1+K2);
save('InputVars.mat','Q_est','Q_val','Exc','UNCpars*')

% Run Ansys
delete('file.*')
% dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "D:\SpaceFrame\SpaceFrameReduced.f" -o "output.txt"');


%% Transient nonlinear analysis  -- Validation set (equispaced grid)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced\')

% Number of experiments
K = 25; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, std_a)
M = 2;

% Equispaced grid
Q(:,1) = kron([-1 -0.5 0 0.5 1],ones(1,5));     % knl
Q(:,2) = kron(ones(1,5),[-0.99 -0.5 0 0.5 0.99]);     % Fmax
Qtrans = 1 + 0.05*Q(:,1)';
BeamProps = Qtrans(:);
save('BeamProps.txt','BeamProps','-ascii')

PDFsigma = ProbDistUnivParam('lognormal',[0.6 0.4]);
xindx = 0:0.001:20;     Fx = pdf(PDFsigma,xindx);
sigma = icdf('lognormal',0.5*Q(:,2) +0.5,0.6,0.4);

% Generation of the random signal excitation
for k=1:K
    rng(k+100)
    Exc((k-1)*NoS+1:k*NoS,1) = sigma(k)*randn(NoS,1);
end
save('AccRandom.txt','Exc','-ascii')

UNCpars = [Qtrans sigma'];
Exc = reshape(Exc,NoS,K);
save('InputVarsValidationSet.mat','Q','Exc','UNCpars')

% Run Ansys
delete('file.*')
dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "D:\SpaceFrame\SpaceFrameReduced.f" -o "output.txt"');



%% Transient nonlinear analysis (various acceleration levels) -- Figures
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced\Responses\EstimationSet')
Dt = 0.02;
Fs = 1/Dt;
K = 25;
figure(1),
clr = colormap(lines(3));
subplot(311),hold on
subplot(312),hold on
subplot(313),hold on

for k = 1:K
    SFRes = load(['SFres',num2str(k),'.txt']);
    time = SFRes(:,1);
    acci = SFRes(:,2);
    sumreact = SFRes(:,3);
    ux = SFRes(:,4:8);
    veloc =  SFRes(:,9:13);
    acc =  SFRes(:,14:18);
    [Py(:,k),F] = pwelch(veloc(:,1),512,480,512,Fs);
    [Txy(:,k)] = tfestimate(acci,veloc(:,1),512,480,512,Fs);
    [Cxy(:,k),F] = mscohere(acci,veloc(:,1),512,480,512,Fs);
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
%  axis([time(1) time(end) -12 12])
subplot(312),ylabel('\Sigma(Reaction forces) (kN)')
set(gca,'Xticklabel',[])
% axis([time(1) time(end) -3000 3000])
subplot(313),ylabel('Top floor displacement (m)')
xlabel('Time (s)')
% axis([time(1) time(end) -0.8 0.8])


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
surf(1:25,F,20*log10(abs(Txy)))
shading interp
zlabel('Magnitude (dB)')
ylabel('Frequency (Hz)')
xlabel('Experiment No.')
title('FRF')
axis tight

figure(5)
surf(1:25,F,Cxy)
shading interp
zlabel('Coherence function')
ylabel('Frequency (Hz)')
xlabel('Experiment No.')
axis tight



%% NARX regressors selection (GA - rss)
% =========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced';
end
cd(local_dir)

% Number of experiments
K = 40; 
% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;

% Simulation experiment
SimExp = 25;
out = load(['Responses/SFres',num2str(SimExp),'.txt']);
Y = out(1:N,9); % 13: acc, 8: velocity
X = out(1:N,2);

% Model orders 
na = 10;     nb = 10;     nd = 0;
% Polynomial regressors
P = 1;      maxterms = 1;
% regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))'});
regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))',...
     'abs(Y(t-2,:))','abs(Y(t-3,:))','abs(Y(t-4,:))','abs(Y(t-5,:))'});
%      ,...
%      'abs(Y(t-6,:))','abs(Y(t-7,:))','abs(Y(t-8,:))','abs(Y(t-9,:))','abs(Y(t-10,:))'});

% Optimization options
options.warnings = 'off';
options.GA.PopulationSize = 100;
options.maxsize = 1e8;

rng(1);
options.focus = 'simulation';
options.criterion = 'rss';
options.parfor = 'n';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);
options.focus = 'prediction';
options.Nr = length(GAreg);
[SEstructure,IndxRem] = SEreduction(Y,X,[na nd nb],GAreg,options);
save('SF_RandomExc_GA_rss.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')


%% NARX regressors selection Gaussian (GA - rss)
% =========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced';
end
cd(local_dir)

% Number of experiments
K = 40; 
% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;

% Simulation experiment
SimExp = 33;
out = load(['ResponsesGaussian/SFres',num2str(SimExp),'.txt']);
Y = out(1:N,9); % 13: acc, 8: velocity
X = out(1:N,2);

% Model orders 
na = 10;     nb = 10;     nd = 0;
% Polynomial regressors
P = 1;      maxterms = 1;
% regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))'});
regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))',...
     'abs(Y(t-2,:))','abs(Y(t-3,:))','abs(Y(t-4,:))','abs(Y(t-5,:))'});
%      ,...
%      'abs(Y(t-6,:))','abs(Y(t-7,:))','abs(Y(t-8,:))','abs(Y(t-9,:))','abs(Y(t-10,:))'});

% Optimization options
options.warnings = 'off';
options.GA.PopulationSize = 100;
options.maxsize = 1e8;

rng(1);
options.focus = 'simulation';
options.criterion = 'rss';
options.parfor = 'n';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);
options.focus = 'prediction';
options.Nr = length(GAreg);
[SEstructure,IndxRem] = SEreduction(Y,X,[na nd nb],GAreg,options);
save('SF_RandomExc_GA_normal_rss.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')



%% NARX regressors selection (GA - mnse)
% =========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced';
end
cd(local_dir)

% Number of experiments
K = 40; 
% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;

% Simulation experiment
SimExp = 13;
out = load(['Responses/SFres',num2str(SimExp),'.txt']);
Y = out(1:N,8); % 13: acc, 8: velocity
X = out(1:N,2);

% Model orders 
na = 10;     nb = 10;     nd = 0;
% Polynomial regressors
P = 1;      maxterms = 1;
% regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))'});
regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))',...
     'abs(Y(t-2,:))','abs(Y(t-3,:))','abs(Y(t-4,:))','abs(Y(t-5,:))',...
     'abs(Y(t-6,:))','abs(Y(t-7,:))','abs(Y(t-8,:))','abs(Y(t-9,:))','abs(Y(t-10,:))'});

% Optimization options
options.warnings = 'off';
options.GA.PopulationSize = 1000;
options.maxsize = 1e8;

rng(100);
options.focus = 'simulation';
options.criterion = 'mnse';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);
options.focus = 'prediction';
options.Nr = length(GAreg);
[SEstructure,IndxRem] = SEreduction(Y,X,[na nd nb],GAreg,options);
save('SF_RandomExc_GA_mnse.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')


%% NARX regressors selection -- Figure
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced')
end
load('SF_RandomExc_GA_rss.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

NoReg = length(GAreg);
close all

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(2,1,1),plot(1:NoReg,SEstructure.values(1:NoReg),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 NoReg],SEstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
set(gca,'yscale','linear','Fontsize',7,'xtick',[1:10:NoReg],'xticklabel',[1:10:NoReg])
% for ii = 1:NoReg
%     if strcmp(SEstructure.removed{ii},'ones(length(t),size(Y,2))')
%         txtstr = '$\mbox{const. term}$'; 
%     elseif length(SEstructure.removed{ii})>13
%         if strcmp(SEstructure.removed{ii}(12),'1');
%             txtstr = ['$',lower(SEstructure.removed{ii}(2)),'[',SEstructure.removed{ii}(4:6),']','\cdot |y[t-1]|$'];
%         else
%             txtstr = ['$',lower(SEstructure.removed{ii}(2)),SEstructure.removed{ii}(11:12),'[',SEstructure.removed{ii}(4:6),']','\cdot |y[t-1]|$'];
%         end
%     else
%         if strcmp(SEstructure.removed{ii}(12),'1');
%             txtstr = ['$',lower(SEstructure.removed{ii}(2)),'[',SEstructure.removed{ii}(4:6),']','$'];
%         else
%             txtstr = ['$',lower(SEstructure.removed{ii}(2)),SEstructure.removed{ii}(11:12),'[',SEstructure.removed{ii}(4:6),']','$'];
%         end
%     end
%     text(ii,0.075,txtstr,'Rotation',65,'Fontsize',7,'Horizontalalignment','right','Interpreter','Latex','Fontsize',9)
% end
% axis([0.75 NoReg+0.25 0.1 0.9])
grid on
text(NoReg/2,-0.35,'Regressors dropped','Fontangle','normal','Fontsize',9,'Horizontalalignment','center')
ylabel('MNSE','Fontangle','normal','Fontsize',9)
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'SF_MSS_StageB']),close;
     result = eps2xxx([write_dir,'SF_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% Local NARX model estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels';
end
cd(local_dir)
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 40;

% Input-Output data
[Y,X] = deal(zeros(NoS,K));
for i = 1:K
    out = load(['Responses/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,9);
    X(:,i) = out(1:NoS,2);
%      ReF(:,i) = out(1:NoS,3);
end
load('Responses/InputVars.mat','Q_est');

% Regressors
SelReg = {'(Y(t-3,:).^1).*abs(Y(t-2,:))';
    '(X(t-10,:).^1).*abs(Y(t-1,:))';
    '(X(t-10,:).^1).*abs(Y(t-2,:))';
    '(Y(t-8,:).^1).*abs(Y(t-2,:))';
    '(Y(t-9,:).^1).*abs(Y(t-2,:))';
    '(X(t-5,:).^1).*abs(Y(t-4,:))';
    '(Y(t-7,:).^1).*abs(Y(t-2,:))';
    '(X(t-5,:).^1).*abs(Y(t-2,:))';
    '(X(t-1,:).^1).*abs(Y(t-3,:))';
    '(X(t-0,:).^1).*abs(Y(t-2,:))';
    '(Y(t-6,:).^1).*abs(Y(t-3,:))';
    '(X(t-3,:).^1).*abs(Y(t-3,:))';
    '(X(t-0,:).^1).*abs(Y(t-1,:))';
    '(Y(t-9,:).^1).*abs(Y(t-5,:))';
    '(Y(t-10,:).^1).*abs(Y(t-5,:))';
    '(X(t-6,:).^1).*abs(Y(t-2,:))';
    '(X(t-6,:).^1).*abs(Y(t-3,:))';
    '(X(t-6,:).^1).*abs(Y(t-1,:))';
    '(X(t-8,:).^1).*abs(Y(t-3,:))';
    '(X(t-4,:).^1).*abs(Y(t-3,:))';
    '(X(t-7,:).^1)';
    '(X(t-3,:).^1).*abs(Y(t-2,:))';
    '(Y(t-3,:).^1).*abs(Y(t-4,:))';
    '(Y(t-7,:).^1).*abs(Y(t-1,:))';
    '(Y(t-5,:).^1).*abs(Y(t-1,:))';
    '(X(t-9,:).^1).*abs(Y(t-4,:))';
    '(Y(t-5,:).^1).*abs(Y(t-3,:))';
    '(Y(t-7,:).^1).*abs(Y(t-3,:))';
    '(X(t-0,:).^1).*abs(Y(t-3,:))';
    '(X(t-6,:).^1).*abs(Y(t-4,:))';
    '(X(t-7,:).^1).*abs(Y(t-5,:))';
    '(X(t-8,:).^1).*abs(Y(t-4,:))';
    '(X(t-7,:).^1).*abs(Y(t-4,:))';
    '(Y(t-4,:).^1).*abs(Y(t-4,:))';
    '(Y(t-5,:).^1).*abs(Y(t-5,:))';
    '(Y(t-8,:).^1).*abs(Y(t-1,:))';
    '(Y(t-1,:).^1).*abs(Y(t-3,:))';
    '(Y(t-3,:).^1).*abs(Y(t-3,:))';
    '(X(t-8,:).^1).*abs(Y(t-5,:))';
    '(Y(t-7,:).^1).*abs(Y(t-5,:))';
    '(Y(t-3,:).^1).*abs(Y(t-5,:))';
    '(X(t-2,:).^1).*abs(Y(t-5,:))';
    '(Y(t-2,:).^1).*abs(Y(t-5,:))';
    '(X(t-5,:).^1).*abs(Y(t-5,:))';
    '(Y(t-4,:).^1).*abs(Y(t-5,:))';
    '(Y(t-10,:).^1)';
    '(X(t-10,:).^1)';
    '(Y(t-8,:).^1)';
    '(X(t-9,:).^1)';
    '(X(t-6,:).^1)';
    '(X(t-2,:).^1)';
    '(X(t-5,:).^1)';
    '(X(t-1,:).^1)';
    '(Y(t-5,:).^1)';
    '(X(t-4,:).^1)';
    '(X(t-3,:).^1)';
    '(Y(t-7,:).^1)';
    '(Y(t-9,:).^1)';
    '(Y(t-4,:).^1)';
    '(Y(t-1,:).^1)';
    '(Y(t-3,:).^1)';
    '(Y(t-2,:).^1)';
    '(X(t-0,:).^1)'};

% Model orders 
na = 10;     nb = 10;     nd = 0;
orders = [na nd nb];

clear *criteria* *res* theta*

% PE method
for k = 1:K
    options.focus = 'prediction';
    options.method = 'ols';
    [thetaPE(:,k),NLresPE{k},PEcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
end

% SE method
options.focus = 'simulation';
options.PEinit = 'y';
options.nlmethod = 'LM';
matlabpool('open',4)
parfor k = 1:K
    disp(k)
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
matlabpool('close')

save('SF_RandomExc_LocalNARX_lsqnonlin.mat','*criteria','NLres*','theta*','SelReg','orders')



%% Local NARX model estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels';
end
cd(local_dir)
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 40;

% Input-Output data
[Y,X] = deal(zeros(NoS,K));
for i = 1:K
    out = load(['ResponsesGaussian/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,9);
    X(:,i) = out(1:NoS,2);
%      ReF(:,i) = out(1:NoS,3);
end
load('ResponsesGaussian/InputVars.mat','Q_est');

% Regressors
SelReg = {   '(X(t-1,:).^1).*abs(Y(t-5,:))';
    '(Y(t-3,:).^1).*abs(Y(t-2,:))';
    '(X(t-9,:).^1).*abs(Y(t-5,:))';
    '(X(t-3,:).^1).*abs(Y(t-2,:))';
    '(X(t-10,:).^1).*abs(Y(t-5,:))';
    '(X(t-8,:).^1).*abs(Y(t-4,:))';
    '(Y(t-3,:).^1).*abs(Y(t-4,:))';
    '(Y(t-4,:).^1).*abs(Y(t-5,:))';
    '(Y(t-5,:).^1).*abs(Y(t-5,:))';
    '(Y(t-10,:).^1).*abs(Y(t-5,:))';
    '(X(t-9,:).^1).*abs(Y(t-2,:))';
    '(X(t-5,:).^1).*abs(Y(t-4,:))';
    '(Y(t-7,:).^1)';
    '(X(t-7,:).^1).*abs(Y(t-3,:))';
    '(X(t-2,:).^1).*abs(Y(t-4,:))';
    '(X(t-7,:).^1).*abs(Y(t-2,:))';
    '(X(t-5,:).^1).*abs(Y(t-1,:))';
    '(Y(t-5,:).^1).*abs(Y(t-2,:))';
    '(X(t-1,:).^1).*abs(Y(t-3,:))';
    '(X(t-6,:).^1).*abs(Y(t-3,:))';
    '(Y(t-6,:).^1).*abs(Y(t-3,:))';
    '(X(t-8,:).^1).*abs(Y(t-3,:))';
    '(Y(t-6,:).^1).*abs(Y(t-5,:))';
    '(X(t-8,:).^1).*abs(Y(t-5,:))';
    '(X(t-0,:).^1).*abs(Y(t-1,:))';
    '(Y(t-7,:).^1).*abs(Y(t-2,:))';
    '(Y(t-2,:).^1).*abs(Y(t-2,:))';
    '(X(t-4,:).^1).*abs(Y(t-5,:))';
    '(X(t-2,:).^1).*abs(Y(t-3,:))';
    '(X(t-3,:).^1).*abs(Y(t-4,:))';
    '(X(t-2,:).^1).*abs(Y(t-1,:))';
    '(X(t-5,:).^1).*abs(Y(t-3,:))';
    '(X(t-1,:).^1).*abs(Y(t-4,:))';
    '(Y(t-7,:).^1).*abs(Y(t-4,:))';
    '(Y(t-7,:).^1).*abs(Y(t-3,:))';
    '(X(t-5,:).^1).*abs(Y(t-5,:))';
    '(X(t-0,:).^1).*abs(Y(t-3,:))';
    '(X(t-1,:).^1).*abs(Y(t-1,:))';
    '(X(t-10,:).^1).*abs(Y(t-2,:))';
    '(X(t-10,:).^1).*abs(Y(t-4,:))';
    '(Y(t-2,:).^1).*abs(Y(t-5,:))';
    '(Y(t-5,:).^1).*abs(Y(t-4,:))';
    '(X(t-9,:).^1).*abs(Y(t-1,:))';
    '(Y(t-10,:).^1).*abs(Y(t-1,:))';
    '(Y(t-4,:).^1).*abs(Y(t-2,:))';
    '(Y(t-3,:).^1).*abs(Y(t-1,:))';
    '(X(t-3,:).^1).*abs(Y(t-3,:))';
    '(Y(t-6,:).^1).*abs(Y(t-4,:))';
    '(X(t-7,:).^1).*abs(Y(t-4,:))';
    '(Y(t-3,:).^1)';
    '(Y(t-4,:).^1)';
    '(X(t-6,:).^1)';
    '(X(t-4,:).^1).*abs(Y(t-4,:))';
    '(X(t-1,:).^1)';
    '(X(t-6,:).^1).*abs(Y(t-4,:))';
    '(X(t-8,:).^1)';
    '(X(t-4,:).^1)';
    '(X(t-9,:).^1)';
    '(X(t-5,:).^1)';
    '(Y(t-6,:).^1)';
    '(X(t-3,:).^1)';
    '(X(t-2,:).^1)';
    '(Y(t-9,:).^1)';
    '(Y(t-1,:).^1)';
    '(Y(t-5,:).^1)';
    '(Y(t-2,:).^1)';
    '(X(t-0,:).^1)'};

% Model orders 
na = 10;     nb = 10;     nd = 0;
orders = [na nd nb];

clear *criteria* *res* theta*

% PE method
for k = 1:K
    options.focus = 'prediction';
    options.method = 'ols';
    [thetaPE(:,k),NLresPE{k},PEcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
end

% SE method
options.focus = 'simulation';
options.PEinit = 'y';
options.nlmethod = 'LM';
matlabpool('open',4)
parfor k = 1:K
    disp(k)
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
matlabpool('close')

save('SF_RandomExc_LocalNARX_normal_lsqnonlin.mat','*criteria','NLres*','theta*','SelReg','orders')


%% Parameters expansion
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced/';
end
cd(local_dir)
load('SF_RandomExc_LocalNARX_lsqnonlin.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 40;
% Input-Output data
[Y,X] = deal(zeros(NoS,K));
for i = 1:K
    out = load(['Responses/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,9);
    X(:,i) = out(1:NoS,2);
end
load('Responses/InputVars.mat','Q_est');

clear basis options
options.basis = 'legen';
options.criterion = 'rss';
Pmax = 7;
for P = 0:Pmax
    % Basis index
    options.Nb = 0;
    [BASISstructure,INDnew] = PCEreduction(Q_est,P,thetaSIM,options);
    R2o(P+1)= BASISstructure.initvalue;
end

% Basis index
P = 3;
options.Nb = 9; Nb = 9;
[BASISstructure,INDnew] = PCEreduction(Q_est,P,thetaSIM,options);


figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(221),plot(0:Pmax,R2o,'-o','Markersize',6,'Linewidth',1.2)
grid on
ylabel('RSS/SSS (%)','Fontsize',11)
xlabel(sprintf('%s\n%s','Total PC basis degree','(complete subspace)'),'Fontsize',10)
set(gca,'yscale','linear','Fontsize',9,'xtick',[0:4])
subplot(222),plot(1:Nb,BASISstructure.values(1:Nb),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 Nb],BASISstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
ylabel('RSS/SSS (%)','Fontsize',11)
set(gca,'yscale','linear','Fontsize',9,'xtick',[1:Nb],'xticklabel',[])
for ii = 1:Nb 
    text(ii,-0.6e-4,['$[',num2str(BASISstructure.removed(ii,1)),',',num2str(BASISstructure.removed(ii,2)),']$'],'Rotation',45,...
        'Horizontalalignment','center','Interpreter','Latex','Fontsize',9)
end
xlim([0.99 9.01])
grid on
text(5,-1.8e-4,sprintf('%s\n%s','PC bases dropped','(multivariable indeces)'),'Fontsize',10,'HorizontalAlignment','center')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'SF_MSS_StageB']),close;
     if ispc
        result = eps2xxx([write_dir,'SF_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
        result = eps2xxx([write_dir,'SF_MSS_StageB.eps'],{'pdf'});
     end
end


%% Parameters expansion - Gaussian
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced/';
end
cd(local_dir)
load('SF_RandomExc_LocalNARX_normal_lsqnonlin.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 40;
% Input-Output data
[Y,X] = deal(zeros(NoS,K));
for i = 1:K
    out = load(['ResponsesGaussian/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,9);
    X(:,i) = out(1:NoS,2);
end
load('ResponsesGaussian/InputVars.mat','Q_est');

clear basis options
options.basis = 'hermi';
options.ortho = 'y';
options.criterion = 'bic';
Pmax = 7;
for P = 0:Pmax
    % Basis index
    options.Nb = 0;
    [BASISstructure,INDnew] = PCEreduction(Q_est',P,thetaSIM,options);
    R2o(P+1)= BASISstructure.initvalue;
end

% Basis index
P = 4;
options.Nb = 15; Nb = 15;
[BASISstructure,INDnew] = PCEreduction(Q_est',P,thetaSIM,options);


figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(221),plot(0:Pmax,R2o,'-o','Markersize',6,'Linewidth',1.2)
grid on
ylabel('RSS/SSS (%)','Fontsize',11)
xlabel(sprintf('%s\n%s','Total PC basis degree','(complete subspace)'),'Fontsize',10)
set(gca,'yscale','linear','Fontsize',9,'xtick',[0:4])
subplot(222),plot(1:Nb,BASISstructure.values(1:Nb),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 Nb],BASISstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
ylabel('RSS/SSS (%)','Fontsize',11)
set(gca,'yscale','linear','Fontsize',9,'xtick',[1:Nb],'xticklabel',[])
for ii = 1:Nb 
    text(ii,-0.6e-4,['$[',num2str(BASISstructure.removed(ii,1)),',',num2str(BASISstructure.removed(ii,2)),']$'],'Rotation',45,...
        'Horizontalalignment','center','Interpreter','Latex','Fontsize',9)
end
% xlim([0.99 9.01])
grid on
text(5,-1.8e-4,sprintf('%s\n%s','PC bases dropped','(multivariable indeces)'),'Fontsize',10,'HorizontalAlignment','center')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'SF_MSS_StageB']),close;
     if ispc
        result = eps2xxx([write_dir,'SF_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
        result = eps2xxx([write_dir,'SF_MSS_StageB.eps'],{'pdf'});
     end
end


%% PC-NARX estimation
% =========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced/';
end
cd(local_dir)
load('SF_RandomExc_LocalNARX_normal_lsqnonlin.mat')

% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 40;

% Input-Output data
[Y,X] = deal(zeros(NoS,K));
for i = 1:K
    out = load(['ResponsesGaussian/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,9);
    X(:,i) = out(1:NoS,2);
end
load('ResponsesGaussian/InputVars.mat','Q_est');

options.basis = 'hermi';

% PC basis index
INDX{1} =   [ 0     0;
     0     2;
     0     1];

% Total number of basis functions
B = size(INDX{1},1);
% Basis constraction: size(basis_i) = p x K
basis = cell(M,1);
% Total PC basis degree 
P = 2;
for m = 1:M
    basis{m} = PCbasis(Q_est(:,m)',P,options);
end

% Regressor matrix (tensor multiplication)
phiB = ones(B,K);
for i = 1:B
    for k = 1:K
        for m = 1:M
            phiB(i,k) = phiB(i,k)*basis{m}(INDX{1}(i,m)+1,k);
        end
    end
end
% Expansion coefficients OLS estimate
THij = (phiB')\(thetaSIM');
options.method = 'ols';
options.focus = 'simulation';
options.PEinit = 'y';
options.nlmethod = 'LM';

% Model orders 
na = 10;     nb = 10;     nd = 0;
orders = [na nd nb];
[TH,res,criteria,output] = pcnarx(Y,X,Q_est',[na nd nb],INDX,SelReg,THij,options);
save('SF_RandomExc_PCNARX_lsqnonlin_normal_rss.mat','TH*','criteria','res','INDX','SelReg','orders')


%% PC-NARX validation
% =========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced/';
end
cd(local_dir)
load('SF_RandomExc_PCNARX_lsqnonlin_normal_rss.mat')

% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 35;

% Input-Output data
[Y,X] = deal(zeros(NoS,K));
for k = 1:K
    out = load(['ResponsesGaussian/SFres',num2str(k+40),'.txt']);
    Y(:,k) = out(1:NoS,9);
    X(:,k) = out(1:NoS,2);
end
load('ResponsesGaussian/InputVars.mat','Q_val');

options.basis = 'hermi';
clear res
noparam = length(INDX{1})*length(SelReg);
for ii = 1:K
    [~,an(ii,:)] = PCparam(Q_val(ii,:)',2,INDX{1},reshape(TH.sim,noparam,1),options);
    Ysim(:,ii) = narxsim(X(:,ii),Y(1:10,ii),10,SelReg,an(ii,:)');
    res(:,ii) = Y(11:end,ii) - Ysim(11:end,ii);
    rss(ii) = 100*norm(res(:,ii))^2/norm(Y(11:end,ii))^2;
    mnse(ii) = mean(abs(res(:,ii))./(1+abs(Y(11:end,ii))));
end

figure(1)
subplot(211),plot(1:K,rss)
xlabel('Experiment number')
ylabel('RSS/SSS %')
subplot(212),plot(1:K,mnse)
xlabel('Experiment number')
ylabel('MNSE')




%% Parameter Surfaces -- [Figures]
% =========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced/';
end
cd(local_dir)
load('SF_RandomExc_PCNARX_lsqnonlin_normal_rss.mat')

pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

PCBASES = INDX{1};
NLREGRESSORS = SelReg;
options.basis = 'hermi';
NoI1 = 200;
NoI2 = 200;
aux = linspace(-1,1,NoI1);
Qsurf(1,:) = kron(ones(1,NoI2),aux);
Qsurf(2,:) = kron(aux,ones(1,NoI2));
noparam = length(INDX{1})*length(SelReg);
[~,an_interp] = PCparam(Qsurf,3,INDX{1},reshape(TH.sim,noparam,1),options);


figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
% subplot(2,2,1)
surf(linspace(0,2500,NoI2),linspace(0,5000,NoI1),reshape(an_interp(:,1),NoI1,NoI2))
shading interp
set(gca,'Fontsize',10)
ylabel('$k_{\mathrm{nl}}$','Interpreter','Latex','Fontsize',15)
xlabel('$F_{\max}$','Interpreter','Latex','Fontsize',15)
zlabel(['$\hat{\theta}_{y^3[t-1]}\ (\xi)$'],'Interpreter','Latex','Fontsize',15)
grid on
box on
% zlim([-0.125 0])
view([125 25])
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'SF_surf']),close;
     result = eps2xxx([write_dir,'SF_surf.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end



%% =========================================================================

%% Figures

%% =========================================================================


%% Uncertain properties -- Figure
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
cd(local_dir)
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
K = 40; 
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 

load('Responses/EstimationSet/InputVars.mat','UNCpars');
% Material 1: Vertical beams 
EX  = 200E9*UNCpars(:,1);        % N/m^2	
sigma = UNCpars(:,2);


figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(6,1,1:3)
hl1 = line(1:K,EX(:,1)/1e9,'Color','b','Marker','o');
ax1 = gca;
set(ax1,'XColor','k','YColor','k','xlim',[1,40],'ylim',[188 212],'YTick',[188:4:212],...
        'xtick',[0:5:40],'xticklabel',[],'Fontname','TimesNewRoman','Fontsize',9)    
ylabel('$E$','Fontsize',12,'Interpreter','Latex')
hleg1 = legend({'$E$'},'Interpreter','Latex','Fontsize',12);
set(hleg1,'Position',[0.105 0.94 0.4 0.05]);
ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','right',...
           'Color','none','XColor','k','YColor','k','XTick',[0:5:40],'YTick',[0.5:6.5],...
           'xlim',[1,40],'ylim',[0.5 6.5],'Fontname','TimesNewRoman','Fontsize',9);
hl2 = line(1:K,sigma(:,1),'Color','r','Marker','d','Linestyle','--','Parent',ax2);
box on
hleg2 = legend({'$\sigma_{\!\! _a}$'},'Interpreter','Latex','Fontsize',12);
set(hleg2,'Position',[0.53 0.94 0.4 0.05]);
xlabel('Experiment number','Fontsize',10,'Fontname','TimesNewRoman')
ylabel('$\sigma_{\!\! _a}$','Fontsize',12,'Interpreter','Latex','Rotation',270)
grid on,hold on
[maxs,Jmax] = max(sigma);
[mins,Jmin] = min(sigma);
% plot(Jmin,mins,'rd',Jmax,maxs,'rd','Linewidth',1.2,'MArkerFacecolor','r')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'SF_unprops']),close;
     result = eps2xxx([write_dir,'SF_unprops.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% Time histories plot (Node 54 acceleration response) - [min(sigma_a) and max(sigma_a) cases]
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
cd(local_dir)
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
N = 40; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;

% load('Responses\EstimationSet\InputVars.mat','UNCpars');
load('ResponsesGaussian\InputVars.mat','UNCpars_est');
sigma = UNCpars_est(:,2);
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
    out = load(['ResponsesGaussian\SFres',num2str(i),'.txt']);
    time = out(:,1);
    acci = out(:,2);
    acc =  out(:,9);
    
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
% axis([0 1/(2*Ts) -120 0])
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
    print(pictype,resolution,[write_dir,'SF_signals']),close;
    result = eps2xxx([write_dir,'SF_signals.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% Simulated responses -- Figures
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
cd(local_dir)
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

% Number of experiments
N = 20; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;

load('Responses\40\InputVars.mat','UNCpars');
sigma = UNCpars(:,2);
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
    out = load(['Responses\40\SFres',num2str(i),'.txt']);
    time = out(:,1);
    resforce(:,p) = out(:,3);
    dspl(:,p) =  out(:,8);
end

figure(1)
clr = colormap(lines(2));
hold on
figindx = [1 1 2 2];
colindx = [1 2 1 2];
linindx = [1 1 1 1];
subplot(2,2,1),plot(dspl(1:250,1)/1000,resforce(1:250,1)/1000000,'o-',...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
hold on
% str = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(Jmin,1))),'\ (N/m^3), \sigma_f = ',num2str(round(Qtrans(ExpLow ,2))),'\ (N)$']};
% legend(str,'Fontsize',9,'Location','Northoutside','Interpreter','Latex')
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Force (MN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Displacement (mm)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,2,5),plot(0.005*(0:999),dspl(1:1000,1)/1000,...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
set(gca,'xticklabel',[])
ylabel('Displacement (mm)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,2,7),plot(0.005*(0:999),resforce(1:1000,1)/1000000,...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
grid on;
ylim([-2 2])
ylabel('Force (MN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(2,2,2),plot(dspl(1:250,2),resforce(1:250,2)/1000000,'-o',...
    'Linewidth',linindx(1),'Color',clr(colindx(2),:))
hold on
grid on
% str = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(ExpHigh,1))),'\ (N/m^3), \sigma_f = ',num2str(round(Qtrans(ExpHigh,2))),'\ (N)$']};
% legend(str,'Fontsize',9,'Location','Northoutside','Interpreter','Latex')
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Force (MN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(4,2,6),plot(0.005*(0:749),dspl(1:750,2),...
    'Linewidth',linindx(2),'Color',clr(colindx(2),:))
set(gca,'xticklabel',[])
ylim([-2 2])
grid on
ylabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(4,2,8),plot(0.005*(0:749),resforce(1:750,2)/1000000,...
    'Linewidth',linindx(2),'Color',clr(colindx(2),:))
% axis tight
ylim([-10 10])
grid on
ylabel('Force (MN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'SF_randexci']),close;
    if ispc
        result = eps2xxx([write_dir,'SF_randexci.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'SF_randexci.eps'],{'pdf'});
    end
end


%% NARX regressors selection -- Figure
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
load('SDOF_Cubic_RandomExc_GA_rss.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

NoReg = length(GAreg);
close all

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(211),plot(1:NoReg,SEstructure.values(1:NoReg),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 NoReg],SEstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
set(gca,'yscale','linear','Fontsize',7,'xtick',[1:NoReg],'xticklabel',[])
for ii = 1:NoReg
    if strcmp(SEstructure.removed{ii},'ones(length(t),size(Y,2))')
        txtstr = '$\mbox{const.term}$';
    else
        if str2num(SEstructure.removed{ii}(12)) > 1
            txtstr = ['$\ ',lower(SEstructure.removed{ii}(2)),'[t -',SEstructure.removed{ii}(6),']',SEstructure.removed{ii}(11:12),'$'];
        else
            txtstr = ['$\ ',lower(SEstructure.removed{ii}(2)),'[t -',SEstructure.removed{ii}(6),']$'];
        end
    end
    text(ii,-10,txtstr,'Rotation',25,'Fontsize',7,'Horizontalalignment','center','Interpreter','Latex','Fontsize',10)
end
axis([0.9 11.1 0 100])
grid on
text(NoReg/2+0.5,-22.5,'Regressors dropped','Fontangle','normal','Fontsize',9,'Horizontalalignment','center')
ylabel('RSS/SSS (%)','Fontangle','normal','Fontsize',9)
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_MSS_StageA']),close;    
    if ispc
        result = eps2xxx([write_dir,'cubic_MSS_StageA.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_MSS_StageA.eps'],{'pdf'});
    end
end


%% Local NARX model estimation -- Figures
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end

K = 20;
figure(1)
clr = colormap(lines(3));
subplot(211),hold on 
xlabel('# Experiments')
ylabel('RSS/SSS %')
subplot(212),hold on 
ylabel('MNSE')

for j  = 1:3
    switch j
        case 1
            load('SDOF_Cubic_RandomExc_LocalNARX_lsqnonlin.mat')
        case 2
            load('SDOF_Cubic_RandomExc_LocalNARX_fminunc.mat')
        case 3
            load('SDOF_Cubic_RandomExc_LocalNARX_fminsearch.mat')
    end
    
    for k=1:K
        rssSIM(j,k) = SIMcriteria{k}.rss_sss;
        mnseSIM(j,k) = SIMcriteria{k}.mnse;
        exitflag(j,k) = SIMcriteria{k}.SIMexitflag;
    end
end
subplot(211),plot(1:K,rssSIM)
subplot(212),plot(1:K,mnseSIM)



%% Parameters expansion -- Figure
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
load('SDOF_Cubic_RandomExc_LHS.mat')
load('SDOF_Cubic_RandomExc_LocalNARX_lsqnonlin.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

clear basis options
options.basis = 'legen';
options.criterion = 'rss';
for P = 0:4
    % Basis index
    options.Nb = 0;
    [BASISstructure,INDnew] = PCEreduction(Q',P,thetaSIM,options);
    R2o(P+1) = BASISstructure.initvalue;
end

options.Nb = 2; Nb = 2;
% Basis index
P = 1;
[BASISstructure,INDnew] = PCEreduction(Q',P,thetaSIM,options);

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(221),plot(0:4,R2o,'-o','Markersize',6,'Linewidth',1.2)
grid on
ylabel('RSS/SSS (%)','Fontsize',11)
xlabel(sprintf('%s\n%s','Total PC basis degree','(complete subspace)'),'Fontsize',10)
set(gca,'yscale','linear','Fontsize',7,'xtick',[0:4])
subplot(222),plot(1:Nb,BASISstructure.values(1:Nb),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 Nb],BASISstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
ylabel('RSS/SSS (%)','Fontsize',11)
set(gca,'yscale','linear','Fontsize',7,'xtick',[1:Nb],'xticklabel',[])
for ii = 1:Nb 
    text(ii,-2e-3,['$[',num2str(BASISstructure.removed(ii,1)),',',num2str(BASISstructure.removed(ii,2)),']$'],'Rotation',0,'Fontsize',7,'Horizontalalignment','center','Interpreter','Latex','Fontsize',11)
end
% axis([0.99 2.01 0.9997 1])
grid on
text(1.5,-5e-3,sprintf('%s\n%s','PC bases dropped','(multivariable indeces)'),'Fontsize',10,'HorizontalAlignment','center')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'cubic_MSS_StageB']),close;
     if ispc
        result = eps2xxx([write_dir,'cubic_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
        result = eps2xxx([write_dir,'cubic_MSS_StageB.eps'],{'pdf'});
     end
end

%% Local NARX vs Global PC-NARX R2 criterion -- Figures
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

K = 10;
load('SDOF_Cubic_RandomExc_LocalNARX_lsqnonlin.mat')
for k=1:K
    rss(1,k) = SIMcriteria{k}.rss_sss;
    mnse(1,k) = SIMcriteria{k}.mnse;
    R2(1,k) = 1- SIMcriteria{k}.rss_sss/100;
    exitflag(1,k) = SIMcriteria{k}.SIMexitflag;
end

load('SDOF_Cubic_RandomExc_LHS.mat')
load('SDOF_Cubic_RandomExc_PCNARX_lsqnonlin.mat')
clear res
options.basis = 'legen';
for ii = 1:length(Q)
    [~,an(ii,:)] = PCparam(Q(ii,:)',1,INDX{1},TH.sim(:),options);
    X = force(1:1000,ii);
    Y = dspl(1:1000,ii);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    res(:,ii) = Y(3:end) - Ysim(3:end);
    rss(2,ii) = 100*norm(res(:,ii))^2/norm(Y(3:end))^2;
    R2(2,ii) = 1 - rss(2,ii)/100;
    mnse(2,ii) = mean(abs(res(:,ii))./(1+abs(Y(3:end))));
end


figure(1)
subplot(4,1,1:2),plot(1:K,R2,'-o')
legend({'Local NARX models','Global PC-NARX model'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('R^2','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
axis([0.9 10.1 0.998 1])
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_local_vs_global']),close;
    if ispc
        result = eps2xxx([write_dir,'cubic_local_vs_global.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_local_vs_global.eps'],{'pdf'});
    end
end



%% Global PC-NARX R2 criterion [estimation vs validation set] -- Figures
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

K = 10;
load('SDOF_Cubic_RandomExc_LHS.mat')
load('SDOF_Cubic_RandomExc_PCNARX_lsqnonlin.mat')
clear res
options.basis = 'legen';
for ii = 1:K
    [~,an(ii,:)] = PCparam(Q(ii,:)',1,INDX{1},TH.sim(:),options);
    X = force(1:1000,ii);
    Y = dspl(1:1000,ii);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    res(:,ii) = Y(3:end) - Ysim(3:end);
    rss(ii) = 100*norm(res(:,ii))^2/norm(Y(3:end))^2;
    R2(ii) = 1 - rss(ii)/100;
    mnse(ii) = mean(abs(res(:,ii))./(1+abs(Y(3:end))));
end

clear res X Y Q
load('SDOF_Cubic_RandomExc_validation.mat')
for ii = 1:K
    [~,an(ii,:)] = PCparam(Q(ii,:)',1,INDX{1},TH.sim(:),options);
    X = force(1:1000,ii);
    Y = dspl(1:1000,ii);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    res(:,ii) = Y(3:end) - Ysim(3:end);
    rss(ii+K) = 100*norm(res(:,ii))^2/norm(Y(3:end))^2;
    R2(ii+K) = 1 - rss(ii+K)/100;
    mnse(ii+K) = mean(abs(res(:,ii))./(1+abs(Y(3:end))));
end


figure(1)
subplot(4,8,[1:4 9:12]),plot(1:K,R2(1:K),'-bo')
legend({'Estimation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
axis([1 10 0.999 1])
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
ylabel('R^2','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,8,[5:8 13:16]),plot(1:K,R2(K+1:2*K),'-rd')
legend({'Validation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
set(gca,'Fontsize',9,'Fontname','TimesNewRoman','yticklabel',[])
axis([1 10 0.999 1])
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_est_vs_val']),close;
    if ispc
        result = eps2xxx([write_dir,'cubic_est_vs_val.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_est_vs_val.eps'],{'pdf'});
    end
end



%% Parameter Surfaces -- [Figures]
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

load('SDOF_Cubic_RandomExc_PCNARX_lsqnonlin.mat')
PCBASES = INDX{1};
NLREGRESSORS = SelReg;
options.basis = 'legen';
NoI1 = 20;
NoI2 = 40;
aux = linspace(-1,1,NoI1);
Qsurf(1,:) = kron(ones(1,NoI2),aux);
Qsurf(2,:) = kron(aux,ones(1,NoI2));
[~,an_interp] = PCparam(Qsurf,1,INDX{1},reshape(TH.sim,8,1),options);


figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
% subplot(2,2,1)
surf(linspace(0,5000,NoI2),linspace(0,5000,NoI1),reshape(an_interp(:,4),NoI1,NoI2))
shading interp
set(gca,'Fontsize',10)
ylabel('$k_{\mathrm{nl}}$','Interpreter','Latex','Fontsize',15)
xlabel('$F_{\max}$','Interpreter','Latex','Fontsize',15)
zlabel(['$\hat{\theta}_{y[t-\ \ 1]^3} \! \! (\xi)$'],'Interpreter','Latex','Fontsize',15)
grid on
box on
zlim([-0.125 0])
view([125 25])
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'cubic_surf']),close;
     if ispc
        result = eps2xxx([write_dir,'cubic_surf.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
        result = eps2xxx([write_dir,'cubic_surf.eps'],{'pdf'});
     end
end


%% ========================================================================

%% ========================================================================

%% ========================================================================


%% Parameters expansion
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels';
end
cd(local_dir)
load('ResponsesGaussian/InputVars.mat','Q');
load('SF_RandomExc_LocalNARX.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

clear basis options

options.basis = 'hermi';
options.Nb = 10;
options.criterion = 'bic';

% Basis index
P = 3;
[BASISstructure,INDnew] = PCEreduction(Q',P,thetaSIM,options);

%% PC-NARX model structure selection (FOLS)
% =========================================================================
clear,clc,close all
local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
cd(local_dir)

% Number of experiments
N = 50; 
% Number of samples
NoS = 500;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;

% Input-Output data
[Y,X] = deal(zeros(NoS,N));
for i = 1:N
    out = load(['Responses50/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:NoS,13);
    X(:,i) = out(1:NoS,2);
end
load('Responses50/InputVars.mat','Q');

% Nonlinear polynomial regressors
P = 1;      q = 1;      maxterms = 1;
% Model orders
na = 10;     nb = 10;     nd = 0;
% Nonlinear regressors (complete vector)
regressors = polyregNARX([na nd nb],P,q,maxterms,'abs');
% FOLS structure selection
options.criterion = 'rss';
PCorder = 3;
indx{1} = combinations(M,PCorder,1);
[Nr,RegSTR,regindx,basisindx,NLres,NLcriteria] = folsPCNARX(Y,X,Q,[na nd nb],indx,[],[],[],regressors,50,[],options);

%% Run candidate models for 1-50 most important regressors
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
plot(1:20,rssPE(1:20),'-bo',1:20,rssSIM(1:20),'-rd')
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
% 9     Velocity response of 1st floor (node 18)
% 10    Velocity response of 2nd floor (node 27)
% 11    Velocity response of 3rd floor (node 36)
% 12    Velocity response of 4th floor (node 45)
% 13    Velocity response of 5th floor (node 54)
% 14     Acceleration response of 1st floor (node 18)
% 15    Acceleration response of 2nd floor (node 27)
% 16    Acceleration response of 3rd floor (node 36)
% 17    Acceleration response of 4th floor (node 45)
% 18    Acceleration response of 5th floor (node 54)


%    SET   TIME/FREQ    LOAD STEP   SUBSTEP  CUMULATIVE
%      1  3.0242             1         1         1
%      2  3.0242             1         2         2
%      4  8.0454             1         4         4
%      5  8.0454             1         5         5
%      7  11.100             1         7         7
%      8  11.100             1         8         8
%     10  17.547             1        10        10
%     11  17.547             1        11        11
%     16  30.584             1        16        16
%     17  30.584             1        17        17