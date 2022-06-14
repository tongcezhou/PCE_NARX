%% Multistorey II -- TESTS (increasing acc/force)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\MultistoreySF2D_II\')

NoS = 1000;
NoExp = 2;
Ts = 0.02;
time = Ts*[1:1000];
% Generation of the input force
ACCinput = [zeros(NoS,1) ; (linspace(0,1,1000)').*(30*sin(time*pi/2.5)')];                           
Finput = [linspace(1e2,3e4,NoS)' ; zeros(NoS,1)];       
TestInput = [ACCinput ; Finput];
save('TestInput.txt','TestInput','-ascii')

% Run Ansys
delete('file.*')


%% Multistorey II - TESTS -- Figures
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')

% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 2;

% Input-Output data
[Y,X,Rc] = deal(zeros(N,K));
for i = 1:2
    AccIn = load(['MultiII_Test_',num2str(i),'.txt']);
    X(:,i) = AccIn(1:N,2);
    Rc(:,i) = AccIn(1:N,7);
    AccOut = load(['MultiII_Test_acc_',num2str(i),'.txt']);
    Y(:,i) = AccOut(1:N,16); 
    q = 0;
    for j = [1 5 9 13 15]
        q = q+1;
        Stress = load(['MultiII_Test_Sx_',num2str(j),'_',num2str(i),'.txt']);
        Sx(:,q) = Stress(1:N,3);
    end
    figure
    plot(Sx)
end



%% Multistorey -- Transient nonlinear analysis (various acc. levels - normal pdf)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\MultistoreySF2D_II\')

% Number of experiments (estimation set)
K1 = 20; 
% Number of experiments (validation set)
K2 = 20; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, std_a)
M = 2; 

% Latin Hypercube Sampling (estimation set)
rng(100)
Q_est = lhsnorm(zeros(M,1),eye(M),K1);
Qtrans_est = 1 + 0.02*Q_est(:,1);
sigma_est = exp(2.5 + 0.5*Q_est(:,2));
% Latin Hypercube Sampling (validation set)
rng(200)
Q_val = lhsnorm(zeros(M,1),eye(M),K2);
Qtrans_val = 1 + 0.02*Q_val(:,1);
sigma_val = exp(2.5 + 0.5*Q_val(:,2));

% Beam properties file
BeamProps = [Qtrans_est(:) ;Qtrans_val(:)];
save('BeamProps20.txt','BeamProps','-ascii')

% Generation of the random signal excitation
for k=1:K1+K2 
    rng(k)
    if k<=K1
        Exc((k-1)*NoS+1:k*NoS,1) = sigma_est(k)*randn(NoS,1);
    else
        Exc((k-1)*NoS+1:k*NoS,1) = sigma_val(k-K1)*randn(NoS,1);
    end
end
save('AccRandom20.txt','Exc','-ascii')

UNCpars_est = [Qtrans_est sigma_est];
UNCpars_val = [Qtrans_val sigma_val];
Exc = reshape(Exc,NoS,K1+K2);
save('InputVars20.mat','Q_est','Q_val','Exc','UNCpars*')

% Run Ansys
delete('file.*')


%% Multistorey II -- Uncertain parameters (norm pdf) -- Figures
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('InputVars.mat')
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
K = 80; 
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 

load('InputVars.mat','UNCpars*');
% Material 1: Vertical beams 
EX  = 200E9*UNCpars_est(:,1);        % N/m^2	
sigma = UNCpars_est(:,2);


figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(6,1,1:3)
hl1 = line(1:K,EX(:,1)/1e9,'Color','b','Marker','o');
ax1 = gca;
set(ax1,'XColor','k','YColor','k','xlim',[1,80],'ylim',[188 212],'YTick',188:4:212,...
        'xtick',0:5:80,'xticklabel',[],'Fontname','TimesNewRoman','Fontsize',9)    
ylabel('$E$','Fontsize',12,'Interpreter','Latex')
hleg1 = legend({'$E$'},'Interpreter','Latex','Fontsize',12);
set(hleg1,'Position',[0.105 0.94 0.4 0.05]);
ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','right',...
           'Color','none','XColor','k','YColor','k','XTick',[0:5:80],'YTick',[0:10:40],...
           'xlim',[1,80],'ylim',[0 40],'Fontname','TimesNewRoman','Fontsize',9);
hl2 = line(1:K,sigma(:,1),'Color','r','Marker','d','Linestyle','--','Parent',ax2);
box on
hleg2 = legend({'$\sigma_{\!\! _a}$'},'Interpreter','Latex','Fontsize',12);
set(hleg2,'Position',[0.53 0.94 0.4 0.05]);
xlabel('Experiment number','Fontsize',10,'Fontname','TimesNewRoman')
ylabel('$\sigma_{\!\! _a}$','Fontsize',12,'Interpreter','Latex','Rotation',270)
grid on,hold on
[maxs,Jmax] = max(sigma);
[mins,Jmin] = min(sigma);
plot(Jmin,mins,'rd',Jmax,maxs,'rd','Linewidth',1.2,'MArkerFacecolor','r')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'SF_unprops']),close;
     result = eps2xxx([write_dir,'SF_unprops.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end




%% Time histories plot - [min(sigma_a) and max(sigma_a) cases]
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
N = 80; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density)
M = 2; 
% Sampling period
Ts =  0.02;

load('InputVars.mat','UNCpars_est');
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
    out = load(['MultiII_',num2str(i),'.txt']);
    time = out(:,1);
    acci = out(:,2);
    out = load(['MultiII_vel_',num2str(i),'.txt']);
    acc =  out(:,1);
    
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
cd('D:\Data\FEM\MultistoreySF2D_II\')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

% Number of experiments
N = 80; 
% Number of samples
NoS = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;

load('InputVars.mat');
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
    out = load(['MultiII_',num2str(i),'.txt']);
    time = out(:,1);
    resforce(:,p) = out(:,7);
    out = load(['MultiII_dspl_',num2str(i),'.txt']);
    dspl(:,p) =  out(:,2)*1000;
end

figure(1)
clr = colormap(lines(2));
hold on
figindx = [1 1 2 2];
colindx = [1 2 1 2];
linindx = [1 1 1 1];
subplot(2,2,1),plot(dspl(1:250,1),resforce(1:250,1)/1000000,'o-',...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
hold on
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Force (MN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Displacement (mm)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,2,5),plot(0.025*(0:999),dspl(1:1000,1),...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
set(gca,'xticklabel',[])
ylabel('Displacement (mm)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,2,7),plot(0.025*(0:999),resforce(1:1000,1)/1000000,...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
grid on;
ylim([-0.2 0.2])
ylabel('Force (MN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(2,2,2),plot(dspl(1:250,2),resforce(1:250,2)/1000000,'-o',...
    'Linewidth',linindx(1),'Color',clr(colindx(2),:))
hold on
grid on
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Force (MN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Displacement (mm)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(4,2,6),plot(0.025*(0:749),dspl(1:750,2),...
    'Linewidth',linindx(2),'Color',clr(colindx(2),:))
set(gca,'xticklabel',[])
ylim([-75 50])
grid on
ylabel('Displacement (mm)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(4,2,8),plot(0.025*(0:749),resforce(1:750,2)/1000000,...
    'Linewidth',linindx(2),'Color',clr(colindx(2),:))
% axis tight
ylim([-2 2])
grid on
ylabel('Force (MN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'SF_randexci']),close;
    result = eps2xxx([write_dir,'SF_randexci.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end



%% NARX regressors selection (GA - rss)
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')

% Number of experiments
K = 80; 
% Number of samples
N = 1000;
% Number of variables (Ex, Density)
M = 2; 
% Sampling period
Ts =  0.02;

% Simulation experiment
SimExp = 26;
AccIn = load(['MultiII_',num2str(SimExp),'.txt']);
X = AccIn(201:N,2);
VelOut = load(['MultiII_vel_',num2str(SimExp),'.txt']);
Y = VelOut(201:N,2);

% for i = 1:16
%     Y = VelOut(201:N,i);
%     [P(:,i),F] = pwelch(Y,256,200,512,40);
% end
% plot(20*log10(abs(P(:,:))))


% Model orders 
na = 14;     nb = 14;     nd = 0;
% Polynomial regressors
P = 3;      maxterms = 1;
% regressors = polyregNARX([na nd nb],P,1,maxterms);
regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))',...
     'abs(Y(t-2,:))','abs(Y(t-3,:))','abs(Y(t-4,:))','abs(Y(t-5,:))',...
     'abs(Y(t-6,:))','abs(Y(t-7,:))'});
 
% Optimization options
options.warnings = 'off';
options.GA.PopulationSize = 200;
options.maxsize = 1e8;

rng(10);
options.focus = 'simulation';
options.criterion = 'rss';
options.parfor = 'n';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);
options.focus = 'prediction';
options.Nr = length(GAreg);
[SEstructure,IndxRem] = SEreduction(Y,X,[na nd nb],GAreg,options);
save('MultiII_RandomExc_GA_rss.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')



%% NARX regressors selection (GA - rss - p = 1)
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')

% Number of experiments
K = 80; 
% Number of samples
N = 1000;
% Number of variables (Ex, Density)
M = 2; 
% Sampling period
Ts =  0.02;

% Simulation experiment
SimExp = 26;
AccIn = load(['MultiII_',num2str(SimExp),'.txt']);
X = AccIn(201:N,2);
VelOut = load(['MultiII_vel_',num2str(SimExp),'.txt']);
Y = VelOut(201:N,2);

% for i = 1:16
%     Y = VelOut(201:N,i);
%     [P(:,i),F] = pwelch(Y,256,200,512,40);
% end
% plot(20*log10(abs(P(:,:))))


% Model orders 
na = 14;     nb = 14;     nd = 0;
% Polynomial regressors
P = 1;      maxterms = 1;
% regressors = polyregNARX([na nd nb],P,1,maxterms);
regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))',...
     'abs(Y(t-2,:))','abs(Y(t-3,:))','abs(Y(t-4,:))','abs(Y(t-5,:))',...
     'abs(Y(t-6,:))','abs(Y(t-7,:))','abs(Y(t-8,:))','abs(Y(t-9,:))','abs(Y(t-10,:))','abs(Y(t-11,:))',...
     'abs(Y(t-12,:))','abs(Y(t-13,:))','abs(Y(t-14,:))'});

% Optimization options
options.warnings = 'off';
options.GA.PopulationSize = 200;
options.maxsize = 1e8;

rng(100);
options.focus = 'simulation';
options.criterion = 'rss';
options.parfor = 'n';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);
options.focus = 'prediction';
options.Nr = length(GAreg);
[SEstructure,IndxRem] = SEreduction(Y,X,[na nd nb],GAreg,options);
save('MultiII_RandomExc_GA_rss_p1.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')



%% NARX regressors selection -- Figure
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('MultiII_RandomExc_GA_rss.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

NoReg = length(GAreg);
close all

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(211),plot(1:NoReg,SEstructure.values(1:NoReg),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 NoReg],SEstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
% set(gca,'yscale','linear','Fontsize',7,'xtick',[1:NoReg],'xticklabel',[])
% for ii = 1:NoReg
%     if strcmp(SEstructure.removed{ii},'ones(length(t),size(Y,2))')
%         txtstr = '$\mbox{const.term}$';
%     else
%         if str2num(SEstructure.removed{ii}(12)) > 1
%             txtstr = ['$\ ',lower(SEstructure.removed{ii}(2)),'[t -',SEstructure.removed{ii}(6),']',SEstructure.removed{ii}(11:12),'$'];
%         else
%             txtstr = ['$\ ',lower(SEstructure.removed{ii}(2)),'[t -',SEstructure.removed{ii}(6),']$'];
%         end
%     end
%     text(ii,-10,txtstr,'Rotation',25,'Fontsize',7,'Horizontalalignment','center','Interpreter','Latex','Fontsize',10)
% end
% axis([0.9 11.1 0 100])
grid on
text(NoReg/2+0.5,-22.5,'Regressors dropped','Fontangle','normal','Fontsize',9,'Horizontalalignment','center')
ylabel('RSS/SSS (%)','Fontangle','normal','Fontsize',9)

subplot(212),plot(2:NoReg,100*abs(diff(SEstructure.values(1:NoReg))./SEstructure.values(1:NoReg-1)),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 NoReg],[1 1],'--r','Markersize',6,'Linewidth',1.2)
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'multi_MSS_StageA']),close;    
    if ispc
        result = eps2xxx([write_dir,'multi_MSS_StageA.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'multi_MSS_StageA.eps'],{'pdf'});
    end
end



%% Local NARX model estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('MultiII_RandomExc_GA_rss.mat')
% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 80;

% Input-Output data
[Y,X] = deal(zeros(N-200,K));
for i = 1:K
    AccIn = load(['MultiII_',num2str(i),'.txt']);
    X(:,i) = AccIn(201:N,2);
    AccOut = load(['MultiII_vel_',num2str(i),'.txt']);
    Y(:,i) = AccOut(201:N,2); 
end
load('InputVars.mat','Q_est');

% Regressors
SelReg = SEstructure.removed(274:end)';

% Model orders 
na = 14;     nb = 14;     nd = 0;
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
options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','iter',...
            'TolFun',1e-3,'TolX',1e-6,'MaxFunEval',10000,'MaxIter',10000);
% matlabpool('open',4)
% parfor k = 1:K
for k = 1:K
    disp(k)
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
% matlabpool('close')

save('MultiII_RandomExc_LocalNARX_lsqnonlin.mat','*criteria','NLres*','theta*','SelReg','orders')

%% Local NARX model estimation (lsqnonlin) -- terms 292,...
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('MultiII_RandomExc_GA_rss.mat')
% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 80;

% Input-Output data
[Y,X] = deal(zeros(N-200,K));
for i = 1:K
    AccIn = load(['MultiII_',num2str(i),'.txt']);
    X(:,i) = AccIn(201:N,2);
    AccOut = load(['MultiII_vel_',num2str(i),'.txt']);
    Y(:,i) = AccOut(201:N,2); 
end
load('InputVars.mat','Q_est');

% Regressors
SelReg = SEstructure.removed(292:end)';

% Model orders 
na = 14;     nb = 14;     nd = 0;
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
options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','final',...
            'TolFun',1e-3,'TolX',1e-6,'MaxFunEval',10000,'MaxIter',10000);
matlabpool('open',4)
parfor k = 1:K
% for k = 1:K
    disp(k)
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
matlabpool('close')

save('MultiII_RandomExc_LocalNARX_lsqnonlin292.mat','*criteria','NLres*','theta*','SelReg','orders')



%% Local NARX model estimation (lsqnonlin) -- 50 terms
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('MultiII_RandomExc_GA_rss.mat')
% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 80;

% Input-Output data
[Y,X] = deal(zeros(N-200,K));
for i = 1:K
    AccIn = load(['MultiII_',num2str(i),'.txt']);
    X(:,i) = AccIn(201:N,2);
    AccOut = load(['MultiII_vel_',num2str(i),'.txt']);
    Y(:,i) = AccOut(201:N,2); 
end
load('InputVars.mat','Q_est');

% Regressors
SelReg = SEstructure.removed(292:end)';

% Model orders 
na = 14;     nb = 14;     nd = 0;
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
options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','final',...
            'TolFun',1e-3,'TolX',1e-6,'MaxFunEval',10000,'MaxIter',10000);
matlabpool('open',4)
parfor k = 1:K
% for k = 1:K
    disp(k)
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
matlabpool('close')

save('MultiII_RandomExc_LocalNARX_lsqnonlin292.mat','*criteria','NLres*','theta*','SelReg','orders')


%% Local NARX model estimation (lsqnonlin) -- P=1
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('MultiII_RandomExc_GA_rss_p1.mat')
% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 80;

% Input-Output data
[Y,X] = deal(zeros(N-200,K));
for i = 1:K
    AccIn = load(['MultiII_',num2str(i),'.txt']);
    X(:,i) = AccIn(201:N,2);
    AccOut = load(['MultiII_vel_',num2str(i),'.txt']);
    Y(:,i) = AccOut(201:N,2); 
end
load('InputVars.mat','Q_est');

% Regressors
SelReg = SEstructure.removed(179:end)';

% Model orders 
na = 14;     nb = 14;     nd = 0;
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
options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','final',...
            'TolFun',1e-3,'TolX',1e-6,'MaxFunEval',10000,'MaxIter',10000);
matlabpool('open',4)
parfor k = 1:K
% for k = 1:K
    disp(k)
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
matlabpool('close')

save('MultiII_RandomExc_LocalNARX_lsqnonlin_p1.mat','*criteria','NLres*','theta*','SelReg','orders')

%% Local NARX model estimation (lsqnonlin) --more terms
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('MultiII_RandomExc_GA_rss.mat')
% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 80;

% Input-Output data
[Y,X] = deal(zeros(N-200,K));
for i = 1:K
    AccIn = load(['MultiII_',num2str(i),'.txt']);
    X(:,i) = AccIn(201:N,2);
    AccOut = load(['MultiII_vel_',num2str(i),'.txt']);
    Y(:,i) = AccOut(201:N,2); 
end
load('InputVars.mat','Q_est');

% Regressors
SelReg = SEstructure.removed(72:end)';

% Model orders 
na = 14;     nb = 14;     nd = 0;
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
options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','iter',...
            'TolFun',1e-3,'TolX',1e-6,'MaxFunEval',10000,'MaxIter',10000);
% matlabpool('open',4)
% parfor k = 1:K
for k = 1:5
    disp(k)
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
% matlabpool('close')



%% Parameters expansion -- GA -- Figure
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
load('MultiII_RandomExc_GA_rss.mat')
load('MultiII_RandomExc_LocalNARX_lsqnonlin292.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

clear basis options
options.basis = 'hermi';
options.ortho = 'y';
options.criterion = 'ar2';


% Number of samples
N = 800;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 80;
load('InputVars.mat','Q_est');

Pmax =5;
options.GA.PopulationSize = 100;
[A,theta,BASISstructure,criteria,PHI,GAoutput] = gaPC(thetaSIM',Q_est',Pmax,1,options);
NoB = length(BASISstructure.values);

close all
figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(211),plot(1:NoB,BASISstructure.values(1:NoB),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 NoB],BASISstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
set(gca,'yscale','linear','Fontsize',7,'xtick',[1:NoB],'xticklabel',[])
for ii = 1:NoB
    txtstr = ['$[',num2str(BASISstructure.removed(ii,1)),',',num2str(BASISstructure.removed(ii,2)),']$'];
    text(ii,-0.2,txtstr,'Rotation',25,'Fontsize',7,'Horizontalalignment','center','Interpreter','Latex','Fontsize',10)
end
axis([1 NoB -0.1 0.6])
grid on
text(NoB/2+0.5,-0.4,'Basis functions dropped','Fontangle','normal','Fontsize',9,'Horizontalalignment','center')
% ylabel('mean(R^2_{adj})','Fontangle','normal','Fontsize',9)
ylabel('$\overline{R^2_{\mbox{adj}}}$','Fontangle','normal','Fontsize',11,'Interpreter','Latex')
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
cd('D:\Data\FEM\MultistoreySF2D_II\')
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
load('MultiII_RandomExc_GA_rss.mat')
load('MultiII_RandomExc_LocalNARX_lsqnonlin292.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

% Number of samples
N = 800;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 80;
% Input-Output data
% [Y,X] = deal(zeros(N-200,K));
% for i = 1:K
%     AccIn = load(['MultiII_',num2str(i),'.txt']);
%     X(:,i) = AccIn(201:N,2);
%     AccOut = load(['MultiII_vel_',num2str(i),'.txt']);
%     Y(:,i) = AccOut(201:N,2); 
% end
load('InputVars.mat','Q_est');

clear basis options
options.basis = 'hermi';
options.ortho = 'y';
options.criterion = 'mne';
Pmax = 10;
for P = 0:Pmax
    % Basis index
    options.Nb = 0;
    [BASISstructure,INDnew] = PCEreduction(Q_est',P,thetaSIM,options);
    R2o(P+1)= BASISstructure.initvalue;
end

% Basis index
P = 7;
options.Nb = 35; Nb = 35;
[BASISstructure,INDnew] = PCEreduction(Q_est',P,thetaSIM,options);


figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(221),plot(0:Pmax,R2o,'-o','Markersize',6,'Linewidth',1.2)
grid on
ylabel('Mean normalized error','Fontsize',11)
xlabel(sprintf('%s\n%s','Total PC basis degree','(complete subspace)'),'Fontsize',10)
set(gca,'yscale','linear','Fontsize',9,'xtick',[0:Pmax])
% axis([0 4 0.044 0.062])
subplot(222),plot(1:Nb,BASISstructure.values(1:Nb),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 Nb],BASISstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
ylabel('Mean normalized error','Fontsize',11)
set(gca,'yscale','linear','Fontsize',9,'xtick',[1:Nb],'xticklabel',[])
for ii = 1:Nb 
    text(ii,0.0435,['$[',num2str(BASISstructure.removed(ii,1)),',',num2str(BASISstructure.removed(ii,2)),']$'],'Rotation',45,...
        'Horizontalalignment','center','Interpreter','Latex','Fontsize',9)
end
% axis([0.99 Nb+0.01 0.045 0.061])
grid on
text(5,0.0405,sprintf('%s\n%s','PC bases dropped','(multivariable indeces)'),'Fontsize',10,'HorizontalAlignment','center')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'SF_MSS_StageB']),close;
     if ispc
        result = eps2xxx([write_dir,'SF_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
        result = eps2xxx([write_dir,'SF_MSS_StageB.eps'],{'pdf'});
     end
end


%% PC-NARX estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('MultiII_RandomExc_GA_rss.mat')

% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;
% Number of experiments
K = 80;

% Input-Output data
[Y,X] = deal(zeros(N-200,K));
load('InputVars.mat','Q_est');
load('MultiII_RandomExc_LocalNARX_lsqnonlin292.mat','*criteria','NLres*','theta*','SelReg','orders')

for i = 1:K
    AccIn = load(['MultiII_',num2str(i),'.txt']);
    X(:,i) = AccIn(201:N,2);
    VelOut = load(['MultiII_vel_',num2str(i),'.txt']);
    Y(:,i) = VelOut(201:N,2); 
    rssPE(i) = PEcriteria{i}.rss_sss;
    rssSIM(i) = SIMcriteria{i}.rss_sss;
    mnseSIM(i) = SIMcriteria{i}.mnse;
end

clear basis options
options.basis = 'hermi';
options.ortho = 'y';

% Basis indx
% INDX{1} = [ 0 0;
%      1     0
%      0     1
%      2     0
%      0     2
%      2     1
%      0     3];
% INDX{1} = [ 0 0;
%      1     0;
%      0     1;
%      2     0;
%      0     2;
%      1     2;
%      0     3;
%      2     2;
%      0     4;
%      0     5;
%      1     4;
%      0     6;
%      2     4;
%      0     7;
%      4     3;
%      6     1;
%      1     6];
     
     
INDX{1} = [ 0 0;     
    1     0;
    0     1;
    0     2;
    0     3;
    1     2;
    0     4;
    2     2;
    5     0;
    0     5;
    4     1;
    1     4;
    3     2]; 


 % Total number of basis functions
B = size(INDX{1},1);
% Basis constraction: size(basis_i) = p x K
basis = cell(M,1);
% Total PC basis degree 
P = 5;
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
% Regressors
SelReg = SEstructure.removed(292:end)';

% Model orders 
na = 14;     nb = 14;     nd = 0;
orders = [na nd nb];
[TH,res,criteria,output] = pcnarx(Y,X,Q_est',[na nd nb],INDX,SelReg,[],options);
% [TH,res,criteria,output] = pcnarx(Y,X,Q_est',[na nd nb],INDX,SelReg,THij,options);

save('SF_RandomExc_PCNARX_lsqnonlin_normal_largePC.mat','TH*','criteria','res','INDX','SelReg','orders')

 

for k =1:20
    figure(k)
    plot(Y(15:end,k))
    hold on 
    plot(Y(15:end,k)-NLresSIM{k},'--r')
end




%% Global PC-NARX R2 criterion [estimation vs validation set] -- Figures
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('MultiII_RandomExc_GA_rss.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

N = 1000;
K = 80;
load('InputVars.mat');
load('SF_RandomExc_PCNARX_lsqnonlin_normal_largePC.mat')
clear res
options.basis = 'hermi';
options.ortho = 'y';

for ii = 1:K
    [~,an(ii,:)] = PCparam(Q_est(ii,:)',5,INDX{1},TH.sim(:),options);
    AccIn = load(['MultiII_',num2str(ii),'.txt']);
    X(:,ii) = AccIn(201:N,2);
    VelOut = load(['MultiII_vel_',num2str(ii),'.txt']);
    Y(:,ii) = VelOut(201:N,2); 
    Ysim = narxsim(X(:,ii),Y(1:14,ii),14,SelReg,an(ii,:)');
    res(:,ii) = Y(15:end,ii) - Ysim(15:end);
    rss(ii) = 100*norm(res(:,ii))^2/norm(Y(15:end,ii))^2;
    R2(ii) = 1 - rss(ii)/100;
    mnse(ii) = mean(abs(res(:,ii))./(1+abs(Y(15:end,ii))));
end

for ii = 81:100
    [~,an(ii,:)] = PCparam(Q_val(ii-80,:)',5,INDX{1},TH.sim(:),options);
    AccIn = load(['MultiII_',num2str(ii),'.txt']);
    X(:,ii) = AccIn(201:N,2);
    VelOut = load(['MultiII_vel_',num2str(ii),'.txt']);
    Y(:,ii) = VelOut(201:N,2); 
    Ysim = narxsim(X(:,ii),Y(1:14,ii),14,SelReg,an(ii,:)');
    res(:,ii) = Y(15:end,ii) - Ysim(15:end);
    rss(ii) = 100*norm(res(:,ii))^2/norm(Y(15:end,ii))^2;
    R2(ii) = 1 - rss(ii)/100;
    mnse(ii) = mean(abs(res(:,ii))./(1+abs(Y(15:end,ii))));
end



figure(1)
subplot(4,8,[1:4 9:12]),plot(1:K,rss(1:K),'-bo')
legend({'Estimation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
% axis([1 20 0.825 1])
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
ylabel('R^2','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,8,[5:8 13:16]),plot(1:20,rss(K+1:100),'-rd')
legend({'Validation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
% set(gca,'Fontsize',9,'Fontname','TimesNewRoman','yticklabel',[])
% axis([1 20 0.825 1])
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'SF_est_vs_val']),close;
    if ispc
        result = eps2xxx([write_dir,'SF_est_vs_val.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'SF_est_vs_val.eps'],{'pdf'});
    end
end



%% Global PC-NARX R2 criterion [estimation vs validation set] -- Figures
% =========================================================================
clear,clc,close all
cd('D:\Data\FEM\MultistoreySF2D_II\')
load('MultiII_RandomExc_GA_rss.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

N = 1000;
K = 80;
load('InputVars.mat');
load('SF_RandomExc_PCNARX_lsqnonlin_normal_largePC.mat')
clear res
options.basis = 'hermi';
options.ortho = 'y';

for ii = 1:K
    [~,an(ii,:)] = PCparam(Q_est(ii,:)',5,INDX{1},TH.sim(:),options);
    AccIn = load(['MultiII_',num2str(ii),'.txt']);
    X(:,ii) = AccIn(201:N,2);
    VelOut = load(['MultiII_vel_',num2str(ii),'.txt']);
    Y(:,ii) = VelOut(201:N,2); 
    Ysim = narxsim(X(:,ii),Y(1:14,ii),14,SelReg,an(ii,:)');
    res(:,ii) = Y(15:end,ii) - Ysim(15:end);
    rss(ii) = 100*norm(res(:,ii))^2/norm(Y(15:end,ii))^2;
    R2(ii) = 1 - rss(ii)/100;
    mnse(ii) = mean(abs(res(:,ii))./(1+abs(Y(15:end,ii))));
end

N = 1000;
K = 49;
load('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\MultistoreySF2D_II\InputVarsVal.mat');

for ii = 1:49
    [~,an(ii,:)] = PCparam(Q_val(ii,:)',5,INDX{1},TH.sim(:),options);
    AccIn = load(['C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\MultistoreySF2D_II\MultiIIval_',num2str(ii),'.txt']);
    X(:,ii) = AccIn(201:N,2);
    VelOut = load(['C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\MultistoreySF2D_II\MultiIIval_vel_',num2str(ii),'.txt']);
    Y(:,ii) = VelOut(201:N,3); 
    Ysim = narxsim(X(:,ii),Y(1:14,ii),14,SelReg,an(ii,:)');
    res(:,ii) = Y(15:end,ii) - Ysim(15:end);
    rss(ii+80) = 100*norm(res(:,ii))^2/norm(Y(15:end,ii))^2;
    R2(ii) = 1 - rss(ii)/100;
    mnse(ii) = mean(abs(res(:,ii))./(1+abs(Y(15:end,ii))));
end

%%
figure(1)
subplot(4,8,[1:4 9:12]),plot(1:80,rss(1:80),'-bo')
legend({'Estimation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
% axis([1 20 0.825 1])
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
ylabel('R^2','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,8,[5:8 13:16]),plot(1:49,rss(81:129),'-rd')
legend({'Validation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
% set(gca,'Fontsize',9,'Fontname','TimesNewRoman','yticklabel',[])
% axis([1 20 0.825 1])
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'SF_est_vs_val']),close;
    if ispc
        result = eps2xxx([write_dir,'SF_est_vs_val.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'SF_est_vs_val.eps'],{'pdf'});
    end
end


%% #########################################################################


%% Steel Properties
% Properties                    Carbon     Alloy        Stainless     Tool
% ==========================================================================
% Density (1000 kg/m3)        7.85          7.85        7.75-8.1    7.72-8.0
% Elastic Modulus (GPa)       190-210       190-210     190-210     190-210
% Poisson's Ratio             0.27-0.3      0.27-0.3    0.27-0.3    0.27-0.3
% Tensile Strength (MPa)      276-1882      758-1882    515-827     640-2000
% Yield Strength (MPa)        186-758       366-1793    207-552     380-440
% ==========================================================================

%% ANSYS file output
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

%   SET   TIME/FREQ    LOAD STEP   SUBSTEP  CUMULATIVE
%     1  2.4250             1         1         1
%     2  4.9991             1         2         2
%     3  5.6430             1         3         3
%     4  8.5153             1         4         4
%     5  11.522             1         5         5
%     6  14.365             1         6         6
%     7  17.286             1         7         7