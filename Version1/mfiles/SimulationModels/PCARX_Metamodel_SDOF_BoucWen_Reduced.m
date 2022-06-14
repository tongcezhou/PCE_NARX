%% NARX regressors selection (GA - rss)
% =========================================================================
clear,clc,close all
if ispc
    cd(pwd);
    write_dir = [pwd,'\Figures\'];
else
    cd(pwd)
end
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')

%Number of samples
N = 1000;
% Simulation experiment
SimExp = 15;

% Load data
X = force(1:N,SimExp);
Y = veloc(1:N,SimExp);

% Model orders 
na = 2;     nb = 2;     nd = 0;

% Polynomial regressors
P = 3;      maxterms = 1; 
regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))'});

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

save('Results/SDOF_BoucWen_RandomExc_GA_rss_reduced.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')


%% Local NARX model estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
if ispc
    cd(pwd)
else
    cd(pwd)
end
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_GA_rss_reduced.mat')
%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = veloc(1:N,:);

% Model orders 
na = 2;     nb = 2;     nd = 0;
% Regressors
SelReg = SEstructure.removed(end-6:end);

clear *criteria* *res* theta*
orders = [na nd nb];
for k = 1:K
    disp(k)
    options.focus = 'prediction';
    options.method = 'ols';
    [thetaPE(:,k),NLresPE{k},PEcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    options.focus = 'simulation';
    options.PEinit = 'y';
    options.nlmethod = 'LM';
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
%     mnseSIM(k) = SIMcriteria{k}.mnse;
end

save('Results/SDOF_BoucWen_RandomExc_LocalNARX_lsqnonlin_reduced.mat','*criteria','NLres*','theta*','SelReg','orders')


%% PC-NARX estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
if ispc
    cd(pwd)
else
    cd(pwd)
end
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_LocalNARX_lsqnonlin_reduced.mat')
%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = veloc(1:N,:);

options.basis = 'legen';
% Basis indx
INDX{1} = [0    0;
     1     0;
     0     1;
     2     2];
     
% Total number of basis functions
B = size(INDX{1},1);
% Total polynomial degree
P = 4;
% Basis constraction: size(basis_i) = p x K
basis = cell(M,1);
Q = Q';
for m = 1:M
    basis{m} = PCbasis(Q(m,:),P,options);
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
% figure
% plot(thetaSIM,'-b'),hold on,plot(THij'*phiB,'--r');
options.focus = 'simulation';
options.PEinit = 'n';
options.nlmethod = 'LM';

% Model orders 
na = 2;     nb = 2;     nd = 0;
% [TH,res,criteria,output] = pcnarx(Y,X,Q,[na nd nb],INDX,SelReg,[],options);
[TH,res,criteria,output] = pcnarx(Y,X,Q,[na nd nb],INDX,SelReg,THij,options);
save('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin_reduced.mat','criteria','res','TH','output','orders','SelReg','INDX','THij','options')


%% =========================================================================
%% Figures
%% =========================================================================

%% Uncertain properties -- Figure
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
load('Results/SDOF_BoucWen_RandomExc_LHS.mat');

figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(6,1,1:3)
hl1 = line(1:K,Qtrans(:,1),'Color','b','Marker','o');
ax1 = gca;
set(ax1,'XColor','k','YColor','k','xlim',[0.5,20.5],'ylim',[0 1],'YTick',[0:0.2:1],...
        'xticklabel',[],'Fontname','TimesNewRoman','Fontsize',9)
ylabel('$\beta$','Fontsize',12,'Interpreter','Latex')
hleg1 = legend({'$\beta$'},'Interpreter','Latex','Fontsize',12);
set(hleg1,'Position',[0.105 0.94 0.4 0.05]);
ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','right',...
           'Color','none','XColor','k','YColor','k','XTick',[0:2:20],'YTick',[0:20:100],...
           'xlim',[0.5,20.5],'ylim',[0 100],'Fontname','TimesNewRoman','Fontsize',9);
hl2 = line(1:K,Qtrans(:,2),'Color','r','Marker','d','Linestyle','--','Parent',ax2);
box on
hleg2 = legend({'$\sigma_{\!\! _x}$'},'Interpreter','Latex','Fontsize',12);
set(hleg2,'Position',[0.53 0.94 0.4 0.05]);
xlabel('Simulation experiment number','Fontsize',10,'Fontname','TimesNewRoman')
ylabel('$\sigma_{\!\! _x}$','Fontsize',12,'Interpreter','Latex','Rotation',270)
grid on
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_uncprops']),close;
     if ispc
        result = eps2xxx([write_dir,'boucwen_uncprops.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
         result = eps2xxx([write_dir,'boucwen_uncprops.eps'],{'pdf'});
     end
end


%% Simulated responses -- Figures
% =========================================================================
clear,clc,close all
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

minexp = 2;
maxexp = 15;
N = 1000;

figure(1)
clr = colormap(lines(2));
hold on
figindx = [1 1 2 2];
colindx = [1 2 1 2];
linindx = [1 1 1 1];
subplot(2,2,1),plot(dspl(1:N,minexp),resforce(1:N,minexp),'-',...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
hold on
str = {['$\beta = ',num2str(roundn(Qtrans(minexp,1),-2)),', \sigma_x = ',num2str(round(Qtrans(minexp,2))),' (N)$']};
legend(str,'Fontsize',9,'Location','Northoutside','Interpreter','Latex')
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Force (N)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,2,5),plot(0.005*(0:N-1),dspl(1:N,minexp),...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
set(gca,'xticklabel',[])
ylabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,2,7),plot(0.005*(0:N-1),resforce(1:N,minexp),...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
grid on;% axis tight
ylabel('Force (N)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(2,2,2),plot(dspl(1:N,maxexp),resforce(1:N,maxexp),'-',...
    'Linewidth',linindx(1),'Color',clr(colindx(2),:))
hold on
grid on
str = {['$\beta = ',num2str(roundn(Qtrans(maxexp,1),-2)),', \sigma_x = ',num2str(round(Qtrans(maxexp,2))),' (N)$']};
legend(str,'Fontsize',9,'Location','Northoutside','Interpreter','Latex')
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Force (kN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(4,2,6),plot(0.005*(0:N-1),dspl(1:N,maxexp),...
    'Linewidth',linindx(2),'Color',clr(colindx(2),:))
set(gca,'xticklabel',[])
% ylim([-1.5 1.5])
grid on
ylabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(4,2,8),plot(0.005*(0:N-1),resforce(1:N,maxexp),...
    'Linewidth',linindx(2),'Color',clr(colindx(2),:))
% axis tight
% ylim([-15 15])
grid on
ylabel('Force (N)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_randexci']),close;
     if ispc
     result = eps2xxx([write_dir,'boucwen_randexci.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
         result = eps2xxx([write_dir,'boucwen_randexci.eps'],{'pdf'});
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
end
load('Results/SDOF_BoucWen_RandomExc_GA_rss_reduced.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

NoReg = length(GAreg);
close all

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(2,1,1),
plot(1:NoReg,SEstructure.values(1:NoReg),'-o','Markersize',5,'Linewidth',1)
hold on 
plot([1 NoReg],SEstructure.initvalue*[1 1],'--r','Markersize',5,'Linewidth',1)
set(gca,'yscale','linear','Fontsize',7,'xtick',[1:NoReg],'xticklabel',[])
% ylabel('NSSE (%)','Fontangle','normal','Fontsize',9)
for ii = 1:NoReg
    if strcmp(SEstructure.removed{ii},'ones(length(t),size(Y,2))')
        txtstr = '$\mbox{const. term}$'; 
    elseif length(SEstructure.removed{ii})>13
        if strcmp(SEstructure.removed{ii}(12),'1');
            txtstr = ['$',lower(SEstructure.removed{ii}(2)),'[',SEstructure.removed{ii}(4:6),']','\cdot |y[t-1]|$'];
        else
            txtstr = ['$',lower(SEstructure.removed{ii}(2)),'[',SEstructure.removed{ii}(4:6),']',SEstructure.removed{ii}(11:12),'\cdot |y[t-1]|$'];
        end
    else
        if strcmp(SEstructure.removed{ii}(12),'1');
            txtstr = ['$',lower(SEstructure.removed{ii}(2)),'[',SEstructure.removed{ii}(4:6),']','$'];
        else
            txtstr = ['$',lower(SEstructure.removed{ii}(2)),SEstructure.removed{ii}(11:12),'[',SEstructure.removed{ii}(4:6),']','$'];
        end
    end
    text(ii,-5,txtstr,'Rotation',65,'Horizontalalignment','right','Interpreter','Latex','Fontsize',8)
end
grid on
axis([0.9 NoReg+0.1 0 105])
grid on
text(9.5,-60,'Regressors dropped','Fontangle','normal','Fontsize',10,'Horizontalalignment','center','Fontname','TimesNewRoman')
ylabel('NSSE (%)','Fontangle','normal','Fontsize',10,'Fontname','TimesNewRoman')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_MSS_StageA']),close;
     result = eps2xxx([write_dir,'boucwen_MSS_StageA.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% Parameters expansion -- GA -- Figure
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
load('Results/SDOF_Boucwen_RandomExc_LHS.mat')
load('Results/SDOF_Boucwen_RandomExc_LocalNARX_lsqnonlin_reduced.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

clear basis options
options.basis = 'legen';
options.criterion = 'ar2';
Pmax = 4;

options.GA.PopulationSize = 100;
[A,theta,BASISstructure,criteria,PHI,GAoutput] = gaPC(thetaSIM',Q',Pmax,1,options);
NoB = length(BASISstructure.values);

close all
figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(211),plot(1:NoB,BASISstructure.values(1:NoB),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 NoB],BASISstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
set(gca,'yscale','linear','Fontsize',8,'xtick',[1:NoB],'xticklabel',[])
for ii = 1:NoB
    txtstr = ['$[',num2str(BASISstructure.removed(ii,1)),',',num2str(BASISstructure.removed(ii,2)),']$'];
    text(ii,-0.15,txtstr,'Rotation',25,'Fontsize',7,'Horizontalalignment','center','Interpreter','Latex','Fontsize',10)
end
axis([1 NoB -0.1 0.4])
grid on
text(NoB/2+0.5,-0.225,'Basis functions dropped','Fontangle','normal','Fontsize',10,'Horizontalalignment','center','Fontname','TimesNewRoman')
ylabel('$\overline{R^2_{\mbox{adj}}}$','Fontangle','normal','Fontsize',11,'Interpreter','Latex')
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'boucwen_MSS_StageB']),close;    
    if ispc
        result = eps2xxx([write_dir,'boucwen_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'boucwen_MSS_StageB.eps'],{'pdf'});
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
print_to_eps = 'n';

K = 20;
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin_reduced.mat')
clear res
options.basis = 'legen';
Pmax = 4;
for ii = 1:K
    [~,an(ii,:)] = PCparam(Q(ii,:)',Pmax,INDX{1},TH.sim(:),options);
    X = force(1:1000,ii);
    Y = veloc(1:1000,ii);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    res(:,ii) = Y(3:end) - Ysim(3:end);
    rss(ii) = 100*norm(res(:,ii))^2/norm(Y(3:end))^2;
    R2(ii) = 1 - rss(ii)/100;
    mnse(ii) = mean(abs(res(:,ii))./(1+abs(Y(3:end))));
end

clear res X Y Q
load('Results/SDOF_BoucWen_RandomExc_validation.mat')
for ii = 1:K
    [~,an(ii,:)] = PCparam(Q(ii,:)',Pmax,INDX{1},TH.sim(:),options);
    X = force(1:1000,ii);
    Y = veloc(1:1000,ii);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    res(:,ii) = Y(3:end) - Ysim(3:end);
    rss(ii+K) = 100*norm(res(:,ii))^2/norm(Y(3:end))^2;
    R2(ii+K) = 1 - rss(ii+K)/100;
    mnse(ii+K) = mean(abs(res(:,ii))./(1+abs(Y(3:end))));
end

figure(1)
subplot(4,8,[1:4 9:12]),plot(1:K,rss(1:K),'-bo')
legend({'Estimation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
% axis([1 20 0.825 1])
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
ylabel('R^2','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,8,[5:8 13:16]),plot(1:K,rss(K+1:2*K),'-rd')
legend({'Validation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
% set(gca,'Fontsize',9,'Fontname','TimesNewRoman','yticklabel',[])
% axis([1 20 0.825 1])
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'boucwen_est_vs_val']),close;
    if ispc
        result = eps2xxx([write_dir,'boucwen_est_vs_val.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'boucwen_est_vs_val.eps'],{'pdf'});
    end
end




%% Global PC-NARX validation runs -- Figures
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Desktop\PCNARXmodels_JournalPaper\mfiles\SimulationModels')
    write_dir = 'C:\Users\sminas\Desktop\PCNARXmodels_JournalPaper\mfiles\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

K = 20;
time = (0:999)*0.005;
load('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin_reduced.mat')
load('Results/SDOF_BoucWen_RandomExc_validation.mat')
clear res
options.basis = 'legen';
Pmax = 4;
SimExp1 = 6;
SimExp2 = 20;
ii=1;
for k = [SimExp1 SimExp2]
    [~,an(ii,:)] = PCparam(Q(k,:)',Pmax,INDX{1},TH.sim(:),options);
    X(:,ii) = force(1:1000,k);
    Y(:,ii) = veloc(1:1000,k);
    Ysim(:,ii) = narxsim(X(:,ii),Y(1:2,ii),2,SelReg,an(ii,:)');
    res(:,ii) = Y(3:end,ii) - Ysim(3:end,ii);
    rss(ii) = 100*norm(res(:,ii))^2/norm(Y(3:end,ii))^2;
    R2(ii) = 1 - rss(ii)/100;
    mnse(ii) = mean(abs(res(:,ii))./(1+abs(Y(3:end,ii))));
    ii=ii+1;
end

figure(1)
subplot(2,1,1),plot(time,Y(:,1),'-b','Linewidth',1)
hold on
plot(time,Ysim(:,1),'--r','Linewidth',1)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$\beta = ',num2str(roundn(Qtrans(SimExp1,1),-2)),', \sigma_x = ',num2str(round(Qtrans(SimExp1,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',12)
legend(str,'Fontsize',11,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',9,'Fontname','TimesNewRoman','xticklabel',[])
% set(h1,'Position',[0.3 0.925 0.4 0.05])
ylabel('Velocity (m/s)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(2,1,2),plot(time,Y(:,2),'-b','Linewidth',1)
hold on
plot(time,Ysim(:,2),'--r','Linewidth',1)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$\beta = ',num2str(roundn(Qtrans(SimExp2,1),-2)),', \sigma_x = ',num2str(round(Qtrans(SimExp2,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',12)
legend(str,'Fontsize',11,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Velocity (m/s)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'boucwen_randexci_val']),close;
    if ispc
        result = eps2xxx([write_dir,'boucwen_randexci_val.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'boucwen_randexci_val.eps'],{'pdf'});
    end
end

% Figure with errors
figure(1)
clr = colormap(lines(2));
hold on
subplot(8,4,[1 2 5 6 9 10]),plot(time,Y(:,1),'-b','Linewidth',0.75)
hold on
plot(time,Ysim(:,1),'--r','Linewidth',0.75)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$\beta = ',num2str(roundn(Qtrans(SimExp1,1),-2)),', \sigma_x = ',num2str(round(Qtrans(SimExp1,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',8)
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','xticklabel',[],'ylim',[-10 10])
% set(h1,'Position',[0.3 0.925 0.4 0.05])
ylabel('Velocity (m/s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on

subplot(8,4,[13 14 17 18]),plot(time(3:end),res(:,1),'-k','Linewidth',0.75)
str = {'PC-NARX metamodel simulation error'};
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','ylim',[-1 1])
ylabel('Error (m)','Fontsize',8,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on

subplot(8,4,[3 4 7 8 11 12]),plot(time,Y(:,2),'-b','Linewidth',0.75)
hold on
plot(time,Ysim(:,2),'--r','Linewidth',0.75)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$\beta = ',num2str(roundn(Qtrans(SimExp2,1),-2)),', \sigma_x = ',num2str(round(Qtrans(SimExp2,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',8)
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','xticklabel',[],'ylim',[-10 10])
grid on

subplot(8,4,[15 16 19 20]),plot(time(3:end),res(:,2),'-k','Linewidth',0.75)
str = {'PC-NARX metamodel simulation error'};
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','ylim',[-1 1])
xlabel('Time (s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'boucwen_randexci_val_error']),close;
    if ispc
        result = eps2xxx([write_dir,'boucwen_randexci_val_error.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'boucwen_randexci_val_error.eps'],{'pdf'});
    end
end





% Figure with errors
figure(1)
clr = colormap(lines(2));
hold on
subplot(18,1,1:6),plot(time,Y(:,1),'-b','Linewidth',0.75)
hold on
plot(time,Ysim(:,1),'--r','Linewidth',0.75)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$\beta = ',num2str(roundn(Qtrans(SimExp1,1),-2)),', \sigma_x = ',num2str(round(Qtrans(SimExp1,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',8)
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','xticklabel',[],'ylim',[-10 10])
% set(h1,'Position',[0.3 0.925 0.4 0.05])
ylabel('Velocity (m/s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on

subplot(18,1,7:8),plot(time(3:end),res(:,1),'-k','Linewidth',0.75)
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','ylim',[-1 1])
ylabel('Error (m/s)','Fontsize',8,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on

subplot(18,1,11:16),plot(time,Y(:,2),'-b','Linewidth',0.75)
hold on
plot(time,Ysim(:,2),'--r','Linewidth',0.75)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$\beta = ',num2str(roundn(Qtrans(SimExp2,1),-2)),', \sigma_x = ',num2str(round(Qtrans(SimExp2,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',8)
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','xticklabel',[],'ylim',[-10 10])
ylabel('Velocity (m/s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on

subplot(18,1,17:18),plot(time(3:end),res(:,2),'-k','Linewidth',0.75)
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','ylim',[-1 1])
ylabel('Error (m/s)','Fontsize',8,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'boucwen_randexci_val_error']),close;
    if ispc
        result = eps2xxx([write_dir,'boucwen_randexci_val_error.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'boucwen_randexci_val_error.eps'],{'pdf'});
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

load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin_reduced.mat')
load('Results/SDOF_BoucWen_RandomExc_LocalNARX_lsqnonlin_reduced.mat')
PCBASES = INDX{1};
NLREGRESSORS = SelReg;
options.basis = 'legen';
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
surf(linspace(0,100,NoI2),linspace(0,1,NoI1),reshape(an_interp(:,3),NoI1,NoI2))
shading interp
% hold on
% plot3(Qtrans(:,2),Qtrans(:,1),thetaSIM(2,:)','o')
set(gca,'Fontsize',10)
ylabel('$\beta$','Interpreter','Latex','Fontsize',20)
xlabel('$\sigma_{x}$','Interpreter','Latex','Fontsize',20)
zlabel(['$\hat{\theta}_{y[t-2]\cdot |y[t-1]|}\quad\ \ (\xi)$'],'Interpreter','Latex','Fontsize',20)
grid on
box on
% zlim([-0.07 0.01])
view([45 25])
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_surf_a']),close;
     if ispc
        result = eps2xxx([write_dir,'boucwen_surf_a.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
         result = eps2xxx([write_dir,'boucwen_surf_a.eps'],{'pdf'});
     end
end


figure(4)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
surf(linspace(0,100,NoI2),linspace(0,1,NoI1),reshape(an_interp(:,7),NoI1,NoI2))
shading interp
set(gca,'Fontsize',10)
ylabel('$\beta$','Interpreter','Latex','Fontsize',20)
xlabel('$\sigma_x$','Interpreter','Latex','Fontsize',20)
zlabel(['$\hat{\theta}_{y[t-2]} \ \ (\xi)$'],'Interpreter','Latex','Fontsize',20)
grid on
box on
% zlim([1.99 2.01])
view([45 25])
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_surf_b']),close;
     if ispc
        result = eps2xxx([write_dir,'boucwen_surf_b.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
         result = eps2xxx([write_dir,'boucwen_surf_b.eps'],{'pdf'});
     end
end