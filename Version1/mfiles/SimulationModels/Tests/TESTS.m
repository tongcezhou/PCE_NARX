%% PC-NARX theoretical model simulations
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
RandomExc = 1;
if RandomExc == 0
    load('SDOF_BoucWen_RandomExc_FineGrid.mat')
    N = 1000;
else
    load('SDOF_BoucWen_SweepSine_FineGrid.mat')
    N = 1000;
end
% Input-Output data
Y = veloc(1:N,:);
X = force(1:N,:);

% Model orders
na = 2;     nb = 2;     nd = 0;
TheorRegstr{1} = '(Y(t-1,:).^1)';
TheorRegstr{2} = '(Y(t-2,:).^1)';
TheorRegstr{3} = '(X(t-1,:).^1)';
TheorRegstr{4} = '(X(t-2,:).^1)';
TheorRegstr{5} = '(X(t-1,:).^3).*abs(Y(t-1,:))';
TheorRegstr{6} = '(Y(t-1,:).^1).*(X(t-1,:).^2).*abs(Y(t-1,:))';
TheorRegstr{7} = '(Y(t-2,:).^1).*(X(t-1,:).^2).*abs(Y(t-1,:))';
TheorRegstr{8} = '(Y(t-1,:).^2).*(X(t-1,:).^1).*abs(Y(t-1,:))';
TheorRegstr{9} = '(Y(t-2,:).^2).*(X(t-1,:).^1).*abs(Y(t-1,:))';
TheorRegstr{10} = '(Y(t-1,:).^1).*(Y(t-2,:).^1).*(X(t-1,:).^1).*abs(Y(t-1,:))';
TheorRegstr{11} = '(Y(t-1,:).^3).*abs(Y(t-1,:))';
TheorRegstr{12} = '(Y(t-2,:).^3).*abs(Y(t-1,:))';
TheorRegstr{13} = '(Y(t-1,:).^2).*(Y(t-2,:).^1).*abs(Y(t-1,:))';
TheorRegstr{14} = '(Y(t-1,:).^1).*(Y(t-2,:).^2).*abs(Y(t-1,:))';

% Random input variables
Fmax = 100:100:500;
BETA = 0:0.2:1;
m = 1;
ki = 50;
Ts = 0.005;

k=0;
for i = 1:indxb
    for j = 1:indxF
        k = k+1;
        beta = BETA(i);
        theta = [2-(Ts^2)*ki/m;
            -1;
            Ts/m;
            -Ts/m;
            (Ts^2)*beta/((ki^2)*m);
            -3*Ts*beta/(ki^2);
            3*Ts*beta/(ki^2);
            3*beta*m/(ki^2);
            3*beta*m/(ki^2);
            6*beta*m/(ki^2);
            -(m^2)*beta/((ki^2)*Ts);
            (m^2)*beta/((ki^2)*Ts);
            3*beta*(m^2)/((ki^2)*Ts);
            -3*beta*(m^2)/((ki^2)*Ts)];
        Ysim(:,k) = narxsim(X(:,k),Y(1:na,k),na,TheorRegstr([1:14]),theta([1:14]));
    end
end



%% PC-NARX model structure selection (GA)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
RandomExc = 0;
if RandomExc == 0
    load('SDOF_BoucWen_RandomExc_FineGrid.mat')
else
    load('SDOF_BoucWen_SweepSine_FineGrid.mat')
end
N = 1000;
% Input-Output data
Y = veloc(1:N,24);
X = force(1:N,24);

% Nonlinear polynomial regressors
P = 3;      q = 1;      maxterms = 2;
% Model orders
na = 2;     nb = 2;     nd = 0;
options.warnings = 'off';
options.GA.PopulationSize = 100 ;
options.maxsize = 1e8;

for P = 1:3
    for maxterms = 1:3
        for absfun = 0:1
            disp([P,maxterms,absfun])
            % Nonlinear regressors (complete vector)
            if absfun == 1
                regressors = polyregNARX([na nd nb],P,q,maxterms,{'abs(Y(t-1,:))'});
            else
                regressors = polyregNARX([na nd nb],P,q,maxterms);
            end
            options.focus = 'simulation';
            [GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);
            rss_sss(P,maxterms,absfun+1) = criteria.rss_sss;
            mse(P,maxterms,absfun+1) = mean(RES.^2);
            sne(P,maxterms,absfun+1) = sum(abs(RES)./(1+abs(Y)));
            nparam(P,maxterms,absfun+1) = length(GAreg);
            options.focus = 'prediction';
            [TH,res,criteria,output] = narx(Y,X,[na nd nb],regressors,[],[],[],[],options);
            Ysim = narxsim(X,Y(1:na),na,regressors,TH);
            E = Y-Ysim;
            rss_sss_all(P,maxterms,absfun+1) = 100*norm(E)^2/norm(Y)^2;
            mse_all(P,maxterms,absfun+1) = mean(E.^2);
            sne_all(P,maxterms,absfun+1) = sum(abs(E)./(1+abs(Y)));
            nparam_all(P,maxterms,absfun+1) = length(regressors);
        end
    end
end






%% Bouc-Wen model NARX regressors selection (GA)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('SDOF_BoucWen_RandomExc_LHS.mat')

options.maxsize = 1e7;
options.focus = 'prediction';
N = 1000;

X = force(1:N,:);
Y = veloc(1:N,:);

% Structure selection for all simulated responses
na = 2;     nb = 2;     nd = 0;
q = 1;
options.warnings = 'off';
options.GA.PopulationSize = 100 ;
options.maxsize = 1e8;
P = 3; maxterms = 2; absfun = 1;


for k = 1:K
    % Nonlinear regressors (complete vector)
    if absfun == 1
        regressors = polyregNARX([na nd nb],P,q,maxterms,{'abs(Y(t-1,:))'});
    else
        regressors = polyregNARX([na nd nb],P,q,maxterms);
    end
    options.focus = 'simulation';
    [GAreg{k},indx{k},RES,criteria,GAoutput] = gaNARX(Y(:,k),X(:,k),[na nd nb],regressors,options);
    rss_sss(k) = criteria.rss_sss;
    mse(k) = mean(RES.^2);
    sne(k) = sum(abs(RES)./(1+abs(Y(:,k))));
    nparam(k) = length(GAreg);
    options.focus = 'prediction';
    [TH,res,criteria,output] = narx(Y(:,k),X(:,k),[na nd nb],regressors,[],[],[],[],options);
    Ysim = narxsim(X(:,k),Y(1:na,k),na,regressors,TH);
    E = Y(:,k)-Ysim;
    rss_sss_all(k) = 100*norm(E)^2/norm(Y(:,k))^2;
    mse_all(k) = mean(E.^2);
    sne_all(k) = sum(abs(E)./(1+abs(Y(:,k))));
    nparam_all(k) = length(regressors);
end



%% %% Structure selection for all simulated responses -- [results]
freqbins = zeros(size(regressors));
for k=1:K
    for j=1:length(indx{k})
        freqbins(indx{k}(j)) = freqbins(indx{k}(j))+1;
    end
end

%% Model orders
na = 2;     nb = 2;     nd = 0;
q = 1;
options.warnings = 'off';
options.GA.PopulationSize = 100 ;
options.maxsize = 1e8;

for P = 1:3
    for maxterms = 1:3
        for absfun = 0:1
            disp([P,maxterms,absfun])
            % Nonlinear regressors (complete vector)
            if absfun == 1
                regressors = polyregNARX([na nd nb],P,q,maxterms,{'abs(Y(t-1,:))'});
            else
                regressors = polyregNARX([na nd nb],P,q,maxterms);
            end
            options.focus = 'simulation';
            [GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);
            rss_sss(P,maxterms,absfun+1) = criteria.rss_sss;
            mse(P,maxterms,absfun+1) = mean(RES.^2);
            sne(P,maxterms,absfun+1) = sum(abs(RES)./(1+abs(Y)));
            nparam(P,maxterms,absfun+1) = length(GAreg);
            options.focus = 'prediction';
            [TH,res,criteria,output] = narx(Y,X,[na nd nb],regressors,[],[],[],[],options);
            Ysim = narxsim(X,Y(1:na),na,regressors,TH);
            E = Y-Ysim;
            rss_sss_all(P,maxterms,absfun+1) = 100*norm(E)^2/norm(Y)^2;
            mse_all(P,maxterms,absfun+1) = mean(E.^2);
            sne_all(P,maxterms,absfun+1) = sum(abs(E)./(1+abs(Y)));
            nparam_all(P,maxterms,absfun+1) = length(regressors);
        end
    end
end


options.focus = 'prediction';
options.tolfun = 1;
[SEstructure,IndxRem] = SEreduction(Y(:,k),X(:,k),[na nd nb],GAreg{k},options);

%%
clear *criteria* *res* theta*
for k=1:K
    X = force(1:N,k);
    Y = veloc(1:N,k);
    options.focus = 'prediction';
    [thetaPE{k},NLresPE{k},PEcriteria{k}] = narx(Y,X,[na nd nb],GAreg(find(IndxRem)),[],[],[],[],options);
    options.focus = 'simulation';
    options.PEinit = 'y';
    [thetaSIM{k},NLresSIM{k},SIMcriteria{k}] = narx(Y,X,[na nd nb],GAreg(find(IndxRem)),[],[],[],[],options);
end

for k =1:K
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    THpe(:,k) = thetaPE{k};
    THsim(:,k) = thetaSIM{k};
end

%% Space frame model - NARX regressors selection (GA)
% =========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SpaceFrameFEMreduced';
else
    local_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/SpaceFrameFEMreduced';
end

cd(local_dir)

% Number of experiments
K = 20; 
% Number of samples
T = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 2; 
% Sampling period
Ts =  0.02;

% Input-Output data
[Y,X] = deal(zeros(T,K));
for i = 1:K
    out = load(['Responses/SFres',num2str(i),'.txt']);
    Y(:,i) = out(1:T,13);
    X(:,i) = out(1:T,2);
end
load('Responses/InputVars.mat','Q');

%%
options.criterion = 'all';
options.NFFT = 512;
options.wind = hamming(256);
options.nsamples = 256;
options.overlap = 80;
options.robust = 'n';
options.percentile = 0;
options.alpha = 0.05;
for k = 1:K
    [NLI(k),TNLI(k),MAXBIC(k),SUMBIC(k),MeanCoh(k)]=  NLmeasures(Y(:,k),X(:,k),options);
    % Cross-corellation(z,z^2)
    CCF(:,k) = ccf(Y(:,k)-mean(Y(:,k)),(Y(:,k)-mean(Y(:,k))).^2,100,0.8);
end
figure(1)
subplot(611),plot(1:K,NLI,'-bo')
subplot(612),plot(1:K,TNLI,'-bo')
subplot(613),plot(1:K,MAXBIC,'-bo')
subplot(614),plot(1:K,SUMBIC,'-bo')
subplot(615),plot(1:K,MeanCoh,'-bo')
subplot(616),plot(1:K,sum(CCF),'-bo')


figure(K+2)
for k = 1:K
    [bic,waxis] = bicoherx (X(:,k),X(:,k),Y(:,k),128,[],128,50,'y',25);
    subplot(4,5,k),surf(bic)
    shading interp
    caxis([0 1])
    zlim([0 1])
end

for k=1:K; 
    der(:,k) = diff(resforce(:,k))./diff(dspl(:,k));
    mder(:,k) = der(:,k) - mean(der(:,k));
    sder(k) = sum(der(:,k));
end


for k=1:K; 
    figure(k),
    subplot(211),plot(dspl(:,k),resforce(:,k)),
    subplot(212),plot(diff(resforce(:,k))./diff(dspl(:,k))),
end

%%
ii = input('Given run with stronger nonlinearities: ');
X = X(1:T,ii);
Y = Y(1:T,ii);
 
na = 20;     nb = 20;     nd = 0;
q = 1;
options.warnings = 'off';
options.GA.PopulationSize = 100 ;
options.maxsize = 1e8;
maxterms = 1;
P = 3;
regressors = polyregNARX([na nd nb],P,q,maxterms);
% regressors = polyregNARX([na nd nb],P,q,maxterms,{'abs(Y(t-1,:))';'abs(Y(t-2,:))';...
%      'abs(Y(t-3,:))';'abs(Y(t-4,:))';'abs(Y(t-5,:))';'abs(Y(t-6,:))';'abs(Y(t-7,:))';...
%      'abs(Y(t-8,:))';'abs(Y(t-9,:))';'abs(Y(t-10,:))'});
%      'abs(Y(t-11,:))';'abs(Y(t-12,:))';...
%      'abs(Y(t-13,:))';'abs(Y()t-14,:))';'abs(Y(t-15,:))';'abs(Y(t-16,:))';'abs(Y(t-17,:))';...
%      'abs(Y(t-18,:))';'abs(Y(t-19,:))';'abs(Y(t-20,:))'});
options.focus = 'simulation';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);


%% 
na = 10;     nb = 10;     nd = 0;
q = 1;
options.warnings = 'off';
options.GA.PopulationSize = 200 ;
options.maxsize = 1e8;
maxterms = 2;
P = 3;

regressors = polyregNARX([na nd nb],P,q,maxterms);
% regressors = polyregNARX([na nd nb],P,q,maxterms,{'abs(Y(t-1,:))';'abs(Y(t-2,:))';...
%       'abs(Y(t-3,:))';'abs(Y(t-4,:))';'abs(Y(t-5,:))'});
% regressors =regressors([1:42 84:329 535:end]);%43:83 330:534
options.focus = 'simulation';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);

%%
options.focus = 'prediction';
options.Nr = length(GAreg)-1;
[SEstructure,IndxRem,IndxRmv] = SEreduction(Y,X,[na nd nb],GAreg,options);

%% 

clear *criteria* *res* theta*
for k=1:K
    disp(k)
    out = load(['Responses/SFres',num2str(k),'.txt']);
    Y = out(1:T,13);
    X = out(1:T,2);
    options.focus = 'prediction';
    [thetaPE{k},NLresPE{k},PEcriteria{k}] = narx(Y,X,[na nd nb],GAreg(find(IndxRem)),[],[],[],[],options);
    options.focus = 'simulation';
    options.PEinit = 'y';
    [thetaSIM{k},NLresSIM{k},SIMcriteria{k}] = narx(Y,X,[na nd nb],GAreg(find(IndxRem)),[],[],[],[],options);
end

for k =1:K
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    THpe(:,k) = thetaPE{k};
    THsim(:,k) = thetaSIM{k};
end

    
%% PC basis
options.maxsize = 1e7;
options.focus = 'prediction';
options.criterion = 'bic';
options.nlreg = 'n';
options.basis = 'legen';
options.common = 'y';
options.GA.PopulationSize = 100;
options.GA.Display = 'iter';
options.GA.PlotFcns = {@gaplotbestf,@gaplotbestindiv};
[RegGA, IndxGA,thetaij,res,criteria,GAoutput]  = gaPCNARX(dspl(1:N,:),force(1:N,:),Q',[na nd nb],regressors(regindx(1:23)),3,1,options);

%
clear INDX
INDX{1} = combinations(2,3,1);
options.focus = 'prediction';
options.method = 'wls';
options.wls_iter = 20;
[THij,res,criteria,output] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],IndxGA,RegGA,[],options);
% SE estimation
options.focus = 'simulation';
options.PEinit = 'n';
[THij_sim,res_sim,criteria_sim] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],INDX,regressors(regindx(1:23)),[],options);


% Start initial vector from the previous estimated model

for j = 1:length(regindx)
    disp(j)
    options.method = 'wls';
    options.wls_iter = 5;
    options.focus = 'prediction';
    [thetaNL,NLresPE{j},PEcriteria] = narx(Y,X,[na nd nb],regressors(regindx(1:j)),[],[],[],[],options);
    rssPE(j) = PEcriteria.rss_sss;
    bicPE(j) = PEcriteria.bic;
    options.PEinit = 'n';
    options.focus = 'simulation';
    if j == 1
        [thetaNLSIM,NLresSIM{j},SIMcriteria] = narx(Y,X,[na nd nb],regressors(regindx(1:j)),[],[],[],thetaNL,options);
    else
        [thetaNLSIM,NLresSIM{j},SIMcriteria] = narx(Y,X,[na nd nb],regressors(regindx(1:j)),[],[],[],[thetaNLSIM ;0],options);
    end
    rssSIM(j) = SIMcriteria.rss_sss;
    %     Ysim = narxsim(X,Y(1:na),na,regressors(regindx(1:j)),thetaNL);
    %     NLresSIM{j} = Y - Ysim;
    %     rssSIM(j) = 100*norm(NLresSIM{j})^2/norm(Y)^2;
end

figure
plot(1:30,rssPE(1:30),'-bo',1:30,rssSIM(1:30),'-rd')
set(gca,'yscale','log')



%% 
% GA optimization
options.warnings = 'off';
options.GA.PopulationSize = 20;
options.focus = 'prediction';
options.maxsize = 1e8;
options.nlreg = 'y';
options.criterion = 'se';
% [RegGA,indx,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,P,maxterms,options);
[REGSTR, INDX,thetaij,res,criteria,GAoutput] = gaPCNARX(Y,X,Q',[na nd nb],regressors,3,1,options);
%%

options.focus = 'prediction';
options.method = 'wls';
options.wls_iter = 5;
[TH,res,criteria,output] = pcnarx(Y,X,Q',[na nd nb],INDX,REGSTR,[],options);

options.focus = 'simulation';
options.method = 'wls';
options.maxsize = 1e9;
options.wls_iter = 5;
[TH,res,criteria,output] = pcnarx(Y,X,Q',[na nd nb],INDX,REGSTR,[],options);


[res,criteria] = pcnarxSIM(TH.wls,Y,X,Q',[na nd nb],INDX,REGSTR,options);




%% SE reduction
P = 3;      q = 1;      maxterms = 1;
% Model orders
na = 2;     nb = 2;     nd = 0;
X = force(1:500,13);
Y = veloc(1:500,13);
regressors = polyregNARX([na nd nb],P,q,maxterms,{'abs(Y(t-1,:))'});
options.focus = 'prediction';
options.PEinit = 'n';
options.maxsize = 1e9;
[SEstructure,CURreg] = SEreduction(Y,X,[na nd nb],regressors,options);



%% Uncrosscorrelatedeness (estimated covariance matrix) --- OK
% =========================================================================
RES = reshape(res_sim,K,N)';

figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
surf(cov(RES))
clr = colormap(jet);
%colormap(flipud(clr));
shading interp
set(gca,'Fontsize',9,'projection','Perspective')
view([60 75])
% xlim([1 50])
% ylim([1 50])
% zlim([-0.005 0.025])
text(100,-40,0,'row index $k_1$','Interpreter','Latex','Fontsize',14,'Rotation',-50,'Horizontalalignment','center')
text(240,80,0,'column index $k_2$','Interpreter','Latex','Fontsize',14,'Rotation',22,'Horizontalalignment','center')
zlabel('cov($\hat{e}_{k_1}\; , \hat{e}_{k_2} \;$)','Interpreter','Latex','Fontsize',14)
grid on



%% Parameter surfaces
clear Qsurf
NoI1 = 200;
NoI2 = 200;
aux = linspace(-1,1,NoI1);
Qsurf(1,:) = kron(ones(1,NoI2),aux);
Qsurf(2,:) = kron(aux,ones(1,NoI2));

[~,an_interp] = PCparam(Qsurf,1,INDX{1},THij_sim,options);

figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(2,2,1),surf(linspace(0,100,NoI2),linspace(0,100,NoI1),reshape(an_interp(:,1),NoI1,NoI2))
shading interp
set(gca,'Fontsize',6)
ylabel('$k_nl$ (\%)','Interpreter','Latex','Fontsize',7)
xlabel('$F_{\max}$ (\%)','Interpreter','Latex','Fontsize',7)
zlabel(['$a_',num2str(1),'(\xi)$'],'Interpreter','Latex','Fontsize',7)
grid on
box on
zlim([1.875 1.925])
subplot(2,2,2),surf(linspace(0,100,NoI2),linspace(0,100,NoI1),reshape(an_interp(:,2),NoI1,NoI2))
shading interp
set(gca,'Fontsize',6)
ylabel('$k_nl$ (\%)','Interpreter','Latex','Fontsize',7)
xlabel('$F_{\max}$ (\%)','Interpreter','Latex','Fontsize',7)
zlabel(['$a_',num2str(1),'(\xi)$'],'Interpreter','Latex','Fontsize',7)
grid on
box on
zlim([-1.005 -.995])
subplot(2,2,3),surf(linspace(0,100,NoI2),linspace(0,100,NoI1),reshape(an_interp(:,3),NoI1,NoI2))
shading interp
set(gca,'Fontsize',6)
ylabel('$k_nl$ (\%)','Interpreter','Latex','Fontsize',7)
xlabel('$F_{\max}$ (\%)','Interpreter','Latex','Fontsize',7)
zlabel(['$a_',num2str(1),'(\xi)$'],'Interpreter','Latex','Fontsize',7)
grid on
box on
zlim([9e-5 10e-5])
subplot(2,2,4),surf(linspace(0,100,NoI2),linspace(0,100,NoI1),reshape(an_interp(:,4),NoI1,NoI2))
shading interp
set(gca,'Fontsize',6)
ylabel('$k_nl$ (\%)','Interpreter','Latex','Fontsize',7)
xlabel('$F_{\max}$ (\%)','Interpreter','Latex','Fontsize',7)
zlabel(['$a_',num2str(1),'(\xi)$'],'Interpreter','Latex','Fontsize',7)
grid on
box on











%% Transient analysis (random excitation) -- LHS -- Chirp
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
% Number of samples
N = 1000;
% Sampling period (sampling freq. = 100 Hz)
Ts = 0.005;                 % Max. freq ~= 5 Hz
% Time index
Tspan = (0:Ts:(N-1)*Ts);
% Number of variables (knl, std(F))
M = 2;
% Number of experiments
K = 20; 

% ODE45 parameters
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
IC = zeros(3,1);
% Model properties
BW.n = 3;    
BW.gamma = 0;    
BW.A = 1;         

Sys.m = 1;
Sys.kf = 0;
Sys.c = 0;
Sys.ki = 5000;

% Latin Hypercube Sampling
Q = 2*(lhsdesign(K,M,'criterion','maximin','iterations',10)-0.5);
Qtrans(:,1) = Q(:,1)*0.5 + 0.5;                                   % beta
Qtrans(:,2) = 2500*Q(:,2)+2500;                                   % Fmax

[dspl,veloc,force] = deal(zeros(1000,K));
for k = 1:K
    % Model properties
    disp(k)
    BW.B = Qtrans(k,1);
    F = Qtrans(k,2)*chirp(Tspan,0,Tspan(end),5);
    [Tsim,Xsim] = ode45(@(t,x) boucwenSDOF(t,x,Tspan,F,BW, Sys),Tspan,IC,options);
    for t = 1:length(Tspan)
        ACC(t,:) = boucwenSDOF(Tspan(t),Xsim(t,:),Tspan(t),F(t),BW,Sys);
    end
    dspl(:,k) = Xsim(1:end,1);
    veloc(:,k) = Xsim(1:end,2);
    acc(:,k) = ACC(1:end,2);
    force(:,k) = F(1:end);
    resforce(:,k) = Xsim(1:end,3);
end


figure(1)
clr = colormap(jet(K));
subplot(411),hold on
subplot(412),hold on
subplot(413),hold on
subplot(414),hold on
for i = 1:K
    [Pf,F] = pwelch(force(:,i),128,120,256,1/(4*Ts));
    [Py,F] = pwelch(veloc(:,i),128,120,256,1/(4*Ts));
    [Txy,F] = tfestimate(force(:,i),veloc(:,i),128,120,256,1/(4*Ts));
    [Cxy,F] = mscohere(force(:,i),veloc(:,i),128,120,256,1/(4*Ts));
    subplot(411),plot(F,20*log10(abs(Pf)),'color',clr(i,:))
    subplot(412),plot(F,20*log10(abs(Py)),'color',clr(i,:))
    subplot(413),plot(F,20*log10(abs(Txy)),'color',clr(i,:))
    subplot(414),plot(F,Cxy,'color',clr(i,:))
end
save('SDOF_BoucWen_RandomExc_LHS_chirp.mat','resforce','force','dspl','veloc','acc','N','Sys','Ts','Tspan','K','M','Q*')






%% Parameter Surfaces -- [Figures]
% =========================================================================
clear Qsurf
close all;
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

PCBASES = {indx{1}(basisindx(1:nreg),:)};
NLREGRESSORS = regressors(regindx(1:nreg));

options.basis = 'legen';
NoI1 = 200;
NoI2 = 200;
aux = linspace(-1,1,NoI1);
Qsurf(1,:) = kron(ones(1,NoI2),aux);
Qsurf(2,:) = kron(aux,ones(1,NoI2));
[~,an_interp] = PCparam(Qsurf,3,PCBASES{1}([5:9 22 23],:),thetaSIM{nreg}.sim([5:9 22 23],:),options);
% [~,an_interp] = PCparam(Qsurf,3,PCBASES{1}([4 14 15],:),thetaSIM{nreg}.sim([4 14 15],:),options);

figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
% subplot(2,2,1)
surf(linspace(0,2500,NoI2),linspace(0,5,NoI1),reshape(an_interp(:,1),NoI1,NoI2))
shading interp
set(gca,'Fontsize',10)
ylabel('$\beta$','Interpreter','Latex','Fontsize',15)
xlabel('$\sigma_f$','Interpreter','Latex','Fontsize',15)
zlabel(['$\hat{\theta}_{y[t-3]^3}\ (\xi)$'],'Interpreter','Latex','Fontsize',15)
grid on
box on
axis tight
view([-125 25])
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_surf']),close;
     result = eps2xxx([write_dir,'boucwen_surf.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end





%% GA optimization


for k = 1:K
    [sg(k,:), sl(k,:)] = glstat(Y(:,k),0.51,128);
    [bic,waxis] = bicoherx(X(:,k),X(:,k),Y(:,k), 256, [], 128, 50);
    maxbic(k) = max(max(bic));
end
figure(1)
subplot(311),plot(1:K,sg(:,3))
subplot(312),plot(1:K,sl(:,1),'-bo',1:K,sl(:,3),'-rd')
legend({'Estimated','Theory'})
subplot(313),plot(1:K,maxbic,'-bo')
ii = input('Given run with stronger nonlinearities: ');
X = force(1:N,ii);
Y = dspl(1:N,ii); 




% Nonlinear polynomial regressors
P = 3;      q = 1;      maxterms = 3;
% Model orders
na = 4;     nb = 4;     nd = 0;
regressors = polyregNARX([na nd nb],P,q,maxterms,{'abs(Y(t-1,:))'});
options.criterion = 'rss';
options.warnings = 'off';
options.GA.PopulationSize = 100;
options.focus = 'simulation';
rng(10)
[RegGA,indx,RES2,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);


clear *criteria* NLres*
for k = 1:20
    X = force(1:N,k);
    Y = veloc(1:N,k);
    Nreg = length(RegGA);
    options.focus = 'prediction';
    [thetaNL,NLresPE{k},PEcriteria] = narx(Y,X,[na nd nb],RegGA,[],[],[],[],options);
    rssPE(k) = PEcriteria.rss_sss;
%     options.PEinit = 'y';
%     options.focus = 'simulation';
%     options.nlmethod = 'LM';
%     [~,NLresSIM{k},SIMcriteria{k}] = narx(Y,X,[na nd nb],RegGA,[],[],[],thetaNL,options);
end


options.focus = 'prediction';
[thetaNL,NLresPE,PEcriteria] = narx(Y,X,[na nd nb],regressors(regindx(1:23)),[],[],[],[],options);
options.focus = 'simulation';
options.PEinit = 'y';
[thetaNL,NLresSIM,SIMcriteria] = narx(Y,X,[na nd nb],regressors(regindx(1:23)),[],[],[],[],options);

%% PC basis
options.maxsize = 1e7;
options.focus = 'prediction';
options.criterion = 'bic';
options.nlreg = 'n';
options.basis = 'legen';
options.common = 'y';
options.GA.PopulationSize = 100;
options.GA.Display = 'iter';
options.GA.PlotFcns = {@gaplotbestf,@gaplotbestindiv};
[RegGA, IndxGA,thetaij,res,criteria,GAoutput]  = gaPCNARX(dspl(1:N,:),force(1:N,:),Q',[na nd nb],regressors(regindx(1:23)),3,1,options);


clear INDX
INDX{1} = combinations(2,3,1);
options.focus = 'prediction';
options.method = 'wls';
options.wls_iter = 20;
[THij,res,criteria,output] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],IndxGA,RegGA,[],options);
% SE estimation
options.focus = 'simulation';
options.PEinit = 'n';
[THij_sim,res_sim,criteria_sim] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],INDX,regressors(regindx(1:23)),[],options);




%% Bouc-Wen (modes-Poincare) 
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
% Number of samples
N = 2000;
% Sampling period (sampling freq. = 100 Hz)
Ts = 0.005;                 % Max. freq ~= 5 Hz
% Time index
Tspan = (0:Ts:(N-1)*Ts);
% Number of variables (knl, std(F))
M = 2;

% ODE45 parameters
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
IC = zeros(3,1);
% Model properties
BW.n = 3;    
BW.gamma = 0;    
BW.A = 1;         
BW.B = 0;

Sys.m = 1;
Sys.c = 0.001;
Sys.kf = 5;
Sys.ki = 50;
omega = 2*2*pi;          % (rad/s)
Fmax = 200;
beta = [0:0.2:1];

figure(1)
clr = colormap(lines);
subplot(311),hold on 
subplot(312),hold on 
subplot(313),hold on 
for p = 1:length(omega)
    for q = 1:length(Fmax)
        for r = 1:length(beta)
            disp([omega(p),Fmax(q),beta(r)])
            BW.B = beta(r);
            U = Fmax(q)*sin(omega(p)*Tspan);
            [Tsim,Ysim{r}] = ode45(@(t,x) boucwenSDOF(t,x,Tspan,U,BW, Sys),Tspan,IC,options);
            subplot(311),plot(Ysim{r}(:,1),Ysim{r}(:,3),'Color',clr(r,:))
            subplot(312),plot(Tsim,cumsum(Ysim{r}(:,3)),'Color',clr(r,:))
            subplot(313),plot(Tsim,cumsum(Ysim{r}(:,3))./cumsum(U'),'Color',clr(r,:))
            sumY(r) = norm(Ysim{r}(:,3));
        end
    end
end




%% Bouc-Wen model NARX regressors selection through CCF
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('SDOF_BoucWen_RandomExc_LHS.mat')

options.maxsize = 1e7;
options.focus = 'prediction';
N = 1000;

X = force(1:N,:);
Y = veloc(1:N,:);

for k = 1:K
    CCF(:,k,1) = ccf(Y(:,k)-mean(Y(:,k)),(Y(:,k)-mean(Y(:,k))).^2,100,0.8);
    CCF(:,k,2) = ccf(Y(:,k)-mean(Y(:,k)),(Y(:,k)-mean(Y(:,k))).^3,100,0.8);
    CCF(:,k,3) = ccf(Y(:,k)-mean(Y(:,k)),(Y(:,k)-mean(Y(:,k))).^4,100,0.8);
    CCF(:,k,4) = ccf(Y(:,k)-mean(Y(:,k)),(Y(:,k)-mean(Y(:,k))).^5,100,0.8);
    CCF(:,k,5) = ccf(Y(:,k)-mean(Y(:,k)),abs(X(:,k)-mean(X(:,k))),100,0.8);
    CCF(:,k,6) = ccf(Y(:,k)-mean(Y(:,k)),sin(Y(:,k)-mean(Y(:,k))),100,0.8);
    CCF(:,k,7) = ccf(Y(:,k)-mean(Y(:,k)),exp(Y(:,k)-mean(Y(:,k))),100,0.8);
    for j = 1:7
        maxCCF(k,j) = max(CCF(:,k,j));
        sumCCF(k,j) = sum(CCF(:,k,j));
    end
end



%% PC-NARX estimation (full space)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('SDOF_Cubic_RandomExc_LHS.mat')

options.maxsize = 1e7;
options.focus = 'prediction';
N = 1000;

% Signals
X = force(1:N,:);
Y = dspl(1:N,:);

% PE estimation
options.maxsize = 1e8;
options.basis = 'legen';


N = 1000;
options.method = 'ols';
options.wls_iter = 5;

for ord = 2:4
    for pow = 1:5
        for pc = 5
            % Model orders
            na = ord;     nb = ord;     nd = 0;
            P = pow;      q = 1;      maxterms = 1;
            regressors = polyregNARX([na nd nb],P,q,maxterms);
            options.criterion = 'rss';
            indx{1} = combinations(2,pc,1);
            indx{2} = indx{1};
            [THij,res,criteria,output] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],indx,regressors ,[],options);
            [res,criteria] = pcnarxSIM(THij.ols,dspl(1:N,:),force(1:N,:),Q',[na 0 na],indx,regressors,options);
            rss(ord,pow,pc) = criteria.rss_sss;
        end
    end
end
