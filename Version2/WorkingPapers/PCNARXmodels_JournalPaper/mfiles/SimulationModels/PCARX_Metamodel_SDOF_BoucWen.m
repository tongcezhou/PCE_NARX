%% NARX approximation of the SDOF Bouc-Wen model
syms m ki b u_i udot v vdot vddot z zdot Ts v_ip1 v_i v_im1 u_ip1 u_im1

% vdot = (v_i - v_im1)/Ts;                  % (1)
% vddot = (v_ip1 - 2*v_i + v_im1)/(Ts^2);   % (2)
% udot = (u_i - u_im1)/Ts;                  % (3)
% vddot = -ki*zdot/m + udot/m;              % (4)
% z = u_i/ki - vdot*m/ki;                   % (5)
% zdot = v_i - b*(z^3)*(abs(v_i));          % (6)

v_ip1 =  + 2*v_i - v_im1 + (Ts^2)*(-ki*(v_i - b*((u_i/ki - ((v_i - v_im1)/Ts)*m/ki)^3)*(abs(v_i)))/m + ((u_i - u_im1)/Ts)/m);
str = expand(v_ip1);
m = 1; b=1; ki = 10; Ts=0.005;
eval(str)


%% Bouc-Wen (modes-Poincare) 
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')


clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
% Number of samples
N = 3000;
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
BW.B = 1;

Sys.m = 1;
Sys.kf = 0;
Sys.c = 0;
Sys.ki = 50;

omega = 2*2*pi;          % (rad/s)
Fmax = [1 200];
beta = [0 1];

for p = 1:length(omega)
    for q = 1:length(Fmax)
        for r = 1:length(beta)
            disp([omega(p),Fmax(q),beta(r)])
            BW.B = beta(r);
            U = Fmax(q)*sin(omega(p)*Tspan);
            [Tsim,Ysim{q,r}] = ode45(@(t,x) boucwenSDOF(t,x,Tspan,U,BW, Sys),Tspan,IC,options);
            A(p,q,r) = (max(abs(Ysim{q,r}(1001:end,1)))/max(abs(U)))^2;
            Tsin = 1/(omega(p)/(2*pi));
            Tstart = round(Tsin/(4*Ts));
            Tperiod = round(Tsin/Ts);
            Tindx = Tstart:Tperiod:N;
            Poincare{p,q,r}(1:length(Tindx),1) = Ysim{q,r}(Tindx,1);
            Poincare{p,q,r}(1:length(Tindx),2) = Ysim{q,r}(Tindx,2);
        end
    end
end
save('Results/SDOF_BoucWen_Sine_Poincare.mat','omega','Fmax','beta','N','Poincare','Sys','BW','Ts','Tspan','A','Ysim')


%% Bouc-Wen (modes-Poincare)  -- Figures
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';
load('Results/SDOF_BoucWen_Sine_Poincare.mat')


figure(1)
clr = colormap(lines(2));
hold on
figindx = [1 1 2 2];
colindx = [1 2 1 2];
linindx = [3 1 3 1];

k = 0;
for q = 1:length(Fmax)
    for r = 1:length(beta)
        k = k+1;
        U = Fmax(q)*sin(omega(1)*Tspan);
        subplot(2,2,figindx(k)),plot(Ysim{q,r}(2001:end,1),Ysim{q,r}(2001:end,3),...
            'Linewidth',linindx(k),'Color',clr(colindx(k),:))
        hold on
        str{r} = ['$\beta$ = ',num2str(beta(r))];
        subplot(4,2,4+figindx(k)),plot(0.01*(0:999),Ysim{q,r}(2001:end,3),...
            'Linewidth',linindx(k),'Color',clr(colindx(k),:))
        hold on
        axis tight
    end
    subplot(2,2,figindx(k)),legend(str,'Fontsize',9,'Location','SouthEast','Interpreter','Latex')
    set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
    ylabel('Force (N)','Fontsize',11,'Fontname','TimesNewRoman')
    xlabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
    grid on
    box on
    subplot(4,2,4+figindx(k))
    set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
    ylabel('Restoring force (N)','Fontsize',11,'Fontname','TimesNewRoman')
    xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
    grid on
    box on
    if q ==1,   ylim=[-1 1];   else ylim = [-10 10];    end
end
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_sinexci']),close;
     result = eps2xxx([write_dir,'boucwen_sinexci.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


figure(2)
linestyle = {'-o','-+'};
clr = colormap(lines(6));
hold on
for j = 1:length(beta)
    for i = 1:length(Fmax)
        plot3(j*ones(1,length(omega)),omega,20*log10(squeeze(A(:,i,j)))',linestyle{i},'Linewidth',1.5,'Color',clr(i,:))
    end
end
set(gca,'xtick',1:length(beta),'xticklabel',beta)
ylabel('Frequency (Hz)')
xlabel('\beta')
zlabel('|A| (dB)')
legend({'F_{max} = 1','F_{max} = 5'})
view([-25 25])
grid on

figure(3)
plot(Poincare{1,2,2}(:,2),Poincare{1,2,2}(:,1),'.')


%% Transient analysis (sweep sine) -- Fine grid
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')

% Number of samples
N = 2000;
% Sampling period (sampling freq. = 100 Hz)
Ts = 0.005;                 % Max. freq ~= 5 Hz
% Time index
Tspan = (0:Ts:(N-1)*Ts);
% Number of variables (knl, std(F))
M = 2;

% ODE45 parameters
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
IC = zeros(3,1);
% Model properties
BW.n = 3;    
BW.gamma = 0;    
BW.A = 1;         
BW.B = 1;

Sys.m = 1;
Sys.kf = 0;
Sys.c = 0;
Sys.ki = 50;

SweepSine= chirp(Tspan,0,Tspan(2000),5);

% Random input variables
Fmax = 50:50:200;
beta = 0:0.25:1;

% Initialization
indxF = length(Fmax);
indxb = length(beta);
K = indxF*indxb;
[dspl,force,resforce] = deal(zeros(2000,K));

k=0;
for i = 1:indxb
    for j = 1:indxF
        k = k+1;
        % Model properties
        disp([i,j])
        BW.B = beta(i);
        rng(k) 
        F = Fmax(j)*SweepSine;
        [Tsim,Xsim] = ode45(@(t,x) boucwenSDOF(t,x,Tspan,F,BW, Sys),Tspan,IC,options);
        dspl(:,k) = Xsim(:,1);
        veloc(:,k) = Xsim(:,2);
        resforce(:,k) = Xsim(:,3);
        force(:,k) = F(:);
    end
end

N = size(dspl,1);
K = size(dspl,2);
save('Results/SDOF_BoucWen_SweepSine_FineGrid.mat','force','veloc','resforce',...
    'beta','Fmax','indxb','indxF','dspl','N','K','Sys','Ts','Tspan')


%% Transient analysis (random excitation) -- Fine grid
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')

% Number of samples
N = 3000;
% Sampling period (sampling freq. = 100 Hz)
Ts = 0.005;                 % Max. freq ~= 5 Hz
% Time index
Tspan = (0:Ts:(N-1)*Ts);
% Number of variables (knl, std(F))
M = 2;

% ODE45 parameters
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
IC = zeros(3,1);
% Model properties
BW.n = 3;    
BW.gamma = 0;    
BW.A = 1;         
BW.B = 1;

Sys.m = 1;
Sys.kf = 0;
Sys.c = 0;
Sys.ki = 50;

% Filter construction
forder = cheb2ord(0.075, 0.1, 0.01, 60);
[Z,P,G] = cheby2(forder,60,0.1);                       % 0.5 Hz
[numf,denf] = zp2tf(Z,P,G);

% Random input variables
Fmax = 50:50:200;
beta = 0:0.25:1;

% Initialization
indxF = length(Fmax);
indxb = length(beta);
K = indxF*indxb;
[dspl,force,resforce] = deal(zeros(1500,K));

k=0;
for i = 1:indxb
    for j = 1:indxF
        k = k+1;
        % Model properties
        disp([i,j])
        BW.B = beta(i);
        rng(k) 
        F = Fmax(j)*randn(N,1);
        F = filtfilt(numf,denf,F);
        [Tsim,Xsim] = ode45(@(t,x) boucwenSDOF(t,x,Tspan,F,BW, Sys),Tspan,IC,options);
        dspl(:,k) = Xsim(1501:3000,1);
        veloc(:,k) = Xsim(1501:3000,2);
        resforce(:,k) = Xsim(1501:3000,3);
        force(:,k) = F(1501:3000);
    end
end

N = size(dspl,1);
K = size(dspl,2);

save('Results/SDOF_BoucWen_RandomExc_FineGrid.mat','force','resforce',...
    'beta','Fmax','indxb','indxF','dspl','veloc','N','K','Sys','Ts','Tspan')




%% Transient analysis (random excitation) -- Fine grid -- Figures
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
load('Results/SDOF_BoucWen_RandomExc_FineGrid.mat')

figure(1)
clr = colormap(jet(K));
subplot(411),hold on
subplot(412),hold on
subplot(413),hold on
subplot(414),hold on
for i = 1:K
    [Pf,F] = pwelch(force(:,i),128,120,256,1/Ts);
    [Py,F] = pwelch(dspl(:,i),128,120,256,1/Ts);
    [Txy,F] = tfestimate(force(:,i),dspl(:,i),128,120,256,1/Ts);
    [Cxy,F] = mscohere(force(:,i),dspl(:,i),128,120,256,1/Ts);
    subplot(411),plot(F,20*log10(abs(Pf)),'color',clr(i,:))
    subplot(412),plot(F,20*log10(abs(Py)),'color',clr(i,:))
    subplot(413),plot(F,20*log10(abs(Txy)),'color',clr(i,:))
    subplot(414),plot(F,Cxy,'color',clr(i,:))
end


figure(2)
clr = colormap(lines(K));
hold on
for i = 1:indxb
    for j = 1:indxF
        k = (i-1)*indxF + j;
        [Py,F] = pwelch(dspl(:,k),128,120,256,1/Ts);
        plot3(beta(i)*ones(length(F),1),F,20*log10(abs(Py)),'color',clr(j,:))
    end
end
view([80 40])






%% Conventional ARX - NARX identification (prediction - O(h) approximation)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
RandomExc = 0;
if RandomExc == 0
    load('Results/SDOF_BoucWen_RandomExc_FineGrid.mat')
    N = 1000;
else
    load('Results/SDOF_BoucWen_SweepSine_FineGrid.mat')
    N = 1000;
end

options.focus = 'prediction';
% Model orders
na = 2;     nb = 2;     nd = 1;
% Constant term: 'ones(length(t),size(Y,2))';
ARXReg{1} = 'ones(length(t),size(Y,2))';
ARXReg{2} = '(Y(t-1,:).^1)';
ARXReg{3} = '(Y(t-2,:).^1)';
ARXReg{4} = '(X(t-1,:).^1)';
ARXReg{5} = '(X(t-2,:).^1)';

% NARX regressors
NARXReg{1} = 'ones(length(t),size(Y,2))';
NARXReg{2} = '(Y(t-1,:).^1)';
NARXReg{3} = '(Y(t-2,:).^1)';
NARXReg{4} = '(X(t-1,:).^1)';
NARXReg{5} = '(X(t-2,:).^1)';
NARXReg{6} = '(X(t-1,:).^3).*abs(Y(t-1,:))';
NARXReg{7} = '(Y(t-1,:).^1).*(X(t-1,:).^2).*abs(Y(t-1,:))';
NARXReg{8} = '(Y(t-2,:).^1).*(X(t-1,:).^2).*abs(Y(t-1,:))';
NARXReg{9} = '(Y(t-1,:).^2).*(X(t-1,:).^1).*abs(Y(t-1,:))';
NARXReg{10} = '(Y(t-1,:).^1).*(Y(t-2,:).^1).*(X(t-1,:).^1).*abs(Y(t-1,:))';
NARXReg{11} = '(Y(t-2,:).^2).*(X(t-1,:).^1).*abs(Y(t-1,:))';
NARXReg{12} = '(Y(t-1,:).^3).*abs(Y(t-1,:))';
NARXReg{13} = '(Y(t-1,:).^2).*(Y(t-2,:).^1).*abs(Y(t-1,:))';
NARXReg{14} = '(Y(t-1,:).^1).*(Y(t-2,:).^2).*abs(Y(t-1,:))';
NARXReg{15} = '(Y(t-2,:).^3).*abs(Y(t-1,:))';


options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','off',...
    'FinDiffType','Central','FinDiffRelStep',sqrt(eps),...
    'TolFun',1e-3,'TolX',1e-8,'MaxFunEval',1000,'MaxIter',1000);

nlopt = 'n';
disp('Linear    Nonlinear')
disp('===================')
m = 1; n =0;
for ii = 1:K
    X = force(1:N,ii);
    Y = veloc(1:N,ii);
    n = n+1;
    if n > indxF
        n = 1;
        m = m +1;
    end
    % ARX estimation
    options.focus = 'prediction';
    [theta(:,ii),res,criteria] = narx(Y,X,[na nd nb],ARXReg,[],[],[],[],options);
    rss_sss(ii) = criteria.rss_sss;
    % NARX estimation
    options.focus = 'prediction';
    [thetaNL(:,ii),NLres,NLcriteria,regstr] = narx(Y,X,[na nd nb],NARXReg,[],[],[],[],options);
    NLrss_sss(ii) = NLcriteria.rss_sss;
    
    if strcmp(nlopt,'y')
        options.focus = 'simulation';
        % ARX estimation
        [~,res,criteria] = narx(Y,X,[na nd nb],ARXReg,[],[],[],theta(:,ii),options);
        rss_sss_sim(ii) = criteria.rss_sss;
        % NARX estimation
        [~,res,criteria] = narx(Y,X,[na nd nb],NARXReg,[],[],[],thetaNL(:,ii),options);
        NLrss_sss_sim(ii) = criteria.rss_sss;
    else
        % ARX estimation
        Ysim = narxsim(X,Y(1:na),na,ARXReg,theta(:,ii));
        rss_sss_sim(ii) = 100*norm(Y-Ysim)^2/norm(Y)^2;
        % NARX estimation
        NLYsim = narxsim(X,Y(1:na),na,NARXReg,thetaNL(:,ii));
        NLrss_sss_sim(ii) = 100*norm(Y-NLYsim)^2/norm(Y)^2;
    end
    
    s = sprintf('%2.5f \t (%2.5f) \t %2.5f \t (%2.5f) [beta = %2.2f, Fmax = %2.4f]',rss_sss(ii),rss_sss_sim(ii),NLrss_sss(ii),NLrss_sss_sim(ii),beta(m),Fmax(n));
    disp(s)
end

close all
figure(1)
plot(rss_sss,'-o')
hold on
plot(NLrss_sss,'r-d')
legend({'ARX(2,1,1)','NARX(2,1,1)'})
xlabel('Simulation experiment')
ylabel('RSS/SSS (%)')

figure(2)
for k=1:10
    q = 0;
    for i = 1:length(beta)
        for j = 1:length(Fmax)
            q = q+1;
            thetasurf(i,j,k) = thetaNL(k,q);
        end
    end
    subplot(2,5,k),surf(Fmax(1:end),beta,thetasurf(:,1:end,k))
    hold on
    ylabel('beta' )
    xlabel('std(F) (N)' )
    shading interp
end


%% Transient analysis (random excitation) -- LHS
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
% Number of samples
N = 3000;
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
Sys.ki = 50;

% Filter construction
forder = cheb2ord(0.15, 0.2, 0.1, 60);
[Z,P,G] = cheby2(forder,60,0.2);                       % 0.5 Hz
[numf,denf] = zp2tf(Z,P,G);

% Latin Hypercube Sampling
rng(500)
Q = 2*(lhsdesign(K,M,'criterion','maximin','iterations',10)-0.5);
Qtrans(:,1) = Q(:,1)*0.5 + 0.5;                                   % beta
Qtrans(:,2) = 50*Q(:,2)+50;                                   % Fmax

[dspl,veloc,force] = deal(zeros(1500,K));
for k = 1:K
    % Model properties
    disp(k)
    BW.B = Qtrans(k,1);
%     stream = RandStream('mt19937ar','seed',1000+k);
%     RandStream.setDefaultStream(stream);
    rng(k+5000)
    F = Qtrans(k,2)*randn(N+1000,1);
    F = filtfilt(numf,denf,F);
    F = F(501:end-500);
    [Tsim,Xsim] = ode45(@(t,x) boucwenSDOF(t,x,Tspan,F,BW, Sys),Tspan,IC,options);
    for t = 1:length(Tspan)
        ACC(t,:) = boucwenSDOF(Tspan(t),Xsim(t,:),Tspan(t),F(t),BW,Sys);
    end
    dspl(:,k) = Xsim(1501:3000,1);
    veloc(:,k) = Xsim(1501:3000,2);
    acc(:,k) = ACC(1501:3000,2);
    force(:,k) = F(1501:3000);
    resforce(:,k) = Xsim(1501:3000,3);
end


figure(1)
clr = colormap(jet(K));
subplot(411),hold on
subplot(412),hold on
subplot(413),hold on
subplot(414),hold on
for i = 1:K
    [Pf,F] = pwelch(force(:,i),128,120,256,1/Ts);
    [Py,F] = pwelch(veloc(:,i),128,120,256,1/Ts);
    [Txy,F] = tfestimate(force(:,i),veloc(:,i),128,120,256,1/Ts);
    [Cxy,F] = mscohere(force(:,i),veloc(:,i),128,120,256,1/Ts);
    subplot(411),plot(F,20*log10(abs(Pf)),'color',clr(i,:))
    subplot(412),plot(F,20*log10(abs(Py)),'color',clr(i,:))
    subplot(413),plot(F,20*log10(abs(Txy)),'color',clr(i,:))
    subplot(414),plot(F,Cxy,'color',clr(i,:))
end
save('Results/SDOF_BoucWen_RandomExc_LHS.mat','resforce','force','dspl','veloc','acc','N','Sys','Ts','Tspan','K','M','Q*')


%% Transient analysis (random excitation) -- LHS -- validation
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
% Number of samples
N = 3000;
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
Sys.ki = 50;

% Filter construction
forder = cheb2ord(0.15, 0.2, 0.1, 60);
[Z,P,G] = cheby2(forder,60,0.2);                       % 0.5 Hz
[numf,denf] = zp2tf(Z,P,G);

% Latin Hypercube Sampling
rng(400)
Q = 2*(lhsdesign(K,M,'criterion','maximin','iterations',10)-0.5);
Qtrans(:,1) = Q(:,1)*0.5 + 0.5;                                   % beta
Qtrans(:,2) = 50*Q(:,2)+50;                                   % Fmax

[dspl,veloc,force] = deal(zeros(1500,K));
for k = 1:K
    % Model properties
    disp(k)
    BW.B = Qtrans(k,1);
%     stream = RandStream('mt19937ar','seed',1000+k);
%     RandStream.setDefaultStream(stream);
    rng(k+10000)
    F = Qtrans(k,2)*randn(N+1000,1);
    F = filtfilt(numf,denf,F);
    F = F(501:end-500);
    [Tsim,Xsim] = ode45(@(t,x) boucwenSDOF(t,x,Tspan,F,BW, Sys),Tspan,IC,options);
    for t = 1:length(Tspan)
        ACC(t,:) = boucwenSDOF(Tspan(t),Xsim(t,:),Tspan(t),F(t),BW,Sys);
    end
    dspl(:,k) = Xsim(1501:3000,1);
    veloc(:,k) = Xsim(1501:3000,2);
    acc(:,k) = ACC(1501:3000,2);
    force(:,k) = F(1501:3000);
    resforce(:,k) = Xsim(1501:3000,3);
end


figure(1)
clr = colormap(jet(K));
subplot(411),hold on
subplot(412),hold on
subplot(413),hold on
subplot(414),hold on
for i = 1:K
    [Pf,F] = pwelch(force(:,i),128,120,256,1/Ts);
    [Py,F] = pwelch(veloc(:,i),128,120,256,1/Ts);
    [Txy,F] = tfestimate(force(:,i),veloc(:,i),128,120,256,1/Ts);
    [Cxy,F] = mscohere(force(:,i),veloc(:,i),128,120,256,1/Ts);
    subplot(411),plot(F,20*log10(abs(Pf)),'color',clr(i,:))
    subplot(412),plot(F,20*log10(abs(Py)),'color',clr(i,:))
    subplot(413),plot(F,20*log10(abs(Txy)),'color',clr(i,:))
    subplot(414),plot(F,Cxy,'color',clr(i,:))
end
save('Results/SDOF_BoucWen_RandomExc_LHS_validation.mat','resforce','force','dspl','veloc','acc','N','Sys','Ts','Tspan','K','M','Q*')



%% Transient analysis (random excitation) -- Validation set (equispaced grid)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
% Number of samples
N = 3000;
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
Sys.ki = 50;

% Filter construction
forder = cheb2ord(0.15, 0.2, 0.1, 60);
[Z,P,G] = cheby2(forder,60,0.2);                       % 0.5 Hz
[numf,denf] = zp2tf(Z,P,G);

% Equispaced grid
Q(:,1) = kron([-1 -0.5 0 0.5 1],ones(1,4));     % beta
Q(:,2) = kron(ones(1,5),[-0.5 0 0.5 1]);        % Fmax
Qtrans(:,1) = Q(:,1)*0.5 + 0.5;                 % beta
Qtrans(:,2) = 50*Q(:,2)+50;                     % Fmax

[dspl,veloc,force] = deal(zeros(1500,K));
for k = 1:K
    % Model properties
    disp(k)
    BW.B = Qtrans(k,1);
%     stream = RandStream('mt19937ar','seed',1000+k);
%     RandStream.setDefaultStream(stream);
    rng(k+6000)
    F = Qtrans(k,2)*randn(N+1000,1);
    F = filtfilt(numf,denf,F);
    F = F(501:end-500);
    [Tsim,Xsim] = ode45(@(t,x) boucwenSDOF(t,x,Tspan,F,BW, Sys),Tspan,IC,options);
    for t = 1:length(Tspan)
        ACC(t,:) = boucwenSDOF(Tspan(t),Xsim(t,:),Tspan(t),F(t),BW,Sys);
    end
    dspl(:,k) = Xsim(1501:3000,1);
    veloc(:,k) = Xsim(1501:3000,2);
    acc(:,k) = ACC(1501:3000,2);
    force(:,k) = F(1501:3000);
    resforce(:,k) = Xsim(1501:3000,3);
end


figure(1)
clr = colormap(jet(K));
subplot(411),hold on
subplot(412),hold on
subplot(413),hold on
subplot(414),hold on
for i = 1:K
    [Pf,F] = pwelch(force(:,i),128,120,256,1/Ts);
    [Py,F] = pwelch(veloc(:,i),128,120,256,1/Ts);
    [Txy,F] = tfestimate(force(:,i),veloc(:,i),128,120,256,1/Ts);
    [Cxy,F] = mscohere(force(:,i),veloc(:,i),128,120,256,1/Ts);
    subplot(411),plot(F,20*log10(abs(Pf)),'color',clr(i,:))
    subplot(412),plot(F,20*log10(abs(Py)),'color',clr(i,:))
    subplot(413),plot(F,20*log10(abs(Txy)),'color',clr(i,:))
    subplot(414),plot(F,Cxy,'color',clr(i,:))
end
save('Results/SDOF_BoucWen_RandomExc_validation.mat','resforce','force','dspl','veloc','acc','N','Sys','Ts','Tspan','K','M','Q*')




%% NARX regressors selection (GA - rss)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
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
P = 3;      maxterms = 3; 
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

save('Results/SDOF_BoucWen_RandomExc_GA_rss.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')


%% NARX regressors selection (GA - mnse)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')

%Number of samples
N = 1000;
% Simulation experiment
SimExp = 8;

% Load data
X = force(1:N,SimExp);
Y = veloc(1:N,SimExp);

% Model orders 
na = 2;     nb = 2;     nd = 0;

% Polynomial regressors
P = 3;      maxterms = 3; 
regressors = polyregNARX([na nd nb],P,1,maxterms,{'abs(Y(t-1,:))'});

% Optimization options
options.warnings = 'off';
options.GA.PopulationSize = 200 ;
options.maxsize = 1e8;
   
rng(100);
options.focus = 'simulation';
options.criterion = 'mnse';
options.parfor = 'n';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,[na nd nb],regressors,options);

options.focus = 'prediction';
options.Nr = length(GAreg);
[SEstructure,IndxRem] = SEreduction(Y,X,[na nd nb],GAreg,options);

save('Results/SDOF_BoucWen_RandomExc_GA_mnse.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')


%% Local NARX model estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_GA_rss.mat')
%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = veloc(1:N,:);

% Model orders 
na = 2;     nb = 2;     nd = 0;
% Regressors
SelReg = SEstructure.removed(end-12:end);

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
    mnseSIM(k) = SIMcriteria{k}.mnse;
end

save('Results/SDOF_BoucWen_RandomExc_LocalNARX_lsqnonlin.mat','*criteria','NLres*','theta*','SelReg','orders')


%% Local NARX model estimation (fminsearch)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_GA_rss.mat')
%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = veloc(1:N,:);

% Model orders 
na = 2;     nb = 2;     nd = 0;
% Regressors
SelReg = SEstructure.removed(end-13:end);

orders = [na nd nb];
clear *criteria* *res* theta*
for k=1:K
    disp(k)
    options.focus = 'prediction';
    options.method = 'ols';
    [thetaPE(:,k),NLresPE{k},PEcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    options.focus = 'simulation';
    options.PEinit = 'y';
    options.nlmethod = 'FMINSEARCH';
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
save('Results/SDOF_BoucWen_RandomExc_LocalNARX_fminsearch.mat','*criteria','NLres*','theta*','SelReg','orders')



%% PC-NARX estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_LocalNARX_lsqnonlin.mat')
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
     1     1;    
     0     2;
     0     3;
     1     2];
     
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
options.PEinit = 'y';
options.nlmethod = 'LM';

% Model orders 
na = 2;     nb = 2;     nd = 0;
[TH,res,criteria,output] = pcnarx(Y,X,Q,[na nd nb],INDX,SelReg,[],options);
% [TH,res,criteria,output] = pcnarx(Y,X,Q,[na nd nb],INDX,SelReg,THij,options);
save('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin.mat','criteria','res','TH','output','orders','SelReg','INDX','THij','options')


%% PC-NARX estimation (fminsearch)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_LocalNARX_fminsearch.mat')
%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = veloc(1:N,:);

options.basis = 'legen';

% Basis indx
INDX{1} = [ 0   0; 
    1     0;     
    0     1;
    1     1;
    2     0;
    0     2];
 
% Total number of basis functions
B = size(INDX{1},1);
% Total polynomial degree
P = 3;
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
options.focus = 'simulation';
options.PEinit = 'n';
options.nlmethod = 'FMINSEARCH';

% Model orders 
na = 2;     nb = 2;     nd = 0;
[TH,res,criteria,output] = pcnarx(Y,X,Q,[na nd nb],INDX,SelReg,THij,options);
save('Results/SDOF_BoucWen_RandomExc_PCNARX_fminsearch.mat','criteria','res','TH','output','orders','SelReg','INDX','THij','options')



%% PC-NARX validation
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end

load('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin.mat')
load('Results/SDOF_BoucWen_RandomExc_LHS_validation.mat')
% load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
clear res
noparam = length(INDX{1})*length(SelReg);
for ii = 1:length(Q)
    [~,an(ii,:)] = PCparam(Q(ii,:)',4,INDX{1},reshape(TH.sim,noparam,1),options);
    X = force(1:1000,ii);
    Y = veloc(1:1000,ii);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    res(:,ii) = Y(3:end) - Ysim(3:end);
    rss(1,ii) = 100*norm(res(:,ii))^2/norm(Y(3:end))^2;
    mnse(1,ii) = mean(abs(res(:,ii))./(1+abs(Y(3:end))));
end

load('Results/SDOF_BoucWen_RandomExc_PCNARX_fminsearch.mat')
clear res an
noparam = length(INDX{1})*length(SelReg);
for ii = 1:length(Q)
    [~,an(ii,:)] = PCparam(Q(ii,:)',4,INDX{1},reshape(TH.sim,noparam,1),options);
    X = force(1:1000,ii);
    Y = veloc(1:1000,ii);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    res(:,ii) = Y(3:end) - Ysim(3:end);
    rss(2,ii) = 100*norm(res(:,ii))^2/norm(Y(3:end))^2;
    mnse(2,ii) = mean(abs(res(:,ii))./(1+abs(Y(3:end))));
end

figure(1)
subplot(211),plot(1:20,rss)
xlabel('Experiment number')
ylabel('RSS/SSS %')
subplot(212),plot(1:20,mnse)
xlabel('Experiment number')
ylabel('MNSE')


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
str = {['$\beta = ',num2str(roundn(Qtrans(minexp,1),-2)),', \sigma_f = ',num2str(round(Qtrans(minexp,2))),' (N)$']};
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
str = {['$\beta = ',num2str(roundn(Qtrans(maxexp,1),-2)),', \sigma_f = ',num2str(round(Qtrans(maxexp,2))),' (N)$']};
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
load('Results/SDOF_BoucWen_RandomExc_GA_rss.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

NoReg = length(GAreg);
close all

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(2,1,1),
plot(1:NoReg,SEstructure.values(1:NoReg),'-o','Markersize',4,'Linewidth',0.8)
hold on 
plot([1 NoReg],SEstructure.initvalue*[1 1],'--r','Markersize',4,'Linewidth',0.8)
set(gca,'yscale','linear','Fontsize',7,'xtick',[1:5:NoReg])
axis([0.9 NoReg+0.1 0 100])
grid on
xlabel('Regressors dropped','Fontangle','normal','Fontsize',9,'Horizontalalignment','center')
ylabel('NSSE (%)','Fontangle','normal','Fontsize',9)
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_MSS_StageA']),close;
     result = eps2xxx([write_dir,'boucwen_MSS_StageA.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end



%% NARX regressors selection -- Figure -- zoomed
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results/SDOF_BoucWen_RandomExc_GA_rss.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

NoReg = length(GAreg);
close all

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(2,5,1:2),
plot(1:NoReg,SEstructure.values(1:NoReg),'-o','Markersize',3,'Linewidth',1)
hold on 
plot([1 NoReg],SEstructure.initvalue*[1 1],'--r','Markersize',3,'Linewidth',1)
set(gca,'yscale','linear','Fontsize',7,'xtick',[0:10:NoReg])
xlim([0 NoReg+1])
grid on 
xlabel('Regressors dropped','Fontangle','normal','Fontsize',9,'Horizontalalignment','center')
ylabel('NSSE (%)','Fontangle','normal','Fontsize',9)


subplot(2,5,3:5),
plot(1:NoReg,SEstructure.values(1:NoReg),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 NoReg],SEstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
set(gca,'yscale','linear','Fontsize',7,'xtick',NoReg-20:NoReg,'xticklabel',[],'yticklabel',[])
axis([NoReg-19.5 NoReg 0 130])
grid on 
% text(NoReg/2+0.5,-55,'Regressors dropped','Fontangle','normal','Fontsize',9,'Horizontalalignment','center')
% ylabel('NSSE (%)','Fontangle','normal','Fontsize',9)
for ii = NoReg-19:NoReg
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
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_MSS_StageA']),close;
     result = eps2xxx([write_dir,'boucwen_MSS_StageA.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
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
            load('Results/SDOF_BoucWen_RandomExc_LocalNARX_lsqnonlin.mat')
        case 2
            load('Results/SDOF_BoucWen_RandomExc_LocalNARX_fminunc.mat')
        case 3
            load('Results/SDOF_BoucWen_RandomExc_LocalNARX_fminsearch.mat')
    end
    
    for k=1:K
        rssSIM(j,k) = SIMcriteria{k}.rss_sss;
        mnseSIM(j,k) = SIMcriteria{k}.mnse;
        exitflag(j,k) = SIMcriteria{k}.SIMexitflag;
    end
end

subplot(211),plot(1:K,rssSIM)
grid on 
legend({'lsqnonlin','fminunc','fminsearch'})
subplot(212),plot(1:K,mnseSIM)
grid on 
legend({'lsqnonlin','fminunc','fminsearch'})



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
load('Results/SDOF_Boucwen_RandomExc_LocalNARX_lsqnonlin.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

clear basis options
options.basis = 'legen';
options.criterion = 'ar2';
Pmax = 3;

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
set(gca,'yscale','linear','Fontsize',7,'xtick',[1:NoB],'xticklabel',[])
for ii = 1:NoB
    txtstr = ['$[',num2str(BASISstructure.removed(ii,1)),',',num2str(BASISstructure.removed(ii,2)),']$'];
    text(ii,-0.2,txtstr,'Rotation',25,'Fontsize',7,'Horizontalalignment','center','Interpreter','Latex','Fontsize',10)
end
axis([1 NoB -0.1 0.9])
grid on
text(NoB/2+0.5,-0.4,'Basis functions dropped','Fontangle','normal','Fontsize',9,'Horizontalalignment','center')
% ylabel('mean(R^2_{adj})','Fontangle','normal','Fontsize',9)
ylabel('$\overline{R^2_{\mbox{adj}}}$','Fontangle','normal','Fontsize',11,'Interpreter','Latex')
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'boucwen_MSS_StageB']),close;    
    if ispc
        result = eps2xxx([write_dir,'boucwen_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'boucwen_MSS_StageB.eps'],{'pdf'});
    end
end



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
load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_LocalNARX_lsqnonlin.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

clear basis options
options.basis = 'legen';
options.criterion = 'ar2';
Pmax = 5;
for P = 0:Pmax
    % Basis index
    options.Nb = 0;
    [BASISstructure,INDnew] = PCEreduction(Q',P,thetaSIM,options);
    R2o(P+1)= BASISstructure.initvalue;
end


options.Nb = 5; Nb = 5;
% Basis index
P = 2;
[BASISstructure,INDnew] = PCEreduction(Q',P,thetaSIM,options);

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(221),plot(0:Pmax,R2o,'-o','Markersize',6,'Linewidth',1.2)
grid on
ylabel('Mean normalized error','Fontsize',11)
xlabel(sprintf('%s\n%s','Total PC basis degree','(complete subspace)'),'Fontsize',10)
set(gca,'yscale','linear','Fontsize',9,'xtick',[0:4])
subplot(222),plot(1:Nb,BASISstructure.values(1:Nb),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 Nb],BASISstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
ylabel('Mean normalized error','Fontsize',11)
set(gca,'yscale','linear','Fontsize',9,'xtick',[1:Nb],'xticklabel',[])
for ii = 1:Nb 
    text(ii,0.25,['$[',num2str(BASISstructure.removed(ii,1)),',',num2str(BASISstructure.removed(ii,2)),']$'],'Rotation',45,...
        'Horizontalalignment','center','Interpreter','Latex','Fontsize',9)
end
xlim([0.99 Nb+.01])
grid on
text(3,-0.25,sprintf('%s\n%s','PC bases dropped','(multivariable indeces)'),'Fontsize',10,'HorizontalAlignment','center')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_MSS_StageB']),close;
     if ispc
        result = eps2xxx([write_dir,'boucwen_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
        result = eps2xxx([write_dir,'boucwen_MSS_StageB.eps'],{'pdf'});
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
print_to_eps = 'n';

K = 20;
load('Results/SDOF_BoucWen_RandomExc_LocalNARX_lsqnonlin.mat')
for k=1:K
    rss(1,k) = SIMcriteria{k}.rss_sss;
    mnse(1,k) = SIMcriteria{k}.mnse;
    R2(1,k) = 1 - SIMcriteria{k}.rss_sss/100;
    exitflag(1,k) = SIMcriteria{k}.SIMexitflag;
end

load('Results/SDOF_BoucWen_RandomExc_LHS.mat')
load('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin.mat')
clear res
options.basis = 'legen';
noparam = size(INDX{1},1)*length(SelReg);
for ii = 1:length(Q)
    [~,an(ii,:)] = PCparam(Q(ii,:)',3,INDX{1},reshape(TH.sim,noparam,1),options);
    X = force(1:1000,ii);
    Y = veloc(1:1000,ii);
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
axis([0.9 20.1 0.985 1])
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'boucwen_local_vs_global']),close;
    if ispc
        result = eps2xxx([write_dir,'boucwen_local_vs_global.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'boucwen_local_vs_global.eps'],{'pdf'});
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
load('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin.mat')
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
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

K = 20;
time = (0:999)*0.005;
load('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin.mat')
load('Results/SDOF_BoucWen_RandomExc_validation.mat')
clear res
options.basis = 'legen';
Pmax = 4;
SimExp1 = 10;
SimExp2 = 18;
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
clr = colormap(lines(2));
hold on
subplot(2,1,1),plot(time,Y(:,1),'-b','Linewidth',1)
hold on
plot(time,Ysim(:,1),'--r','Linewidth',1)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$\beta = ',num2str(roundn(Qtrans(SimExp1,1),-2)),', \sigma_x = ',num2str(round(Qtrans(SimExp1,2))),'\ (N)$']};
title(strT,'Interpreter','Latex')
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
title(strT,'Interpreter','Latex')
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
load('Results/SDOF_BoucWen_RandomExc_PCNARX_lsqnonlin.mat')
load('Results/SDOF_BoucWen_RandomExc_LocalNARX_lsqnonlin.mat')
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
surf(linspace(0,200,NoI2),linspace(0,1,NoI1),reshape(an_interp(:,11),NoI1,NoI2))
shading interp
% hold on
% plot3(Qtrans(:,2),Qtrans(:,1),thetaSIM(2,:)','o')
set(gca,'Fontsize',10)
ylabel('$\beta$','Interpreter','Latex','Fontsize',15)
xlabel('$\sigma_{x}$','Interpreter','Latex','Fontsize',15)
zlabel(['$\hat{\theta}_{x[t]}\ (\xi)$'],'Interpreter','Latex','Fontsize',15)
grid on
box on
% zlim([-0.07 0.01])
view([125 25])
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
surf(linspace(0,200,NoI2),linspace(0,1,NoI1),reshape(an_interp(:,6),NoI1,NoI2))
shading interp
set(gca,'Fontsize',10)
ylabel('$\beta$','Interpreter','Latex','Fontsize',15)
xlabel('$\sigma_x$','Interpreter','Latex','Fontsize',15)
zlabel(['$\hat{\theta}_{y[t-1]^2 \cdot y[t-2] \cdot | y[t-1]|}\qquad \ (\xi)$'],'Interpreter','Latex','Fontsize',15)
grid on
box on
% zlim([1.99 2.01])
view([125 25])
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'boucwen_surf_b']),close;
     if ispc
        result = eps2xxx([write_dir,'boucwen_surf_b.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
         result = eps2xxx([write_dir,'boucwen_surf_b.eps'],{'pdf'});
     end
end


%% Abs function approximation
x = -10:0.1:10;
c2 = [ones(size(x)); x.^2]'\abs(x)';
c4 = [ones(size(x)); x.^2 ; x.^4]'\abs(x)';
c6 = [ones(size(x)); x.^2 ; x.^4 ; x.^6]'\abs(x)';
plot(x,abs(x),'-b',x,c2(1)+c2(2)*x.^2,'-r',...
    x,c4(1)+c4(2)*x.^2+c4(3)*x.^4,'-g',x,c6(1)+c6(2)*x.^2+c6(3)*x.^4+c6(4)*x.^6,'-m')
legend({'abs(x)','c0 + c2*x^2','c0 + c2*x^2 + c4*x^4','c0 + c2*x^2 + c4*x^4 + c6*x^6'})
title('Absolute function approximation')
xlabel('x')
ylabel('f(x)')