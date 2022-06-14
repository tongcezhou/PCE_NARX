%% Duffing oscillator (modes-Poincare) 
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')

% Number of samples
N = 5000;

% Sampling period (sampling freq. = 100 Hz)
Ts = 0.005;                 % Max. freq ~= 5 Hz

% Time index
Tspan = (0:Ts:(N-1)*Ts);

% Number of variables (knl, std(F))
M = 2;

% ODE45 parameters
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
IC = zeros(2,1);

% Model properties
Sys.m = 1;
Sys.k = 5000;
Sys.c = 10;
% Limiting cases
omega = 5*2*pi;          % (rad/s)
Fmax = [1 5000];
knl = [0 5000];

for p = 1:length(omega)
    for q = 1:length(Fmax)
        for r = 1:length(knl)
            disp([omega(p),Fmax(q),knl(r)])
            Sys.knl = knl(r);
            U = Fmax(q)*sin(omega(p)*Tspan);
            [Tsim,Ysim{q,r}] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,U,Sys),Tspan,IC,options);
            A(p,q,r) = (max(abs(Ysim{q,r}(1:end,1)))/max(abs(U)))^2;
            Tsin = 1/(omega(p)/(2*pi));
            Tstart = round(Tsin/(4*Ts));
            Tperiod = round(Tsin/Ts);
            Tindx = Tstart:Tperiod:N;
            Poincare{p,q,r}(1:length(Tindx),1) = Ysim{q,r}(Tindx,1);
            Poincare{p,q,r}(1:length(Tindx),2) = Ysim{q,r}(Tindx,2);
        end
    end
end
save('Results\SDOF_Cubic_Poincare.mat','A','omega','Fmax','knl','N','Poincare','Sys','Ts','Tspan','Ysim')




%% Transient analysis (random excitation) -- Fine grid
% =========================================================================
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
IC = zeros(2,1);
% Model properties
Sys.m = 1;
Sys.k = 5000;
Sys.c = 10;

% Filter construction
forder = cheb2ord(0.15, 0.2, 0.1, 60);
[Z,P,G] = cheby2(forder,60,0.2);                       % 0.5 Hz
[numf,denf] = zp2tf(Z,P,G);

% Random input variables
Fmax = 1000:1000:5000;
knl = 0:1000:5000;

% Initialization
indxF = length(Fmax);
indxk = length(knl);
K = indxF*indxk;
[dspl,force] = deal(zeros(1500,K));

k=0;
for i = 1:indxk
    for j = 1:indxF
        k = k+1;
        % Model properties
        disp([i,j])
        Sys.knl = knl(i);     
        rng(k) 
        F = Fmax(j)*randn(N,1);
        F = filtfilt(numf,denf,F);
        [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F,Sys),Tspan,IC,options); 
        dspl(:,k) = Xsim(501:2000,1);
        force(:,k) = F(501:2000);
    end
end

N = size(dspl,1);
K = size(dspl,2);

save('Results\SDOF_Cubic_RandomExc_FineGrid.mat','force','K','knl','Fmax','indxk','indxF','dspl','N','Sys','Ts','Tspan')



%% Transient analysis (sine sweep) -- Fine grid
% =========================================================================
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
IC = zeros(2,1);
% Model properties
Sys.m = 1;
Sys.k = 5000;
Sys.c = 10;
SweepSine = chirp(Tspan,0,Tspan(2000),5);

% Random input variables
Fmax = 1000:1000:5000;
knl = 0:1000:5000;

% Initialization
indxF = length(Fmax);
indxk = length(knl);
K = indxF*indxk;
[dspl,force] = deal(zeros(2000,K));

k=0;
for i = 1:indxk
    for j = 1:indxF
        k = k+1;
        % Model properties
        disp([i,j])
        Sys.knl = knl(i);     
        F = Fmax(j)*SweepSine;
        [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F,Sys),Tspan,IC,options); 
        dspl(:,k) = Xsim(1:2000,1);
        velo(:,k) = Xsim(1:2000,2);
        force(:,k) = F(1:2000);
        resforce(:,k) = Sys.k*dspl(:,k)+Sys.knl*(dspl(:,k).^3);
    end
end

N = size(dspl,1);
K = size(dspl,2);
save('Results\SDOF_Cubic_SweepSine_FineGrid.mat','force','velo','resforce','K','knl','Fmax','indxk','indxF','dspl','N','Sys','Ts','Tspan')




%% Conventional ARX - NARX identification (prediction - O(h) approximation)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')

RandomExc = 0;
if RandomExc == 1
    load('Results\SDOF_Cubic_RandomExc_FineGrid.mat')
    N = 1000;
else
    load('Results\SDOF_Cubic_SweepSine_FineGrid.mat')
    N = 1000;
end

options.focus = 'prediction';
   
% Model orders
na = 2;     nb = 1;     nd = 1;
% Constant term: 'ones(length(t),size(Y,2))';
ARXReg{1} = 'ones(length(t),size(Y,2))';
ARXReg{2} = '(Y(t-1,:).^1)';
ARXReg{3} = '(Y(t-2,:).^1)';
ARXReg{4} = '(X(t-1,:).^1)';

NARXReg{1} = 'ones(length(t),size(Y,2))';
NARXReg{2} = '(Y(t-1,:).^1)';
NARXReg{3} = '(Y(t-2,:).^1)';
NARXReg{4} = '(X(t-1,:).^1)';
NARXReg{5} = '(Y(t-1,:).^3)';

disp('Linear    Nonlinear')
disp('===================')
for ii = 1:K
    X = force(1:N,ii);
    X = X-mean(X);
    Y = dspl(1:N,ii);
    Y = Y-mean(Y);
    % ARX estimation
    [theta(:,ii),res,criteria] = narx(Y,X,[na nd nb],ARXReg,[],[],[],[],options);
    rss_sss(ii) = criteria.rss_sss;
    % NARX estimation
    [thetaNL(:,ii),NLres,NLcriteria,regstr] = narx(Y,X,[na nd nb],NARXReg,[],[],[],[],options);
    NLrss_sss(ii) = NLcriteria.rss_sss;
    s = sprintf('%2.5f \t %2.5f',rss_sss(ii),NLrss_sss(ii));
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

% Theoretical curves
theta0 = ones(length(knl),length(Fmax))*(Sys.m + Sys.c)/(Sys.m/Ts^2 + Sys.c/(2*Ts));
theta1 = ones(length(knl),length(Fmax))*(2*Sys.m/Ts^2 - Sys.k)/(Sys.m/Ts^2 + Sys.c/(2*Ts));
theta2 = ones(length(knl),length(Fmax))*(-Sys.m/Ts^2 + Sys.c/(2*Ts))/(Sys.m/Ts^2 + Sys.c/(2*Ts));
theta3 = ones(length(knl),length(Fmax))*(1/(Sys.m/Ts^2+Sys.c/(2*Ts)));
theta4 = -kron(knl'/(Sys.m/Ts^2+Sys.c/(2*Ts)),ones(1,length(Fmax)));


figure(2)
for k=1:size(thetaNL,1)
    q = 0;
    for i = 1:length(knl)
        for j = 1:length(Fmax)
            q = q+1;
            thetasurf(i,j,k) = thetaNL(k,q);
        end
    end
    subplot(2,3,k),surf(Fmax(1:end),knl,thetasurf(:,1:end,k))
    hold on
    %if k>1
        surf(Fmax(1:end),knl,eval(['theta',num2str(k-1)]))
    %end
    ylabel('k_{nl} (N/m)' )
    xlabel('std(F) (N)' )
    shading interp
end


%% Conventional ARX - NARX identification (prediction - O(h^2) approximation)
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
RandomExc = 0;
if RandomExc == 1
    load('Results\SDOF_Cubic_RandomExc_FineGrid.mat')
    N = 1000;
else
    load('Results\SDOF_Cubic_SweepSine_FineGrid.mat')
    N = 2000;
end

options.focus = 'prediction';
   
% Model orders
na = 4;     nb = 2;     nd = 2;
ARXReg{1} = 'ones(length(t),size(Y,2))';
ARXReg{2} = '(Y(t-1,:).^1)';
ARXReg{3} = '(Y(t-2,:).^1)';
ARXReg{4} = '(Y(t-3,:).^1)';
ARXReg{5} = '(Y(t-4,:).^1)';
ARXReg{6} = '(X(t-2,:).^1)';

NARXReg{1} = 'ones(length(t),size(Y,2))';
NARXReg{2} = '(Y(t-1,:).^1)';
NARXReg{3} = '(Y(t-2,:).^1)';
NARXReg{4} = '(Y(t-3,:).^1)';
NARXReg{5} = '(Y(t-4,:).^1)';
NARXReg{6} = '(X(t-2,:).^1)';
NARXReg{7} = '(Y(t-2,:).^3)';

disp('Linear    Nonlinear')
disp('===================')
for ii = 1:K
    X = force(1:N,ii);
    X = X-mean(X);
    Y = dspl(1:N,ii);
    Y = Y-mean(Y);
    % ARX estimation
    [theta(:,ii),res,criteria] = narx(Y,X,[na nd nb],ARXReg,[],[],[],[],options);
    rss_sss(ii) = criteria.rss_sss;
    % NARX estimation
    [thetaNL(:,ii),NLres,NLcriteria,regstr] = narx(Y,X,[na nd nb],NARXReg,[],[],[],[],options);
    NLrss_sss(ii) = NLcriteria.rss_sss;
    s = sprintf('%2.5f \t %2.5f',rss_sss(ii),NLrss_sss(ii));
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

% Theoretical curves
Denom = -Sys.m/(12*Ts^2)- Sys.c/(12*Ts);
theta1 = ones(length(knl),length(Fmax))*(-16*Sys.m/(12*Ts^2) - 8*Sys.c/(12*Ts))/Denom;
theta2 = ones(length(knl),length(Fmax))*(30*Sys.m/(12*Ts^2) - Sys.k)/Denom;
theta3 = ones(length(knl),length(Fmax))*(-16*Sys.m/(12*Ts^2) + 8*Sys.c/(12*Ts))/Denom;
theta4 = ones(length(knl),length(Fmax))*(Sys.m/(12*Ts^2) - Sys.c/(12*Ts))/Denom;
theta5 = ones(length(knl),length(Fmax))/Denom;
theta6 = -kron(knl'/Denom,ones(1,length(Fmax)));

figure(2)
for k=1:size(thetaNL,1)
    q = 0;
    for i = 1:length(knl)
        for j = 1:length(Fmax)
            q = q+1;
            thetasurf(i,j,k) = thetaNL(k,q);
        end
    end
    subplot(2,4,k),surf(Fmax(1:end),knl,thetasurf(:,1:end,k))
    hold on
    if k>1
        surf(Fmax(1:end),knl,eval(['theta',num2str(k-1)]))
    end
    ylabel('k_{nl} (N/m)' )
    xlabel('std(F) (N)' )
    shading interp
end



%% Conventional ARX - NARX identification (simulation - O(h))
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
RandomExc = 0;
if RandomExc == 1
    load('Results\SDOF_Cubic_RandomExc_FineGrid.mat')
    N = 500;
else
    load('Results\SDOF_Cubic_SweepSine_FineGrid.mat')
    N = 500;
end

options.focus = 'simulation';
options.PEinit = 'y';
options.nlopts = optimset('Algorithm','Levenberg-Marquardt' ,'Display','none',...
        'TolFun',1e-3,'TolX',1e-8,'MaxFunEval',5000,'MaxIter',5000);
    
% Model orders
na = 2;     nb = 1;     nd = 1;
ARXReg{1} = 'ones(length(t),size(Y,2))';
ARXReg{2} = '(Y(t-1,:).^1)';
ARXReg{3} = '(Y(t-2,:).^1)';
ARXReg{4} = '(X(t-1,:).^1)';

NARXReg{1} = 'ones(length(t),size(Y,2))';
NARXReg{2} = '(Y(t-1,:).^1)';
NARXReg{3} = '(Y(t-2,:).^1)';
NARXReg{4} = '(X(t-1,:).^1)';
NARXReg{5} = '(Y(t-1,:).^3)';


disp('Linear    Nonlinear')
disp('===================')
m = 1; n =0;
for ii = 1:K
    if RandomExc == 1
        X = force(1001:1000+N,ii);
        Y = dspl(1001:1000+N,ii);
    else
        X = force(1:N,ii);
        Y = dspl(1:N,ii);
    end
    X = X-mean(X);
    Y = Y-mean(Y);
    n = n+1;
    if n > indxF
        n = 1;
        m = m +1;
    end
    % ARX estimation
    [theta(:,ii),res{ii},criteria{ii}] = narx(Y,X,[na nd nb],ARXReg,[],[],[],[],options);
    rss_sss(ii) = criteria{ii}.rss_sss;
    % TH0 = [theta1(m,n) theta2(m,n) theta3(m,n) theta4(m,n)]';
    % NARX estimation
    [thetaNL(:,ii),NLres{ii},NLcriteria{ii}] = narx(Y,X,[na nd nb],NARXReg,[],[],[],[],options);
    NLrss_sss(ii) = NLcriteria{ii}.rss_sss;
    s = sprintf('%2.5f \t %2.5f \t (knl = %2.4f, Fmax = %2.4f)',rss_sss(ii),NLrss_sss(ii),knl(m),Fmax(n));
    disp(s)    
end

if RandomExc == 1
    save('Results\SDOF_Cubic_RandomExc_FineGrid_NARX_simulation.mat','theta',...
        'thetaNL','criteria','NLcriteria','res','NLres','N','options','ARXReg','NARXReg')
else
    save('Results\SDOF_Cubic_SweepSine_FineGrid_NARX_simulation.mat','theta',...
        'thetaNL','criteria','NLcriteria','res','NLres','N','options','ARXReg','NARXReg')
end



%% Transient analysis (random excitation) -- LHS
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
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
IC = zeros(2,1);
% Model properties
Sys.m = 1;
Sys.k = 5000;
Sys.c = 10;

% Filter construction
forder = cheb2ord(0.15, 0.2, 0.1, 60);
[Z,P,G] = cheby2(forder,60,0.2);                       % 0.5 Hz
[numf,denf] = zp2tf(Z,P,G);

% Latin Hypercube Sampling
rng(1)
Q = 2*(lhsdesign(K,M,'criterion','maximin','iterations',10)-0.5);
Qtrans(:,1) = 2500*Q(:,1) + 2500;           % knl
Qtrans(:,2) = 2500*Q(:,2) + 2500;           % Fmax

[dspl,force,acc] = deal(zeros(1500,K));
for k = 1:K
    % Model properties
    disp(k)
    Sys.knl = Qtrans(k,1);
%     stream = RandStream('mt19937ar','seed',1000+k);
%     RandStream.setDefaultStream(stream);
    rng(k+100)
    F = Qtrans(k,2)*randn(N+1000,1);
    F = filtfilt(numf,denf,F);
    F = F(501:end-500);
    [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F,Sys),Tspan,IC,options);
    for t = 1:length(Tspan)
        ACC(t,:) = cubicstiffnessSDOF(Tspan(t),Xsim(t,:),Tspan(t),F(t),Sys);
    end
    dspl(:,k) = Xsim(1501:3000,1);
    velo(:,k) = Xsim(1501:3000,2);
    acc(:,k) = ACC(1501:3000,2);
    force(:,k) = F(1501:3000);
    resforce(:,k) = Sys.k*dspl(:,k)+Sys.knl*(dspl(:,k).^3);
end

figure(1)
clr = colormap(jet(K));
subplot(411),hold on
subplot(412),hold on
subplot(413),hold on
subplot(414),hold on
for i = 1:K
    [Pf,F] = pwelch(force(:,i)-mean(force(:,i)),128,120,256,1/Ts);
    [Py,F] = pwelch(dspl(:,i)-mean(dspl(:,i)),128,120,256,1/Ts);
    [Txy,F] = tfestimate(force(:,i)-mean(force(:,i)),dspl(:,i)-mean(dspl(:,i)),128,120,256,1/Ts);
    [Cxy,F] = mscohere(force(:,i)-mean(force(:,i)),dspl(:,i)-mean(dspl(:,i)),128,120,256,1/Ts);
    subplot(411),plot(F,20*log10(abs(Pf)),'color',clr(i,:))
    subplot(412),plot(F,20*log10(abs(Py)),'color',clr(i,:))
    subplot(413),plot(F,20*log10(abs(Txy)),'color',clr(i,:))
    subplot(414),plot(F,Cxy,'color',clr(i,:))
end
save('Results\SDOF_Cubic_RandomExc_LHS.mat','force','acc','dspl','velo','resforce','N','Sys','Ts','Tspan','K','M','Q*')


%% Transient analysis (random excitation) -- LHS -- validation
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
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
IC = zeros(2,1);
% Model properties
Sys.m = 1;
Sys.k = 5000;
Sys.c = 10;

% Filter construction
forder = cheb2ord(0.15, 0.2, 0.1, 60);
[Z,P,G] = cheby2(forder,60,0.2);                       % 0.5 Hz
[numf,denf] = zp2tf(Z,P,G);

% Latin Hypercube Sampling
rng(2)
Q = 2*(lhsdesign(K,M,'criterion','maximin','iterations',10)-0.5);
Qtrans(:,1) = 2500*Q(:,1) + 2500;           % knl
Qtrans(:,2) = 2500*Q(:,2) + 2500;           % Fmax

[dspl,force,acc] = deal(zeros(1500,K));
for k = 1:K
    % Model properties
    disp(k)
    Sys.knl = Qtrans(k,1);
%     stream = RandStream('mt19937ar','seed',1000+k);
%     RandStream.setDefaultStream(stream);
    rng(k+1000)
    F = Qtrans(k,2)*randn(N+1000,1);
    F = filtfilt(numf,denf,F);
    F = F(501:end-500);
    [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F,Sys),Tspan,IC,options);
    for t = 1:length(Tspan)
        ACC(t,:) = cubicstiffnessSDOF(Tspan(t),Xsim(t,:),Tspan(t),F(t),Sys);
    end
    dspl(:,k) = Xsim(1501:3000,1);
    velo(:,k) = Xsim(1501:3000,2);
    acc(:,k) = ACC(1501:3000,2);
    force(:,k) = F(1501:3000);
    resforce(:,k) = Sys.k*dspl(:,k)+Sys.knl*(dspl(:,k).^3);
end

figure(1)
clr = colormap(jet(K));
subplot(411),hold on
subplot(412),hold on
subplot(413),hold on
subplot(414),hold on
for i = 1:K
    [Pf,F] = pwelch(force(:,i)-mean(force(:,i)),128,120,256,1/Ts);
    [Py,F] = pwelch(dspl(:,i)-mean(dspl(:,i)),128,120,256,1/Ts);
    [Txy,F] = tfestimate(force(:,i)-mean(force(:,i)),dspl(:,i)-mean(dspl(:,i)),128,120,256,1/Ts);
    [Cxy,F] = mscohere(force(:,i)-mean(force(:,i)),dspl(:,i)-mean(dspl(:,i)),128,120,256,1/Ts);
    subplot(411),plot(F,20*log10(abs(Pf)),'color',clr(i,:))
    subplot(412),plot(F,20*log10(abs(Py)),'color',clr(i,:))
    subplot(413),plot(F,20*log10(abs(Txy)),'color',clr(i,:))
    subplot(414),plot(F,Cxy,'color',clr(i,:))
end
save('SDOF_Cubic_RandomExc_LHS_validation.mat','force','acc','dspl','velo','resforce','N','Sys','Ts','Tspan','K','M','Q*')






%% Transient analysis (random excitation) -- LHS
% =========================================================================
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
K = 100; 

% ODE45 parameters
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
IC = zeros(2,1);
% Model properties
Sys.m = 1;
Sys.k = 5000;
Sys.c = 10;

% Filter construction
forder = cheb2ord(0.15, 0.2, 0.1, 60);
[Z,P,G] = cheby2(forder,60,0.2);                       % 0.5 Hz
[numf,denf] = zp2tf(Z,P,G);

% Latin Hypercube Sampling
rng(1000)
Q = 2*(lhsdesign(K,M,'criterion','maximin','iterations',10)-0.5);
Qtrans(:,1) = 2500*Q(:,1) + 2500;           % knl
Qtrans(:,2) = 2500*Q(:,2) + 2500;           % Fmax

[dspl,force,acc] = deal(zeros(1000,K));
for k = 1:K
    % Model properties
    disp(k)
    Sys.knl = Qtrans(k,1);
%     stream = RandStream('mt19937ar','seed',1000+k);
%     RandStream.setDefaultStream(stream);
    rng(k+100)
    F = Qtrans(k,2)*randn(N,1);
    F = filtfilt(numf,denf,F);
    tic;
    [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F,Sys),Tspan,IC,options);
    cputime(k) = toc;
    for t = 1:length(Tspan)
        ACC(t,:) = cubicstiffnessSDOF(Tspan(t),Xsim(t,:),Tspan(t),F(t),Sys);
    end
    dspl(:,k) = Xsim(:,1);
    velo(:,k) = Xsim(:,2);
    acc(:,k) = ACC(:,2);
    force(:,k) = F(:);
    resforce(:,k) = Sys.k*dspl(:,k)+Sys.knl*(dspl(:,k).^3);
end
save('Results\SDOF_Cubic_Computational_Time.mat','cputime','force','dspl','N','Sys','Ts','Tspan','K','M','Q*')



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
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
IC = zeros(2,1);
% Model properties
Sys.m = 1;
Sys.k = 5000;
Sys.c = 10;

% Filter construction
forder = cheb2ord(0.15, 0.2, 0.1, 60);
[Z,P,G] = cheby2(forder,60,0.2);                       % 0.5 Hz
[numf,denf] = zp2tf(Z,P,G);

% Equispaced grid
Q(:,1) = kron([-1 -0.5 0 0.5 1],ones(1,4));  % knl
Q(:,2) = kron(ones(1,5),[-0.5 0 0.5 1]);  % Fmax
Qtrans(:,1) = 2500*Q(:,1) + 2500;           
Qtrans(:,2) = 2500*Q(:,2) + 2500;           

[dspl,force,acc] = deal(zeros(1500,K));
for k = 1:K
    % Model properties
    disp(k)
    Sys.knl = Qtrans(k,1);
%     stream = RandStream('mt19937ar','seed',1000+k);
%     RandStream.setDefaultStream(stream);
    rng(k+2000)
    F = Qtrans(k,2)*randn(N+1000,1);
    F = filtfilt(numf,denf,F);
    F = F(501:end-500);
    [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F,Sys),Tspan,IC,options);
    for t = 1:length(Tspan)
        ACC(t,:) = cubicstiffnessSDOF(Tspan(t),Xsim(t,:),Tspan(t),F(t),Sys);
    end
    dspl(:,k) = Xsim(1501:3000,1);
    velo(:,k) = Xsim(1501:3000,2);
    acc(:,k) = ACC(1501:3000,2);
    force(:,k) = F(1501:3000);
    resforce(:,k) = Sys.k*dspl(:,k)+Sys.knl*(dspl(:,k).^3);
end

figure(1)
clr = colormap(jet(K));
subplot(411),hold on
subplot(412),hold on
subplot(413),hold on
subplot(414),hold on
for i = 1:K
    [Pf,F] = pwelch(force(:,i)-mean(force(:,i)),128,120,256,1/Ts);
    [Py,F] = pwelch(dspl(:,i)-mean(dspl(:,i)),128,120,256,1/Ts);
    [Txy,F] = tfestimate(force(:,i)-mean(force(:,i)),dspl(:,i)-mean(dspl(:,i)),128,120,256,1/Ts);
    [Cxy,F] = mscohere(force(:,i)-mean(force(:,i)),dspl(:,i)-mean(dspl(:,i)),128,120,256,1/Ts);
    subplot(411),plot(F,20*log10(abs(Pf)),'color',clr(i,:))
    subplot(412),plot(F,20*log10(abs(Py)),'color',clr(i,:))
    subplot(413),plot(F,20*log10(abs(Txy)),'color',clr(i,:))
    subplot(414),plot(F,Cxy,'color',clr(i,:))
end
save('Results\SDOF_Cubic_RandomExc_validation.mat','force','acc','dspl','velo','resforce','N','Sys','Ts','Tspan','K','M','Q*')


%% NARX regressors selection (GA - rss)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results\SDOF_Cubic_RandomExc_LHS.mat')

%Number of samples
N = 1000;
% Simulation experiment
SimExp = 19;
% Load data
X = force(1:N,:);
Y = dspl(1:N,:);
% Model orders 
na = 2;     nb = 2;     nd = 0;
% Polynomial regressors
P = 3;      maxterms = 1;
regressors = polyregNARX([na nd nb],P,1,maxterms);
% Optimization options
options.warnings = 'off';
options.GA.PopulationSize = 100 ;
options.maxsize = 1e8;
   
rng(100);
options.focus = 'simulation';
options.criterion = 'rss';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y(1:N,SimExp),X(1:N,SimExp),[na nd nb],regressors,options);

options.focus = 'prediction';
options.Nr = length(GAreg);
[SEstructure,IndxRem] = SEreduction(Y(1:N,SimExp),X(1:N,SimExp),[na nd nb],GAreg,options);

save('Results\SDOF_Cubic_RandomExc_GA_rss.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')



%% NARX regressors selection (GA - mnse)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results\SDOF_Cubic_RandomExc_LHS.mat')

%Number of samples
N = 1000;
% Simulation experiment
SimExp = 19;
% Load data
X = force(1:N,:);
Y = dspl(1:N,:);
% Model orders 
na = 2;     nb = 2;     nd = 0;
% Polynomial regressors
P = 3;      maxterms = 1;
regressors = polyregNARX([na nd nb],P,1,maxterms);
% Optimization options
options.warnings = 'off';
options.GA.PopulationSize = 100 ;
options.maxsize = 1e8;
   
rng(100);
options.focus = 'simulation';
options.criterion = 'mnse';
[GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y(1:N,SimExp),X(1:N,SimExp),[na nd nb],regressors,options);

options.focus = 'prediction';
options.Nr = length(GAreg);
[SEstructure,IndxRem] = SEreduction(Y(1:N,SimExp),X(1:N,SimExp),[na nd nb],GAreg,options);

save('Results\SDOF_Cubic_RandomExc_GA_mnse.mat','SEstr*','GAreg','indx','GAoutput','regressors','options')



%% Local NARX model estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results\SDOF_Cubic_RandomExc_LHS.mat')

%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = dspl(1:N,:);
% Model orders 
na = 2;     nb = 2;     nd = 0;
% Polynomial regressors
P = 3;      maxterms = 1;
regressors = polyregNARX([na nd nb],P,1,maxterms);
SelReg = regressors([2 3 5 12]);
orders = [na nd nb];

clear *criteria* *res* theta*
for k=1:K
    disp(k)
    options.focus = 'prediction';
    options.method = 'ols';
    [thetaPE(:,k),NLresPE{k},PEcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    options.focus = 'simulation';
    options.PEinit = 'y';
    options.nlmethod = 'LM';
    options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','final',...
            'TolFun',1e-6,'TolX',1e-9,'MaxFunEval',10000,'MaxIter',10000);
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
save('Results\SDOF_Cubic_RandomExc_LocalNARX_lsqnonlin.mat','*criteria','NLres*','theta*','SelReg','orders')



%% PC-NARX estimation (lsqnonlin)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results\SDOF_Cubic_RandomExc_LHS.mat')
load('Results\SDOF_Cubic_RandomExc_LocalNARX_lsqnonlin.mat')
%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = dspl(1:N,:);

options.basis = 'legen';
% Polynomial regressors
P = 4;      maxterms = 1;
% INDX{1} = [0     0;
%      0     1;
%      1     0;
%      1     1];

INDX{1} = [0     0;
        0     1;
        1     0;    
        0     2;
        3     0;
        4     0;
        3     1];
 
% Total number of basis functions
B = size(INDX{1},1);
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
options.focus = 'prediction';
options.method = 'ols';
options.focus = 'simulation';
options.PEinit = 'n';
options.nlmethod = 'LM';

% Model orders 
na = 2;     nb = 2;     nd = 0;
% NL regressors
regressors = polyregNARX([na nd nb],P,1,maxterms);
% SelReg = regressors([2 3 5 12]);
orders = [na nd nb];
[TH,res,criteria,output] = pcnarx(Y,X,Q,[na nd nb],INDX,SelReg,THij,options);
save('Results\SDOF_Cubic_RandomExc_PCNARX_lsqnonlin.mat','criteria','res','TH','output','orders','SelReg','INDX','THij','options')


%% Local NARX model estimation (fminsearch)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results\SDOF_Cubic_RandomExc_LHS.mat')

%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = dspl(1:N,:);
% Model orders 
na = 2;     nb = 2;     nd = 0;
% Polynomial regressors
P = 3;      maxterms = 1;
regressors = polyregNARX([na nd nb],P,1,maxterms);
SelReg = regressors([2 3 5 12]);
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
save('Results\SDOF_Cubic_RandomExc_LocalNARX_fminsearch.mat','*criteria','NLres*','theta*','SelReg','orders')


%% PC-NARX estimation (fminsearch)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results\SDOF_Cubic_RandomExc_LHS.mat')
load('Results\SDOF_Cubic_RandomExc_LocalNARX_fminsearch.mat')
%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = dspl(1:N,:);

options.basis = 'legen';
% Polynomial regressors
P = 3;      maxterms = 1;
INDX{1} = [0     0;
     0     1;
     1     0;
     1     1];
% Total number of basis functions
B = size(INDX{1},1);
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
options.focus = 'prediction';
options.method = 'ols';
options.focus = 'simulation';
options.PEinit = 'n';
options.nlmethod = 'FMINSEARCH';

% Model orders 
na = 2;     nb = 2;     nd = 0;
% NL regressors
regressors = polyregNARX([na nd nb],P,1,maxterms);
SelReg = regressors([2 3 5 12]);
orders = [na nd nb];

[TH,res,criteria,output] = pcnarx(Y,X,Q,[na nd nb],INDX,SelReg,THij,options);
save('Results\SDOF_Cubic_RandomExc_PCNARX_fminsearch.mat','criteria','res','TH','output','orders','SelReg','INDX','THij','options')



%% Local NARX model estimation (fminunc)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('Results\SDOF_Cubic_RandomExc_LHS.mat')

%Number of samples
N = 1000;
% Load data
X = force(1:N,:);
Y = dspl(1:N,:);
% Model orders 
na = 2;     nb = 2;     nd = 0;
% Polynomial regressors
P = 3;      maxterms = 1;
regressors = polyregNARX([na nd nb],P,1,maxterms);
SelReg = regressors([2 3 5 12]);
orders = [na nd nb];

clear *criteria* *res* theta*
for k=1:K
    disp(k)
    options.focus = 'prediction';
    options.method = 'ols';
    [thetaPE(:,k),NLresPE{k},PEcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    options.focus = 'simulation';
    options.PEinit = 'y';
    options.nlmethod = 'FMINUNC';
    [thetaSIM(:,k),NLresSIM{k},SIMcriteria{k}] = narx(Y(:,k),X(:,k),[na nd nb],SelReg ,[],[],[],[],options);
    rssPE(k) = PEcriteria{k}.rss_sss;
    rssSIM(k) = SIMcriteria{k}.rss_sss;
    mnseSIM(k) = SIMcriteria{k}.mnse;
end
save('Results\SDOF_Cubic_RandomExc_LocalNARX_fminunc.mat','*criteria','NLres*','theta*','SelReg','orders')



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
load('Results/SDOF_Cubic_RandomExc_LHS.mat');

figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(6,1,1:3)
hl1 = line(1:K,Qtrans(:,1),'Color','b','Marker','o');
ax1 = gca;
set(ax1,'XColor','k','YColor','k','xlim',[0.9,20.1],'ylim',[0 5000],'YTick',[0:1000:5000],...
        'xtick',[0:4:20],'xticklabel',[],'Fontname','TimesNewRoman','Fontsize',9)
ylabel('$k_{\mathrm{nl}}$','Fontsize',14,'Interpreter','Latex')
hleg1 = legend({'$k_{\mathrm{nl}}$'},'Interpreter','Latex','Fontsize',12);
set(hleg1,'Position',[0.105 0.94 0.4 0.05]);
ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','right',...
           'Color','none','XColor','k','YColor','k','XTick',[0:4:20],'YTick',[0:1000:5000],...
           'xlim',[0.9,20.1],'ylim',[0 5000],'Fontname','TimesNewRoman','Fontsize',9);
hl2 = line(1:K,Qtrans(:,2),'Color','r','Marker','d','Linestyle','--','Parent',ax2);
box on
hleg2 = legend({'$\sigma_{\!\! _x}$'},'Interpreter','Latex','Fontsize',12);
set(hleg2,'Position',[0.53 0.94 0.4 0.05]);
xlabel('Simulation experiment number','Fontsize',10,'Fontname','TimesNewRoman')
ylabel('$\sigma_{\!\! _x}$','Fontsize',14,'Interpreter','Latex','Rotation',270)
grid on
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'cubic_uncprops']),close;
     if ispc
        result = eps2xxx([write_dir,'cubic_uncprops.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
        result = eps2xxx([write_dir,'cubic_uncprops.eps'],{'pdf'});
     end
end



%% Simulated responses -- Figures
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
load('Results/SDOF_Cubic_RandomExc_LHS.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';

ExpLow = 18;
ExpHigh = 19;

figure(1)
clr = colormap(lines(2));
hold on
figindx = [1 1 2 2];
colindx = [1 2 1 2];
linindx = [1 1 1 1];
subplot(2,2,1),plot(dspl(1:1000,ExpLow),resforce(1:1000,ExpLow)/1000,'.',...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
hold on
str = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(ExpLow,1))),'\ (N/m^3), \sigma_x = ',num2str(round(Qtrans(ExpLow ,2))),'\ (N)$']};
legend(str,'Fontsize',11,'Location','Northoutside','Interpreter','Latex')
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Force (kN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,2,5),plot(0.005*(0:999),dspl(1:1000,ExpLow),...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
set(gca,'xticklabel',[])
ylabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,2,7),plot(0.005*(0:999),resforce(1:1000,ExpLow)/1000,...
    'Linewidth',linindx(1),'Color',clr(colindx(1),:))
grid on;
ylabel('Force (kN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(2,2,2),plot(dspl(1:1000,ExpHigh),resforce(1:1000,ExpHigh)/1000,'.',...
    'Linewidth',linindx(1),'Color',clr(colindx(2),:))
hold on
grid on
str = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(ExpHigh,1))),'\ (N/m^3), \sigma_x = ',num2str(round(Qtrans(ExpHigh,2))),'\ (N)$']};
legend(str,'Fontsize',11,'Location','Northoutside','Interpreter','Latex')
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Force (kN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(4,2,6),plot(0.005*(0:999),dspl(1:1000,ExpHigh),...
    'Linewidth',linindx(2),'Color',clr(colindx(2),:))
set(gca,'xticklabel',[])
ylim([-2 2])
grid on
ylabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
subplot(4,2,8),plot(0.005*(0:999),resforce(1:1000,ExpHigh)/1000,...
    'Linewidth',linindx(2),'Color',clr(colindx(2),:))
% axis tight
ylim([-25 25])
grid on
ylabel('Force (kN)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_randexci']),close;
    if ispc
        result = eps2xxx([write_dir,'cubic_randexci.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_randexci.eps'],{'pdf'});
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
load('Results/SDOF_Cubic_RandomExc_GA_rss.mat')
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
axis([1 10 0 100])
grid on
text(NoReg/2+0.5,-22.5,'Regressors dropped','Fontangle','normal','Fontsize',10,'Horizontalalignment','center','Fontname','TimesNewRoman')
ylabel('NSSE (%)','Fontangle','normal','Fontsize',10,'Fontname','TimesNewRoman')
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
            load('Results/SDOF_Cubic_RandomExc_LocalNARX_lsqnonlin.mat')
        case 2
            load('Results/SDOF_Cubic_RandomExc_LocalNARX_fminunc.mat')
        case 3
            load('Results/SDOF_Cubic_RandomExc_LocalNARX_fminsearch.mat')
    end
    
    for k=1:K
        rssSIM(j,k) = SIMcriteria{k}.rss_sss;
        mnseSIM(j,k) = SIMcriteria{k}.mnse;
        exitflag(j,k) = SIMcriteria{k}.SIMexitflag;
    end
end
subplot(211),plot(1:K,rssSIM),grid on
legend({'lsqnonlin','fminunc','fminsearch'})
subplot(212),plot(1:K,mnseSIM),grid on




%% Parameters expansion -- Trial and Error -- Figure
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end
load('Results/SDOF_Cubic_RandomExc_LHS.mat')
load('Results/SDOF_Cubic_RandomExc_LocalNARX_lsqnonlin.mat')
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

clear basis options
options.basis = 'legen';
options.criterion = 'ar2';
Pmax = 4;
for P = 0:Pmax
    % Basis index
    options.Nb = 0;
    [BASISstructure,INDnew] = PCEreduction(Q',P,thetaSIM,options);
    R2o(P+1) = BASISstructure.initvalue;
end

options.Nb = 9; Nb = 9;
% Basis index
P = 4;
[BASISstructure,INDnew] = PCEreduction(Q',P,thetaSIM,options);

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(221),plot(0:Pmax,R2o,'-o','Markersize',6,'Linewidth',1.2)
grid on
ylabel('Mean normalized error','Fontsize',11)
xlabel(sprintf('%s\n%s','Total PC basis degree','(complete subspace)'),'Fontsize',10)
set(gca,'yscale','linear','Fontsize',7,'xtick',[0:4])
subplot(222),plot(1:Nb,BASISstructure.values(1:Nb),'-o','Markersize',6,'Linewidth',1.2)
hold on 
plot([1 Nb],BASISstructure.initvalue*[1 1],'--r','Markersize',6,'Linewidth',1.2)
ylabel('Mean normalized error','Fontsize',11)
set(gca,'yscale','linear','Fontsize',7,'xtick',[1:Nb],'xticklabel',[])
for ii = 1:Nb 
    text(ii,-0.5e-3,['$[',num2str(BASISstructure.removed(ii,1)),',',num2str(BASISstructure.removed(ii,2)),']$'],'Interpreter','Latex','Rotation',45,...
        'Horizontalalignment','center','Interpreter','Latex','Fontsize',9)
end
% axis([0.99 2.01 0.9997 1])
grid on
text(1.5,-0.8e-3,sprintf('%s\n%s','PC bases dropped','(multivariable indeces)'),'Fontsize',10,'HorizontalAlignment','center')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'cubic_MSS_StageB']),close;
     if ispc
        result = eps2xxx([write_dir,'cubic_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
        result = eps2xxx([write_dir,'cubic_MSS_StageB.eps'],{'pdf'});
     end
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
load('Results/SDOF_Cubic_RandomExc_LHS.mat')
load('Results/SDOF_Cubic_RandomExc_LocalNARX_lsqnonlin.mat')
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
    text(ii,-0.2,txtstr,'Rotation',25,'Fontsize',7,'Horizontalalignment','center','Interpreter','Latex','Fontsize',10)
end
axis([1 NoB -0.1 0.8])
grid on
text(NoB/2+0.5,-0.3,'Basis functions dropped','Fontangle','normal','Fontsize',10,'Horizontalalignment','center','Fontname','TimesNewRoman')
% ylabel('mean(R^2_{adj})','Fontangle','normal','Fontsize',9)
ylabel('$\overline{R^2_{\mbox{adj}}}$','Fontangle','normal','Fontsize',11,'Interpreter','Latex')
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_MSS_StageB']),close;    
    if ispc
        result = eps2xxx([write_dir,'cubic_MSS_StageB.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_MSS_StageB.eps'],{'pdf'});
    end
end


%% Local NARX vs Global PC-NARX -- Figures
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
load('Results/SDOF_Cubic_RandomExc_LocalNARX_lsqnonlin.mat')
for k=1:K
    rss(1,k) = SIMcriteria{k}.rss_sss;
    mnse(1,k) = SIMcriteria{k}.mnse;
    R2(1,k) = 1- SIMcriteria{k}.rss_sss/100;
    exitflag(1,k) = SIMcriteria{k}.SIMexitflag;
end

load('Results/SDOF_Cubic_RandomExc_LHS.mat')
load('Results/SDOF_Cubic_RandomExc_PCNARX_lsqnonlin.mat')
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
axis([0.9 20.1 0.999 1])
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_local_vs_global']),close;
    if ispc
        result = eps2xxx([write_dir,'cubic_local_vs_global.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_local_vs_global.eps'],{'pdf'});
    end
end



%% Global PC-NARX [estimation vs validation set] -- Figures
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
print_to_eps = 'n';

K = 20;
load('Results/SDOF_Cubic_RandomExc_LHS.mat')
load('Results/SDOF_Cubic_RandomExc_PCNARX_lsqnonlin.mat')
clear res
options.basis = 'legen';
Pmax = 4;
for ii = 1:K
    [~,an(ii,:)] = PCparam(Q(ii,:)',Pmax,INDX{1},TH.sim(:),options);
    X = force(1:1000,ii);
    Y = dspl(1:1000,ii);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    [Yp,E(:,ii)] = narxpred(Y,X,2,SelReg,an(ii,:));
    res(:,ii) = Y(3:end) - Ysim(3:end);
    rss(ii) = 100*norm(res(:,ii))^2/norm(Y(3:end))^2;
    R2(ii) = 1 - rss(ii)/100;
    mnse(ii) = mean(abs(res(:,ii))./(1+abs(Y(3:end))));
end

clear res X Y Q
load('Results/SDOF_Cubic_RandomExc_validation.mat')
for ii = 1:K
    [~,an(ii,:)] = PCparam(Q(ii,:)',Pmax,INDX{1},TH.sim(:),options);
    X = force(1:1000,ii);
    Y = dspl(1:1000,ii);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    [Yp,E(:,ii+K)] = narxpred(Y,X,2,SelReg,an(ii,:));
    res(:,ii) = Y(3:end) - Ysim(3:end);
    rss(ii+K) = 100*norm(res(:,ii))^2/norm(Y(3:end))^2;
    R2(ii+K) = 1 - rss(ii+K)/100;
    mnse(ii+K) = mean(abs(res(:,ii))./(1+abs(Y(3:end))));
end

figure(1)
subplot(4,8,[1:4 9:12]),plot(1:K,R2(1:K),'-bo')
legend({'Estimation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
axis([1 20 0.997 1])
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
ylabel('R^2','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(4,8,[5:8 13:16]),plot(1:K,R2(K+1:2*K),'-rd')
legend({'Validation set'},'Fontsize',9,'Location','Northoutside','Orientation','Horizontal')
xlabel('Simulation experiment','Fontsize',11,'Fontname','TimesNewRoman')
set(gca,'Fontsize',9,'Fontname','TimesNewRoman','yticklabel',[])
axis([1 20 0.99 1])
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_est_vs_val']),close;
    if ispc
        result = eps2xxx([write_dir,'cubic_est_vs_val.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_est_vs_val.eps'],{'pdf'});
    end
end



%% Global PC-NARX [estimation vs validation set] -- Figures
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
    write_dir = '/home/minas/Dropbox/MATLAB/PCNARXjournal/Figures/';
end

K = 100;
load('Results/SDOF_Cubic_Computational_Time.mat')
load('Results/SDOF_Cubic_RandomExc_PCNARX_lsqnonlin.mat')
clear res
options.basis = 'legen';
Pmax = 4;
for ii = 1:K
    ii
    X = force(1:1000,ii);
    Y = dspl(1:1000,ii);
    tic
    [~,an(ii,:)] = PCparam(Q(ii,:)',Pmax,INDX{1},TH.sim(:),options);
    Ysim = narxsim(X,Y(1:2),2,SelReg,an(ii,:)');
    cputimePCNARX(ii) = toc;
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
load('Results/SDOF_Cubic_RandomExc_PCNARX_lsqnonlin.mat')
load('Results/SDOF_Cubic_RandomExc_validation.mat')
clear res
options.basis = 'legen';
Pmax = 4;
ii=1;
for k = [15 20]
    [~,an(ii,:)] = PCparam(Q(k,:)',Pmax,INDX{1},TH.sim(:),options);
    X(:,ii) = force(1:1000,k);
    Y(:,ii) = dspl(1:1000,k);
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
strT = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(15,1))),'\ (N/m^3), \sigma_x = ',num2str(round(Qtrans(15,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',12)
legend(str,'Fontsize',11,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',9,'Fontname','TimesNewRoman','xticklabel',[])
% set(h1,'Position',[0.3 0.925 0.4 0.05])
ylabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
subplot(2,1,2),plot(time,Y(:,2),'-b','Linewidth',1)
hold on
plot(time,Ysim(:,2),'--r','Linewidth',1)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(20,1))),'\ (N/m^3), \sigma_x = ',num2str(round(Qtrans(20,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',12)
legend(str,'Fontsize',11,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
ylabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',11,'Fontname','TimesNewRoman')
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_randexci_val']),close;
    if ispc
        result = eps2xxx([write_dir,'cubic_randexci_val.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_randexci_val.eps'],{'pdf'});
    end
end

% Figure with errors
figure(1)
hold on
subplot(8,4,[1 2 5 6 9 10]),plot(time,Y(:,1),'-b','Linewidth',0.75)
hold on
plot(time,Ysim(:,1),'--r','Linewidth',0.75)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(15,1))),'\ (N/m^3), \sigma_x = ',num2str(round(Qtrans(15,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',8)
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','xticklabel',[],'ylim',[-1.5 2])
% set(h1,'Position',[0.3 0.925 0.4 0.05])
ylabel('Displacement (m)','Fontsize',8,'Fontname','TimesNewRoman')
grid on

subplot(8,4,[13 14 17 18]),plot(time(3:end),res(:,1),'-k','Linewidth',0.75)
str = {'PC-NARX metamodel simulation error'};
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','ylim',[-0.075 0.075])
ylabel('Error (m)','Fontsize',8,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on

subplot(8,4,[3 4 7 8 11 12]),plot(time,Y(:,2),'-b','Linewidth',0.75)
hold on
plot(time,Ysim(:,2),'--r','Linewidth',0.75)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(20,1))),'\ (N/m^3), \sigma_x = ',num2str(round(Qtrans(20,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',8)
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','xticklabel',[])
grid on

subplot(8,4,[15 16 19 20]),plot(time(3:end),res(:,2),'-k','Linewidth',0.75)
str = {'PC-NARX metamodel simulation error'};
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','ylim',[-0.075 0.075])
xlabel('Time (s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_randexci_val_error']),close;
    if ispc
        result = eps2xxx([write_dir,'cubic_randexci_val_error.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_randexci_val_error.eps'],{'pdf'});
    end
end


% Figure with errors
figure(1)
subplot(18,1,1:6),plot(time,Y(:,1),'-b','Linewidth',0.75)
hold on
plot(time,Ysim(:,1),'--r','Linewidth',0.75)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(15,1))),'\ (N/m^3), \sigma_x = ',num2str(round(Qtrans(15,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',8)
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','xticklabel',[],'ylim',[-1.5 2])
ylabel('Displacement (m)','Fontsize',8,'Fontname','TimesNewRoman')
grid on
subplot(18,1,7:8),plot(time(3:end),res(:,1),'-k','Linewidth',0.75)
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','ylim',[-0.075 0.075])
ylabel('Error (m)','Fontsize',8,'Fontname','TimesNewRoman')
xlabel('Time (s)','Fontsize',8,'Fontname','TimesNewRoman')
grid on
subplot(18,1,11:16),plot(time,Y(:,2),'-b','Linewidth',0.75)
hold on
plot(time,Ysim(:,2),'--r','Linewidth',0.75)
str = {'Numerical model','PC-NARX metamodel'};
strT = {['$k_{\mathrm{nl}} = ',num2str(round(Qtrans(20,1))),'\ (N/m^3), \sigma_x = ',num2str(round(Qtrans(20,2))),'\ (N)$']};
title(strT,'Interpreter','Latex','Fontsize',8)
legend(str,'Fontsize',7,'Orientation','Horizontal','Location','NorthEast','Fontname','TimesNewRoman');
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','xticklabel',[])
ylabel('Displacement (m)','Fontsize',8,'Fontname','TimesNewRoman')
grid on
subplot(18,1,17:18),plot(time(3:end),res(:,2),'-k','Linewidth',0.75)
set(gca,'Fontsize',7,'Fontname','TimesNewRoman','ylim',[-0.075 0.075])
xlabel('Time (s)','Fontsize',8,'Fontname','TimesNewRoman')
ylabel('Error (m)','Fontsize',8,'Fontname','TimesNewRoman')
grid on
if print_to_eps=='y';
    print(pictype,resolution,[write_dir,'cubic_randexci_val_error']),close;
    if ispc
        result = eps2xxx([write_dir,'cubic_randexci_val_error.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    else
        result = eps2xxx([write_dir,'cubic_randexci_val_error.eps'],{'pdf'});
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

load('Results/SDOF_Cubic_RandomExc_PCNARX_lsqnonlin.mat')
PCBASES = INDX{1};
NLREGRESSORS = SelReg;
options.basis = 'legen';
NoI1 = 20;
NoI2 = 40;
aux = linspace(-1,1,NoI1);
Qsurf(1,:) = kron(ones(1,NoI2),aux);
Qsurf(2,:) = kron(aux,ones(1,NoI2));
[~,an_interp] = PCparam(Qsurf,4,INDX{1},reshape(TH.sim,28,1),options);


figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
% subplot(2,2,1)
surf(linspace(0,5000,NoI2),linspace(0,5000,NoI1),reshape(an_interp(:,4),NoI1,NoI2))
shading interp
set(gca,'Fontsize',10)
ylabel('$k_{\mathrm{nl}}$','Interpreter','Latex','Fontsize',20)
xlabel('$\sigma_{x}$','Interpreter','Latex','Fontsize',20)
zlabel(['$\hat{\theta}_{y[t-\ \ 1]^3} (\xi)$'],'Interpreter','Latex','Fontsize',20)
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


% Theoretical parameter Surfaces -- [Figures]
% =========================================================================
% Random input variables
Fmax = 0:250:5000;
knl = 0:250:5000;
% Model properties
Sys.m = 1;
Sys.k = 5000;
Sys.c = 10;
% Sampling period (sampling freq. = 100 Hz)
Ts = 0.005;                 % Max. freq ~= 5 Hz

na = 2;     nb = 1;     nd = 1;
NARXReg{1} = 'ones(length(t),size(Y,2))';
NARXReg{2} = '(Y(t-1,:).^1)';
NARXReg{3} = '(Y(t-2,:).^1)';
NARXReg{4} = '(X(t-1,:).^1)';
NARXReg{5} = '(Y(t-1,:).^3)';

% Theoretical curves
theta0 = -ones(length(knl),length(Fmax))*(Sys.m + Sys.c)/(Sys.m/Ts^2 + Sys.c/(2*Ts));
theta1 = -ones(length(knl),length(Fmax))*(2*Sys.m/Ts^2 - Sys.k)/(Sys.m/Ts^2 + Sys.c/(2*Ts));
theta2 = -ones(length(knl),length(Fmax))*(-Sys.m/Ts^2 + Sys.c/(2*Ts))/(Sys.m/Ts^2 + Sys.c/(2*Ts));
theta3 = -ones(length(knl),length(Fmax))*(1/(Sys.m/Ts^2+Sys.c/(2*Ts)));
theta4 = -kron(knl'/(Sys.m/Ts^2+Sys.c/(2*Ts)),ones(1,length(Fmax)));


figure(4)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
% subplot(2,2,1)
surf(Fmax,knl,theta4)
shading interp
set(gca,'Fontsize',10)
ylabel('$k_{\mathrm{nl}}$','Interpreter','Latex','Fontsize',20)
xlabel('$\sigma_{x}$','Interpreter','Latex','Fontsize',20)
zlabel(['${\theta}_{y[t-\ \ 1]^3} (\xi)$'],'Interpreter','Latex','Fontsize',20)
grid on
box on
zlim([-0.125 0])
view([125 25])
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'cubic_surf_theor']),close;
     if ispc
        result = eps2xxx([write_dir,'cubic_surf_theor.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
     else
        result = eps2xxx([write_dir,'cubic_surf_theor.eps'],{'pdf'});
     end
end




%% Conventional ARX - NARX identification (simulation) -- Figures
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
RandomExc = 1;
if RandomExc == 1
    load('Results/SDOF_Cubic_RandomExc_FineGrid.mat','K','knl','Fmax','Sys','Ts')
    load('Results/SDOF_Cubic_RandomExc_FineGrid_NARX_simulation.mat')
else
    load('Results/SDOF_Cubic_SweepSine_FineGrid.mat','K','knl','Fmax','Sys','Ts')
    load('Results/SDOF_Cubic_SweepSine_FineGrid_NARX_simulation.mat')
end
K = length(theta);
for ii=1:K
    rss_sss(ii) = criteria{ii}.rss_sss;
    extflg(ii) = criteria{ii}.PEexitflag;
    NLrss_sss(ii) = NLcriteria{ii}.rss_sss;
    NLextflg(ii) = NLcriteria{ii}.PEexitflag;
end

figure(1)
plot(rss_sss,'-o')
hold on
plot(NLrss_sss,'r-d')
legend({'ARX(2,1,1)','NARX(2,1,1)'})
xlabel('Simulation experiment')
ylabel('RSS/SSS (%)')

% Theoretical curves
theta1 = ones(length(knl),length(Fmax))*(2*Sys.m/Ts^2 - Sys.k)/(Sys.m/Ts^2 + Sys.c/(2*Ts));
theta2 = ones(length(knl),length(Fmax))*(-Sys.m/Ts^2 + Sys.c/(2*Ts))/(Sys.m/Ts^2 + Sys.c/(2*Ts));
theta3 = ones(length(knl),length(Fmax))*(1/(Sys.m/Ts^2+Sys.c/(2*Ts)));
theta4 = -kron(knl'/(Sys.m/Ts^2+Sys.c/(2*Ts)),ones(1,length(Fmax)));


figure(2)
for k=1:size(thetaNL,1)
    q = 0;
    for i = 1:length(knl)
        for j = 1:length(Fmax)
            q = q+1;
            thetasurf(i,j,k) = thetaNL(k,q);
        end
    end
    subplot(2,3,k),surf(Fmax(1:end),knl,thetasurf(:,1:end,k))
    hold on
    if k>1
        surf(Fmax(1:end),knl,eval(['theta',num2str(k-1)]))
    end
    ylabel('k_{nl} (N/m)' )
    xlabel('std(F) (N)' )
    shading interp
end


%% Duffing oscillator (modes-Poincare) -- Figures
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
load('Results/SDOF_Cubic_Poincare.mat')
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'n';

figure(1)
clr = colormap(lines(2));
hold on
figindx = [1 1 2 2];
colindx = [1 2 1 2];
linindx = [3 1 3 1];
k = 0;
for q = 1:length(Fmax)
    for r = 1:length(knl)
        k = k+1;
        U = Fmax(q)*sin(omega(1)*Tspan);
        subplot(2,2,figindx(k)),plot(Ysim{q,r}(2001:end,1),Sys.k*Ysim{q,r}(2001:end,1) + Sys.knl*Ysim{q,r}(2001:end,1).^3,...
            'Linewidth',linindx(k),'Color',clr(colindx(k),:))
        hold on
        str{r} = ['$k_{\mathrm{nl}}$ = ',num2str(knl(r))];
    end
    legend(str,'Fontsize',9,'Location','SouthEast','Interpreter','Latex')
    set(gca,'Fontsize',9,'Fontname','TimesNewRoman')
    ylabel('Force (N)','Fontsize',11,'Fontname','TimesNewRoman')
    xlabel('Displacement (m)','Fontsize',11,'Fontname','TimesNewRoman')
    grid on
    box on
    if q ==1,   ylim=[-1 1];   else ylim = [-10 10];    end
end

if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'cubic_sinexci']),close;
     result = eps2xxx([write_dir,'cubic_sinexci.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end



figure(2)
linestyle = {'-o','-+'};
clr = colormap(lines(6));
hold on
for j = 1:length(knl)
    for i = 1:length(Fmax)
        plot3(j*ones(1,length(omega)),omega,20*log10(squeeze(A(:,i,j)))',linestyle{i},'Linewidth',1.5,'Color',clr(i,:))
    end
end
set(gca,'xtick',1:length(knl),'xticklabel',knl)
ylabel('Frequency (Hz)')
xlabel('k_{nl}')
zlabel('|A| (dB)')
legend({'F_{max} = 1','F_{max} = 10'})
view([-25 25])
grid on

figure(3)
plot(Poincare{1,2,2}(:,2),Poincare{1,2,2}(:,1),'o')



%% Transient analysis-- Fine grid -- Figures
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Desktop\PCNARXmodels_JournalPaper\mfiles\SimulationModels')

RandomExc = 1;
if RandomExc == 1
    load('Results/SDOF_Cubic_RandomExc_FineGrid.mat')
else
    load('Results/SDOF_Cubic_SweepSine_FineGrid.mat')
end

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
for i = 1:indxk
    for j = 1:indxF
        k = (i-1)*indxF + j;
        [Py,F] = pwelch(dspl(:,k),128,120,256,1/Ts);
        plot3(knl(i)*ones(length(F),1),F,20*log10(abs(Py)),'color',clr(j,:))
    end
end
view([80 40])






















%% ========================================================================
%% ========================================================================

%% Conventional ARX - NARX identification (simulation - O(h^2))
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
RandomExc = 1;
if RandomExc == 1
    load('SDOF_Cubic_RandomExc_FineGrid.mat')
    N = 1000;
else
    load('SDOF_Cubic_SweepSine_FineGrid.mat')
    N = 2000;
end

options.PEinit = 'y';
options.focus = 'simulation';
options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','none',...
        'TolFun',1e-3,'TolX',1e-8,'MaxFunEval',5000,'MaxIter',5000);
    
% Model orders
na = 4;     nb = 2;     nd = 2;
% Constant term: 'ones(length(t),size(Y,2))';
ARXReg{1} = 'ones(length(t),size(Y,2))';
ARXReg{2} = '(Y(t-1,:).^1)';
ARXReg{3} = '(Y(t-2,:).^1)';
ARXReg{4} = '(Y(t-3,:).^1)';
ARXReg{5} = '(Y(t-4,:).^1)';
ARXReg{6} = '(X(t-2,:).^1)';

NARXReg{1} = 'ones(length(t),size(Y,2))';
NARXReg{2} = '(Y(t-1,:).^1)';
NARXReg{3} = '(Y(t-2,:).^1)';
NARXReg{4} = '(Y(t-3,:).^1)';
NARXReg{5} = '(Y(t-4,:).^1)';
NARXReg{6} = '(X(t-2,:).^1)';
NARXReg{7} = '(Y(t-2,:).^3)';


disp('Linear    Nonlinear')
disp('===================')
m = 1; n =0;
for ii = 1:K
    X = force(1:N,ii);
    X = X-mean(X);
    Y = dspl(1:N,ii);
    Y = Y-mean(Y);
    n = n+1;
    if n > indxF
        n = 1;
        m = m +1;
    end
    % ARX estimation
    [theta(:,ii),res{ii},criteria{ii}] = narx(Y,X,[na nd nb],ARXReg,[],[],[],[],options);
    rss_sss(ii) = criteria{ii}.rss_sss;
    % NARX estimation
    [thetaNL(:,ii),NLres{ii},NLcriteria{ii}] = narx(Y,X,[na nd nb],NARXReg,[],[],[],[],options);
    NLrss_sss(ii) = NLcriteria{ii}.rss_sss;
    s = sprintf('%2.5f \t %2.5f \t (knl = %2.4f, Fmax = %2.4f)',rss_sss(ii),NLrss_sss(ii),knl(m),Fmax(n));
    disp(s)    
end

if RandomExc == 1
    save('SDOF_Cubic_RandomExc_FineGrid_NARX_simulation_Oh2.mat','theta',...
        'thetaNL','criteria','NLcriteria','res','NLres','N','options','ARXReg','NARXReg')
else
    save('SDOF_Cubic_SweepSine_FineGrid_NARX_simulation_Oh2.mat','theta',...
        'thetaNL','criteria','NLcriteria','res','NLres','N','options','ARXReg','NARXReg')
end


%% Conventional ARX - NARX identification (simulation - O(h^2)) -- Figures
% =========================================================================
RandomExc = 1;
if RandomExc == 1
    load('SDOF_Cubic_RandomExc_FineGrid_NARX_simulation_Oh2.mat')
else
    load('SDOF_Cubic_SweepSine_FineGrid_NARX_simulation_Oh2.mat')
end

figure(1)
plot(rss_sss,'-o')
hold on
plot(NLrss_sss,'r-d')
legend({'ARX(2,1,1)','NARX(2,1,1)'})
xlabel('Simulation experiment')
ylabel('RSS/SSS (%)')


NARXReg{1} = '(Y(t-1,:).^1)';
NARXReg{2} = '(Y(t-2,:).^1)';
NARXReg{3} = '(Y(t-3,:).^1)';
NARXReg{4} = '(Y(t-4,:).^1)';
NARXReg{5} = '(X(t-2,:).^1)';
NARXReg{6} = '(Y(t-2,:).^3)';

% Theoretical curves
Denom = -Sys.m/(12*Ts^2)- Sys.c/(12*Ts);
theta1 = ones(length(knl),length(Fmax))*(-16*Sys.m/(12*Ts^2) - 8*Sys.c/(12*Ts))/Denom;
theta2 = ones(length(knl),length(Fmax))*(30*Sys.m/(12*Ts^2) - Sys.k)/Denom;
theta3 = ones(length(knl),length(Fmax))*(-16*Sys.m/(12*Ts^2) + 8*Sys.c/(12*Ts))/Denom;
theta4 = ones(length(knl),length(Fmax))*(Sys.m/(12*Ts^2) - Sys.c/(12*Ts))/Denom;
theta5 = ones(length(knl),length(Fmax))/Denom;
theta6 = -kron(knl'/Denom,ones(1,length(Fmax)));

figure(2)
for k=1:size(thetaNL,1)
    q = 0;
    for i = 1:length(knl)
        for j = 1:length(Fmax)
            q = q+1;
            thetasurf(i,j,k) = thetaNL(k,q);
        end
    end
    subplot(2,3,k),surf(Fmax(1:end),knl,thetasurf(:,1:end,k))
    hold on
    surf(Fmax(1:end),knl,eval(['theta',num2str(k)]))
    ylabel('k_{nl} (N/m)' )
    xlabel('std(F) (N)' )
    shading interp
end



%% PC-NARX model structure selection (FOLS)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('SDOF_Cubic_RandomExc_LHS.mat')

% Nonlinear polynomial regressors
P = 3;      q = 1;      maxterms = 1;
% Model orders
na = 2;     nb = 2;     nd = 0;
N = 1000;

% Input-Output data
Y = dspl(1:N,:);
X = force(1:N,:);
% Nonlinear regressors (complete vector)
regressors = polyregNARX([na nd nb],P,q,maxterms);
% FOLS structure selection
options.criterion = 'rss';
PCorder = 3;
indx{1} = combinations(M,PCorder,1);
[Nr,RegSTR,regindx,basisindx,NLres,NLcriteria] = folsPCNARX(Y,X,Q',[na nd nb],indx,P,q,maxterms,[],20,[],options);

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
    [thetaPE,resPE{j},PEcriteria{j}] = pcnarxuncommon(Y,X,Q',[na nd nb],NLREGRESSORS,PCBASES,[],options);
    rssPE(j) = PEcriteria{j}.rss_sss;
    bicPE(j) = PEcriteria{j}.bic;
    options.focus = 'simulation';
    TH0 = thetaPE.wls;
%     if j == 1
%         TH0 = thetaPE.wls;
%     else
%         TH0 = [thetaSIM{j-1}.sim;0];
%     end
    % [resSIM{j},SIMcriteria{j}] = pcnarxuncommonSIM(thetaPE.wls,Y,X,Q',[na nd nb],NLREGRESSORS,PCBASES,options);
    [thetaSIM{j},resSIM{j},SIMcriteria{j}] = pcnarxuncommon(Y,X,Q',[na nd nb],NLREGRESSORS,PCBASES,TH0,options);
    rssSIM(j) = SIMcriteria{j}.rss_sss;
    bicSIM(j) = SIMcriteria{j}.bic;
end

figure
plot(1:20,rssPE(1:20),'-bo',1:20,rssSIM(1:20),'-rd')
set(gca,'yscale','log')

save('SDOF_Cubic_RandomExc_LHS_FOLS.mat','P','q','maxterms','na','nb','nd','N','regressors','PCorder',...
     'indx','regindx','basisindx','options','res*','rss*','bic*','theta*','*criteria')

 

 %% PC-NARX model structure selection (FOLS) -- [Figure]
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';
load('SDOF_Cubic_RandomExc_LHS_FOLS.mat')

figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(211),plot(1:length(regindx),rssPE,'-o',1:length(regindx),rssSIM,'-rd')
set(gca,'yscale','log')
grid on
set(gca,'xlim',[0.5,20.5],'ylim',[10^-6 10^3],'YTick',10.^[-6:2:3],'Fontname','TimesNewRoman','Fontsize',9)
ylabel(sprintf('%s\n%s','Normalized sum','of squared error'),'Fontsize',11,'Fontname','TimesNewRoman')
legend({'prediction error','simulation error'},'Fontname','TimesNewRoman','Fontsize',9);
xlabel('# of regressors','Fontsize',11,'Fontname','TimesNewRoman')
grid on
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'cubic_FOLS']),close;
     result = eps2xxx([write_dir,'cubic_FOLS.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end

%% PC-NARX model with 8 regressors -- Validation -- [Figure]
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';
load('SDOF_Cubic_RandomExc_LHS.mat')
load('SDOF_Cubic_RandomExc_LHS_FOLS.mat')

Y = dspl;
X = force;
ch = 7;
time = 0.01*(0:N-1);

figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(2,2,1),hold on,plot(time,Y(1:N,ch),'-b',time,Y(1:N,ch)-resPE{8}(:,ch),'--r',time,Y(1:N,ch)-resSIM{8}(:,ch),'-.g','Linewidth',1.2)
set(gca,'xlim',[0,3],'ylim',[-0.8,0.8],'YTick',[-0.8:0.4:0.8],'Fontname','TimesNewRoman','Fontsize',8)
ylabel('Displacement (m)','Fontsize',9,'Fontname','TimesNewRoman')
h1 = legend({'$y_7[t]$','$\hat{y}_7[t|t-1]$','$\bar{y}_7[t]$'},'Fontname','TimesNewRoman','Fontsize',11,'Orientation','Horizontal','Interpreter','Latex');
set(h1,'Position',[0.072 0.94 0.45 0.05])
xlabel('Time (s)','Fontsize',9,'Fontname','TimesNewRoman')
box on
grid on
subplot(2,2,2),hold on,plot(1:20,std(resPE{8}),'-o',1:20,std(resSIM{8}),'-rd','Markersize',4)
set(gca,'yscale','log','xlim',[0.5,20.5],'ylim',[10^-5,10^-2],'YTick',10.^(-5:1:-2),'Fontname','TimesNewRoman','Fontsize',8)
ylabel('$\mathrm{std}(e[t])$','Fontsize',9,'Fontname','TimesNewRoman','Interpreter','Latex')
legend({'prediction error','simulation error'},'Fontname','TimesNewRoman','Fontsize',7,'Orientation','Horizontal');
xlabel('Simulation experiment','Fontsize',9,'Fontname','TimesNewRoman')
box on
grid on
% subplot(3,2,3),hold on,acf(resPE{8}(:,ch),50,0.8)
% set(gca,'xlim',[0,100],'ylim',[-1,1],'YTick',[-1:0.4:1],'Fontname','TimesNewRoman','Fontsize',8)
% box on
% grid on
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'cubic_rss']),close;
     result = eps2xxx([write_dir,'cubic_rss.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end

figure(2),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(2,1,1),hold on,surf(corr(resPE{8}))
clr = colormap(gray);
colormap(flipud(clr));
%h2 = colorbar;
colorbar;
set(gca,'xlim',[1,20],'ylim',[1,20],'Fontname','TimesNewRoman','Fontsize',8)
% set(h2,'Position',[0.92 0.41 0.02 0.216])
ylabel('Column index','Fontsize',9,'Fontname','TimesNewRoman')
xlabel('Row index','Fontsize',9,'Fontname','TimesNewRoman')
view([90 90])
box on
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'cubic_validation']),close;
     result = eps2xxx([write_dir,'cubic_validation.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end
 

%% NARX regressors selection (PE estimation method)
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

X = force(1:N,:);
Y = acc(1:N,:); 

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
Y = acc(1:N,ii); 

% Nonlinear polynomial regressors
P = 3;      q = 1;      maxterms = 1;
% Model orders
na = 2;     nb = 2;     nd = 0;
regressors = polyregNARX([na nd nb],P,q,maxterms);
options.criterion = 'rss';
[Nr,RegFOLS,regindx,NLres,NLcriteria] = folsNARX(Y,X,[na nd nb],regressors,[],[],[],20,[],options);



%% Nonlinearity measures
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

% Load data
X = force(1:N,:);
Y = dspl(1:N,:);

% NL measure criteria
options.criterion = 'all';
options.NFFT = 256;
options.wind = hamming(256);
options.nsamples = 256;
options.overlap = 50;
% Robust bicoherence estimation [1]
options.robust = 'y';
options.percentile = 5;
options.alpha = 0.05;
for k = 1:K
    disp(k);
    [NLI(k),TNLI(k),MAXBIC(k),SUMBIC(k),MeanCoh(k)]=  NLmeasures(Y(:,k),X(:,k),options);
    % Cross-corellation(z,z^2)
    [CCF(:,k),limit] = ccf(Y(:,k)-mean(Y(:,k)),(Y(:,k)-mean(Y(:,k))).^2,100,0.8);
    % Counts instances of CCF above 95% limits
    Exc(k) = sum(CCF(:,k)>limit | CCF(:,k)<-limit);
    % Second derivative
    der(:,k) = diff(resforce(:,k))./diff(dspl(:,k));
    mder(:,k) = der(:,k) - mean(der(:,k));
    sder(k) = sum(der(:,k));
    % Time reversibility criterion
    stat(k) = timerev(Y(:,k),100,'n');
end

figure(1)
subplot(811),plot(1:K,NLI,'-bo')
ylabel('NLI')
subplot(812),plot(1:K,TNLI,'-bo')
ylabel('TNLI')
subplot(813),plot(1:K,MAXBIC,'-bo')
ylabel('max(BIC)')
subplot(814),plot(1:K,SUMBIC,'-bo')
ylabel('sum(BIC)')
subplot(815),plot(1:K,MeanCoh,'-bo')
ylabel('mean(Coh)')
subplot(816),plot(1:K,sum(abs(CCF)),'-bo')
ylabel('sum(CCF)')
subplot(817),plot(1:K,sder,'-bo')
ylabel('2nd der.')
subplot(818),plot(1:K,stat,'-bo')
ylabel('Time rev.')

figure(2)
for k = 1:K
    [bic,waxis] = bicoherx (X(:,k),X(:,k),Y(:,k),options.NFFT,options.wind,options.nsamples,options.overlap,'y',options.percentile);
    subplot(4,5,k),surf(bic)
    shading interp
    caxis([0 1])
    zlim([0 1])
end

for k=1:K; 
    figure(k+100),
    subplot(211),plot(dspl(:,k),resforce(:,k)),
    xlabel('Displacement')
    ylabel('Restoring force')
    subplot(212),plot(diff(resforce(:,k))./diff(dspl(:,k))),
    xlabel('Time')
    ylabel('df_{res}^2/d^2x')
end



%% FOLS-NARX (prediction simulation criteria)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('SDOF_Cubic_RandomExc_LHS_Ts0p005.mat')

options.maxsize = 1e7;
options.focus = 'prediction';
N = 1000;

for k = 1:K
    [Cxy,F] = mscohere(force(:,k),dspl(:,k),128,120,256,1/Ts);
    NLmeas(k) = 1 - mean(Cxy(F<40));
end
plot(1:K,NLmeas,'-o');
ii = input('Given run with stronger nonlinearities: ');
X = force(1:N,ii);
Y = dspl(1:N,ii); 

% Nonlinear polynomial regressors
P = 5;      q = 1;      maxterms = 1;
% Model orders
na = 10;     nb = 10;     nd = 0;
regressors = polyregNARX([na nd nb],P,q,maxterms);
options.criterion = 'rss';
[Nr,RegSTR,regindx,NLres,NLcriteria] = folsNARX(Y,X,[na nd nb],[],P,q,maxterms,20,[],options);

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

save('SDOF_Cubic_RandomExc_LHS_FOLS_results.mat','regressors',...
    'NLres*','regindx','rss*','bic*','na','nb','nd','P','q','maxterms','N');


%% FOLS-NARX (prediction simulation criteria) -- Figures
% =========================================================================
load('SDOF_Cubic_RandomExc_LHS_FOLS_results.mat');
n = 20;
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


%% GA PC functional subspace selection
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
options.criterion = 'bic';
options.basis = 'legen';
options.common = 'n';
options.GA.PopulationSize = 100;
options.GA.Display = 'iter';
options.GA.PlotFcns = {@gaplotbestf,@gaplotbestindiv};

FinalRegstr{1} = '(Y(t-1,:).^1)'; 
FinalRegstr{2} = '(Y(t-2,:).^1)';
FinalRegstr{3} = '(X(t-1,:).^1)'; 
FinalRegstr{4} = '(Y(t-1,:).^3)';


% PC basis characterisitcs
p = 5;      q = 1;      
% Model orders
na = 2;     nb = 1;     nd = 1;
N = 200;

P =3; q = 1;
maxterms = 2;
% [INDX{1}] = combinations(2,P,q);
% INDX{2} = INDX{1};
% [ms,criteria] = folsPCNARX(dspl(1:N,:),force(1:N,:),Q',[na nd nb],INDX,P,q,maxterms,[],20,[],options)
[INDX,thetaij,res,criteria,GAoutput] = gaPCNARX(dspl(1:N,:),force(1:N,:),Q',[na nd nb],FinalRegstr,p,q,options);
save('SDOF_Cubic_GAPCNARXresults.mat','INDX','FinalRegstr')




%% PC-NARX model structure selection (GA)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('SDOF_Cubic_RandomExc_LHS_Ts0p01.mat')

% Nonlinear polynomial regressors
P = 3;      q = 1;      maxterms = 1;
% Model orders
na = 6;     nb = 6;     nd = 0;
N = 1000;

% Input-Output data
Y = dspl(1:N,:);
X = force(1:N,:);
% Nonlinear regressors (complete vector)
regressors = polyregNARX([na nd nb],P,q,maxterms);
% FOLS structure selection
options.criterion = 'se';
options.maxsize = 10e7;
options.nlreg = 'y';
options.GA.PopulationSize = 50;
PCorder = 3;
indx{1} = combinations(M,PCorder,1);
[REGSTR,INDX,thetaij,res,criteria,GAoutput] = gaPCNARX(Y,X,Q',[na nd nb],regressors,P,q,options);


%% Run candidate models for 1-50 most important regressors
options.maxsize = 1e9;
options.wls_iter = 5;
options.method = 'wls';
options.focus = 'prediction';
for j = 1:length(regindx)
    disp(j)
    PCBASES = {indx{1}(basisindx(1:j),:)};
    NLREGRESSORS = regressors(regindx(1:j));
    [thetaNL,NLresPE{j},PEcriteria] = pcnarxuncommon(Y,X,Q',[na nd nb],NLREGRESSORS,PCBASES,[],options);
    rssPE(j) = PEcriteria.rss_sss;
    bicPE(j) = PEcriteria.bic;
    [NLresSIM{j},SIMcriteria] = pcnarxuncommonSIM(thetaNL.wls,Y,X,Q',[na nd nb],NLREGRESSORS,PCBASES,options);
    rssSIM(j) = SIMcriteria.rss_sss;
    bicSIM(j) = SIMcriteria.bic;
end


%% PC-NARX estimation (full space)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('SDOF_Cubic_RandomExc_LHS_Ts0p005.mat')

% PE estimation
options.maxsize = 1e8;
options.basis = 'legen';


N = 1000;
options.method = 'wls';
options.focus = 'prediction';
options.wls_iter = 5;

for ord = 2:8
    for pow = 3
        for pc = 1:3
            % Model orders
            na = ord;     nb = ord;     nd = 0;
            P = pow;      q = 1;      maxterms = 1;
            regressors = polyregNARX([na nd nb],P,q,maxterms);
            options.criterion = 'rss';
            indx{1} = combinations(2,pc,1);
            indx{2} = indx{1};
            [THij,res,criteria,output] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],indx,regressors ,[],options);
            rss(ord,pow,pc) = criteria.rss_sss;
        end
    end
end

%% SE estimation
na = 4;     nb = 4;     nd = 0;
P = 3;      q = 1;      maxterms = 1;
regressors = polyregNARX([na nd nb],P,q,maxterms);
options.focus = 'prediction';
options.method = 'wls';
indx{1} = combinations(2,1,1);
indx{2} = indx{1};
[THij,res,criteria,output] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],indx,regressors ,[],options);

options.focus = 'simulation';
[THij,res_sim,criteria_sim] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],indx,regressors,THij.wls,options);

save('SDOF_Cubic_PCNARXmodel.mat','TH*','criter*','res*');



%% PC-NARX estimation
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('SDOF_Cubic_RandomExc_LHS.mat')
load('SDOF_Cubic_GAPCNARXresults.mat')

% PE estimation
options.maxsize = 1e8;
options.basis = 'legen';

% Model orders
na = 2;     nb = 1;     nd = 1;

N = 1000;
options.method = 'wls';
options.focus = 'prediction';
options.wls_iter = 20;
[THij,res,criteria,output] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],INDX,FinalRegstr,[],options);
% SE estimation
options.focus = 'simulation';
[THij,res_sim,criteria_sim] = pcnarx(dspl(1:N,:),force(1:N,:),Q',[na 0 na],INDX,FinalRegstr,[],options);

save('SDOF_Cubic_PCNARXmodel.mat','TH*','criter*','res*');






%% NARX regressors selection (forward backward - rss)
% =========================================================================
clear,clc,close all
if ispc
    cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
    write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
else
    cd('/home/minas/Dropbox/MATLAB/PCNARXjournal/SimulationModels')
end
load('SDOF_Cubic_RandomExc_LHS.mat')

%Number of samples
N = 1000;
% Simulation experiment
SimExp = 19;
% Load data
X = force(1:N,:);
Y = dspl(1:N,:);
% Model orders 
na = 2;     nb = 2;     nd = 0;
% Polynomial regressors
P = 3;      maxterms = 1;
regressors = polyregNARX([na nd nb],P,1,maxterms);
   
rng(100);
options.focus = 'prediction';
options.criterion = 'rss';
options.Nr = length(regressors);
[SEstructure,IndxRem] = MSSforback(Y(1:N,SimExp),X(1:N,SimExp),[na nd nb],regressors,options);
