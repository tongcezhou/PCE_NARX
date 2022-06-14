%% Duffing oscillator (modes-Poincare) 
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')

% Number of samples
N = 2000;
% Sampling period (sampling freq. = 40 Hz)
Ts = 0.2;                 % Max. freq ~= 1 Hz
% Time index
Tspan = (0:Ts:(N-1)*Ts);

% ODE45 parameters
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
IC = zeros(4,1);
% Model properties
Sys.m1 = 200;        Sys.c1 = 0.5;     Sys.k1 = 1e4;
Sys.m2 =  20;        Sys.c2 = 0.5;     Sys.k2 = 1e2;

omega = 0.005:0.005:2;
Fmax = [0.1 1 10 100];
[Poincare1,Poincare2] = deal(cell(length(omega),length(Fmax)));
[A1,A2] = deal(zeros(length(omega),length(Fmax)));

p = 0;
for omega_cur = omega
    p = p + 1;
    q = 0;
    for Fmax_cur = Fmax
        q = q+1;
        disp([omega_cur,Fmax_cur])
        U = Fmax_cur*sin(2*pi*omega_cur*Tspan);
        [Tsim,Ysim] = ode45(@(t,x) duffingTDOF(t,x,Tspan,U,Sys),Tspan,IC,options);     
        A1(p,q) = (max(abs(Ysim(:,1)))/max(abs(U)))^2;
        A2(p,q) = (max(abs(Ysim(:,3)))/max(abs(U)))^2;
        Tsin = 1/omega_cur;
        Tstart = ceil(Tsin/(4*Ts));
        Tperiod = round(Tsin/Ts);
        Tindx = Tstart:Tperiod:N;
        Poincare1{p,q}(1:length(Tindx),1) = Ysim(Tindx,1);
        Poincare1{p,q}(1:length(Tindx),2) = Ysim(Tindx,2);
        Poincare2{p,q}(1:length(Tindx),1) = Ysim(Tindx,3);
        Poincare2{p,q}(1:length(Tindx),2) = Ysim(Tindx,4);
    end
end
save('Poincare.mat','A1','A2','Fmax','N','Poincare*','Sys','Ts','Tstart','omega')



%% Duffing oscillator (modes-Poincare) -- Figures
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
load('Poincare.mat','A1','A2','Fmax','N','Poincare*','Sys','Ts','Tstart','omega')

figure(1)
subplot(211),plot(omega,A1,'Linewidth',2)
legend({'F_{max} = 0.1', 'F_{max} = 1', 'F_{max} = 10','F_{max} = 100'})
set(gca,'yscale','log')
xlabel('Frequency (rad/s)')
ylabel('|A_1|')
grid on 
box on 
subplot(212),plot(omega,A2,'Linewidth',2)
set(gca,'yscale','log')
xlabel('Frequency (rad/s)')
ylabel('|A_2|')
grid on 
box on 

% Linear system eigenfrequencies
A = [0 1 0 0; -(Sys.k1+Sys.k2)/Sys.m1  -(Sys.c1+Sys.c2)/Sys.m1 Sys.k2/Sys.m1 Sys.c2/Sys.m1;...
    0 0 0 1 ; Sys.k2/Sys.m2 Sys.c2/Sys.m2 -Sys.k2/Sys.m2 -Sys.c2/Sys.m2];
wn = abs(eig(A));
zeta = -real(eig(A))./wn;
wn = wn/2/pi;

% Poincare plots
windx = 400;
figure(2)
subplot(221),plot(Poincare1{windx,1}(:,1),Poincare2{windx,1}(:,1),'.')
axis([-2.5e-3 2.5e-3 -0.6 +0.6])
subplot(222),plot(Poincare1{windx,2}(:,1),Poincare2{windx,2}(:,1),'.')
axis([-2.5e-3 2.5e-3 -0.6 +0.6])
subplot(223),plot(Poincare1{windx,3}(:,1),Poincare2{windx,3}(:,1),'.')
axis([-2.5e-3 2.5e-3 -0.6 +0.6])
subplot(224),plot(Poincare1{windx,4}(:,1),Poincare2{windx,4}(:,1),'.')
axis([-2.5e-3 2.5e-3 -0.6 +0.6])


%% Duffing oscillator (sweep sine) 
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')

% Number of samples
N = 20000;
% Sampling period (sampling freq. = 40 Hz)
Ts = 0.2;                 % Max. freq ~= 1 Hz
% Time index
Tspan = (0:Ts:(N-1)*Ts);

% ODE45 parameters
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
IC = zeros(4,1);
% Model properties
Sys.m1 = 200;        Sys.c1 = 0.5;     Sys.k1 = 1e4;
Sys.m2 =  20;        Sys.c2 = 0.5;     Sys.k2 = 1e2;

Fmax = [0.1 1 10 100];

q = 0;
for Fmax_cur = Fmax
    q = q+1;
    disp(Fmax_cur)
    U = Fmax_cur*chirp(Tspan,0.001,Tspan(end),2);
    [Tsim,Ysim{q}] = ode45(@(t,x) duffingTDOF(t,x,Tspan,U,Sys),Tspan,IC,options);
end

%% Duffing oscillator with uncertain input and excitation parameters 
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels')
load('Poincare.mat','Sys')
load('C:\Users\sminas\Dropbox\MATLAB\EarthquakeSimulation\SimulatedEQs.mat','Q*','EQsim','K','N','M','Ts')
% K : Number of experiments
% N: Number of samples
% M: Number of variables (beta, gamma, n)
% Ts : Sampling period
Ts = 0.2;
% Time index
Tspan = (0:Ts:(N-1)*Ts);
% ODE45 parameters
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
IC = zeros(4,1);

% Convert synthetic earthquakes in m/s^2
EQsim = EQsim*9.81;

% Initialization
[acc1,acc2,force] = deal(zeros((N),K));
for k = 1:K
    disp(['Iteration: ',num2str(k)]);
    [Tsim,Ysim] = ode45(@(t,x) duffingTDOF(t,x,Tspan,EQsim(:,k),Sys),Tspan,IC,options);     
    for t=1:N
        Ysimdot(t,:) = duffingTDOF(Tspan(t),Ysim(t,:),Tspan(t),EQsim(t,k),Sys);
    end
    acc1(:,k) = Ysimdot(1:end,2);
    acc2(:,k) = Ysimdot(1:end,4);
    force(:,k) = EQsim(:,k);
end

for k = 1:K
    [C(:,k),F] = mscohere(force(:,k),acc1(:,k),512,400,512,1/Ts);
    [Pxx(:,k),F] = pwelch(acc1(:,k),512,400,512,1/Ts);
    [Txy(:,k),F] = tfestimate(force(:,k),acc1(:,k),512,400,512,1/Ts);
end

figure(1)
plot(acc1(:,k),force(:,k))
figure(2)
plot(F,20*log10(abs(Pxx)))
figure(3)
plot(F,20*log10(abs(Txy)))
surf(1:K,F,20*log10(abs(Txy)))
shading interp
figure(4)
surf(1:K,F,C)
shading interp
save('Duffing_SimEQs.mat','Q*','acc*','force','K','N','M','Ts')



%% PC-ARX modelling
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels');
load('Duffing_SimEQs.mat','Q*','acc1','force','K','N','M','Ts')
cd('C:\Users\sminas\Dropbox\MATLAB\PCE\m-files\');

options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';
% AR modelling
% Q = Qtrans';
na = 4;
q = 1;

for p = 0:3
    disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
    indx{1} = combinations(M,p,q);
    indx{2} = indx{1};
    B = length(indx{1});
    [aij,bij,res,criteria] = pcarx(acc1,force,Q',[na 0 na],indx,[],options);
    rss_sss(p+1) = criteria.rss_sss;
    connum(p+1) = criteria.connum;
    spp(p+1) = criteria.spp;
    sigma = var(res);
    bic(p+1) = log(sigma) + ((na+na+1)*B)*log(K*N)/(K*N);
end


%% Conventional ARX - NARX identification
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels');
load('Duffing_SimEQs.mat','Q*','acc1','force','K','N','M','Ts')

na = 4;
for k = 1:K
    IdData = iddata(acc1(:,k),force(:,k),Ts);
    disp(['Experiment number: ',num2str(k)]);
    model = arx(IdData,[na na 0]);
    res = resid(model,IdData);
    rss_sss(k) = 100*norm(res.outputdata)^2/norm(acc1(:,k))^2;
    sigma(k) = var(res.outputdata);
    bic(k) = log(sigma(k)) + (2*na+1)*log(N)/N;
    nlmodel = nlarx(IdData,[na na 0],'linear','CustomReg',{'y1(t-1)','y1(t-1)^2','y1(t-1)^3',...
        'y1(t-2)','y1(t-2)^2','y1(t-2)^3','y1(t-3)','y1(t-3)^2','y1(t-3)^3','y1(t-4)','y1(t-4)^2','y1(t-4)^3',...
        'u1(t)','u1(t-1)','u1(t-2)','u1(t-3)','u1(t-4)',...
        'u1(t)^2','u1(t-1)^2','u1(t-2)^2','u1(t-3)^2','u1(t-4)^2',...
        'u1(t)^3','u1(t-1)^3','u1(t-2)^3','u1(t-3)^3','u1(t-4)^3'});
    NLres = resid(nlmodel,IdData);
    NLrss_sss(k) = 100*norm(NLres.outputdata)^2/norm(acc1(:,k))^2;
    NLsigma(k) = var(NLres.outputdata);
    NLbic(k) = log(NLsigma(k)) + (2*na+1)*log(N)/N;
end





%% Conventional ARX - NARX identification (my routine)
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels');
load('Duffing_SimEQs.mat','Q*','acc1','force','K','N','M','Ts')
cd('C:\Users\sminas\Dropbox\MATLAB\PCE\m-files\');

options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';

P = 3;
q = 1;
maxterms = [];

disp('Linear    Nonlinear')
disp('===================')
for ii = 1:K
    Y = acc1(:,ii);
    X = force(:,ii);
    na_max = 4;
    for na = 4:na_max
        [theta,res,criteria,regstr] = narx(Y,X,[na 0 na],[],1,1,[],[],options);
        rss_sss(ii,na) = criteria.rss_sss;
        sigma(na) = var(res);
        bic(na) = log(sigma(na)) + (2*na+1)*log(N)/N;
        [theta,NLres,NLcriteria,regstr] = narx(Y,X,[na 0 na],[],P,q,maxterms,[],options);
        NLrss_sss(ii,na) = NLcriteria.rss_sss;
        NLsigma(na) = var(NLres);
        NLbic(na) = log(NLsigma(na)) + (2*na+1)*log(N)/N;
        s = sprintf('%2.5f \t %2.5f',rss_sss(ii,na),NLrss_sss(ii,na));
        disp(s)
    end
    if ii==K
        % Cross-corellation(z,z^2)
        ccf(res-mean(res),(res-mean(res)).^2,100,0.8);
        acf(res-mean(res),100,0.8);
        ccf(NLres-mean(NLres),(NLres-mean(NLres)).^2,100,0.8);
        acf(NLres-mean(NLres),100,0.8);
    end
end




%% GA - NARX
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\SimulationModels');
load('Duffing_SimEQs.mat','Q*','acc1','force','K','N','M','Ts')
cd('C:\Users\sminas\Dropbox\MATLAB\PCE\m-files\');

options.focus = 'prediction';
options.GA.PopulationSize = 100 ;
options.GA.PlotFcns = {};
options.GA.Display = 'off';
P = 3;
q = 1;
maxterms = 2;

disp('Linear    Nonlinear')
disp('===================')
for ii=1:K
    Y = acc1(:,ii);
    X = force(:,ii);
    na = 4;
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
for ii=1:K
    TotIndx = [TotIndx  regindx{ii}];
end
N = hist(TotIndx,length(InitRegstr));
FinalRegstr  = InitRegstr(N==K);

% Cross-corellation(x^2,y)
figure(1)
ccf(X.^2,Y-mean(Y),100,0.8);
figure(2)
ccf(Y-mean(Y),(Y-mean(Y)).^2,100,0.8);
% Cross-corellation(z,z^2)
figure(1)
ccf(res-mean(res),(res-mean(res)).^2,100,0.8);
figure(2)
ccf(NLres-mean(NLres),(NLres-mean(NLres)).^2,100,0.8);




%% PC-NARX modelling
clear rss* connum spp bic
options.maxsize = 1e8;
options.basis = 'legen';
options.focus = 'prediction';
% AR modelling
% Q = Qtrans';
n_max = 8;
Nr = length(regstr);
q = 1;
for na = 2 : n_max
    for p = 0:3
        disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
        indx{1} = combinations(M,p,q);
        indx{2} = indx{1};
        B = length(indx{1});
        [~,~,~,InitRegstr] = narx(acc1(:,1),force(:,1),[na 0 na],[],3,1,2,[],options);
        [aij,res,criteria] = pcnarx(acc1,force,Q',[na 0 na],indx,InitRegstr,[],options);
        rss_sss(na,p+1) = criteria.rss_sss;
        connum(na,p+1) = criteria.connum;
        spp(na,p+1) = criteria.spp;
        sigma = var(res);
        bic(na,p+1) = log(sigma) + (length(InitRegstr)*B)*log(K*N)/(K*N);
    end
end



%% PC-ARX modelling (Testing FOLS routine)
cd('J:\Documents\MATLABroutines\PCE\m-files')

options.maxsize = 1e8;
options.basis = 'legen';
options.focus = 'prediction';
options.criterion = 'RSS';
options.GA.PopulationSize = 100 ;
% AR modelling
% Q = Qtrans';
q = 1;
p = 3;
indx{1} = combinations(M,p,q);
indx{2} = indx{1};
B = length(indx{1});
na = 2;
% [aij,bij,res,criteria] = pcarx(y,x,Q,[na 0 na],indx,[],options);
% [ms,criteria] = folsPCARX(y,x,Q,[na 0 na],indx,100,[],options);

[INDX,aij,bij,res,criteria,GAoutput] = gaPCARX(y,x,Q',[na 0 na],3,1,options);


%% GA PC-NARX modelling
cd('J:\Documents\MATLABroutines\PCE\m-files')

options.maxsize = 1e9;
options.basis = 'legen';
options.focus = 'prediction';
options.criterion = 'RSS';
options.GA.PopulationSize = 100 ;
options.GA.Display = 'iter';
options.GA.PlotFcns = {@gaplotbestf,@gaplotbestindiv};
% AR modelling
% Q = Qtrans';
q = 1;
p = 3;
indx{1} = combinations(M,p,q);
B = length(indx{1});
na = 2;
[INDX,aij,res,criteria,GAoutput] = gaPCNARX(y,x,Q',[na 0 na],regstr,p,1,options);



%% LHS experiments
clear,clc,close all
cd('J:\Documents\MATLABroutines\EarthquakeSimulation\');
load('J:\Documents\MATLABroutines\EarthquakeSimulation\FittedDist')

% Number of experiments
K = 10; 
% Number of samples
N = 1000;
% Number of variables (beta, gamma, n)
M = 8;
% Sampling period
Ts = 0.5;

% System constant properties
Sys.m1 = 1;         Sys.l1 = 1;
Sys.m2 = 0.05;      Sys.l2 = 0.05;
Sys.c1 = 0.0015;

% Latin Hypercube Sampling
Q = lhsdesign(K,M,'criterion','maximin','iterations',10);

Qtrans(:,1) = icdf(PD_wmid,Q(:,1))*2*pi;    % omega_mid
Qtrans(:,2) = icdf(PD_wdot,Q(:,2))*2*pi;    % omega_dot
Qtrans(:,3) = icdf(PD_zeta,Q(:,3));         % zeta
Qtrans(:,4) = icdf(PD_a2,Q(:,4))+1;         % a2
Qtrans(:,5) = icdf(PD_a3,Q(:,5));           % a3
Qtrans(:,6) = 0.1 + 4*Q(:,6);               % k1
Qtrans(:,7) = 0.5 + 10*Q(:,7);              % k2
Qtrans(:,8) = 0.0005 + 0.002*Q(:,8);        % c2

cd('J:\Documents\MATLABroutines\EarthquakeSimulation\');
for i = 1:length(Q)
    u = randn(N,1);
    wmid = Qtrans(i,1);
    wdot = Qtrans(i,2);
    zeta = Qtrans(i,3);
    a2 = Qtrans(i,4);
    a3 = Qtrans(i,5);
    a1 = COEF(1);
    for j = 2:size(Regressors,1);
        a1 = a1 + COEF(j)*eval(Regressors(j,:));
    end
    a1 = exp(a1);
    q = modulfun([a1 a2 a3]',N,Ts);
    [~,~,tmid(i)] = EQarias(q,Ts);
    z(:,i) = syntheticEQirf([wmid wdot zeta a1 a2 a3],u,Ts,tmid(i),0.25);
end
%%
cd('J:\Documents\MATLABroutines\PCE\')

% Time index
Tspan = (0:Ts:(N-1)*Ts);
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
IC = zeros(4,1);

for i = 1:K
    % Model properties
    Sys.k1 = Qtrans(i,6);
    Sys.k2 = Qtrans(i,7);
    Sys.c2 = Qtrans(i,8);

    [~,Xsim] = ode45(@(t,x) duffingTDOF(t,x,Tspan,z(:,i),Sys),Tspan,IC,options);
    for t=1:N
        Xsimdot(t,:) = duffingTDOF(Tspan(t),Xsim(t,:),Tspan(t),z(t,i),Sys);
    end
    acc1(:,i) = Xsimdot(:,2);
    acc2(:,i) = Xsimdot(:,4);
end
mscohere(z(:,10),acc2(:,10))

y = zeros(N,K);
x = zeros(N,K);
for i = 1:K
    y(:,i) = acc2(:,i);
    x(:,i) = z(:,i);
    figure(i)
    subplot(211),pwelch(acc2(:,i));
    subplot(212),mscohere(z(:,i),acc2(:,i))
end





%% Filter construction
forder = cheb2ord(0.2, 0.25, 0.01, 60);
[Z,P,G] = cheby2(forder,60,0.25);       % 4 Hz
[sos,g] = zp2sos(Z,P,G);