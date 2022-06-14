%% Transient analysis through ANSYS (sinusoidal excitation) 
clear,clc,close all
cd('J:\Documents\MATLABroutines\PCNARXjournal\SimulationModels\')

% Number of samples
N = 5000;
% Sampling period
Ts = 0.5;
% Time index
Tspan = (0:Ts:(N-1)*Ts);

% Filter construction
% [z,p,k] = cheby2(20,40,0.475);
% [sos,g] = zp2sos(z,p,k);
% fsig = filtfilt(sos,g,sig);

options = odeset('RelTol',1e-3,'AbsTol',1e-6);
IC = zeros(4,1);
power = -3:2;

figure(1),
clr = colormap(lines(5));
hold on
figure(2),
clr = colormap(lines(5));
hold on
figure(3),
clr = colormap(lines(5));
hold on
figure(4),
clr = colormap(lines(5));
hold on
for i = 1:5
    disp(i)
    % Model properties
    Sys.m1 = 1;         Sys.m2 = 0.25;
	Sys.k1 = 0.75;      Sys.k2 = 0.25;
    Sys.knl1 = 0.5;     Sys.knl2 = 0.05;    
    U1 = zeros(1,N);
    U2 = (10^(power(i)))*rand(1,N);
    % U2 = linspace(0,10^(power(i)),N).*sin(linspace(-10*pi,10*pi,N));
    [Tsim,Ysim] = ode45(@(t,x) cubicstiffnessTDOF(t,x,Tspan,[U1;U2],Sys),Tspan,IC,options); 
    for t=1:N
        Ysimdot(t,:) = cubicstiffnessTDOF(Tspan(t),Ysim(t,:),Tspan(t),[0;U2(t)],Sys);
    end
    dspl(:,i) = Ysimdot(1001:end,4);
    force(:,i) = U2(1001:end);
    [C,F] = mscohere(force(:,i),dspl(:,i),512,400,512,1/Ts);
    figure(1)
    plot(dspl(:,i),force(:,i),'Color',clr(i,:))
    [Pxx,F] = pwelch(dspl(:,i),512,400,512,1/Ts);
    [Txy,F] = tfestimate(force(:,i),dspl(:,i),512,400,512,1/Ts);
    figure(2)
    plot(F,20*log10(abs(Pxx)),'Color',clr(i,:))
    figure(3)
    plot(F,20*log10(abs(Txy)),'Color',clr(i,:))
    figure(4)
    plot(F,C,'Color',clr(i,:))
end

A = [0 1 0 0; -(Sys.k1+Sys.k2)/Sys.m1 0 Sys.k2/Sys.m1 0;...
    0 0 0 1 ; Sys.k2/Sys.m2 0 -Sys.k2/Sys.m2 0;];
wn = abs(eig(A));
zeta = -real(eig(A))./wn;
wn = wn/2/pi;


%% LHS experiments (linear system k - c)
clear,clc,close all
if ispc
    cd('J:\Documents\MATLABroutines\EarthquakeSimulation\');
    load('J:\Documents\MATLABroutines\EarthquakeSimulation\FittedDist')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/EarthquakeSimulation/');
    load('/media/minas/USB DISK/Documents/MATLABroutines/EarthquakeSimulation/FittedDist/')
end

% Number of experiments
K = 10; 
% Number of samples
N = 5000;
% Number of variables (beta, gamma, n)
M = 3;
% Sampling period
Ts = 1;

% Latin Hypercube Sampling
Q = -1+2*lhsdesign(K,M,'criterion','maximin','iterations',10);
Qtrans(:,1) = 0.5*(0.5*Q(:,1)+0.5);                  % k
Qtrans(:,2) = 0.01 + 0.09*(0.5*Q(:,2)+0.5);            % c
Qtrans(:,3) = 10.^(-2+4*(0.5*Q(:,3)+0.5));             % sigma
% System constant properties
Sys.m = 1;
Sys.knl = 0;

if ispc
    cd('J:\Documents\MATLABroutines\PCE\')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/PCE/')
end
% Time index
Tspan = (0:Ts:(N-1)*Ts);
options = odeset('RelTol',1e-3,'AbsTol',1e-6);
IC = zeros(2,1);

for i = 1:K
    % Model properties
    disp(i)
    F(:,i) = randn(N,1);
    Sys.k = Qtrans(i,1);     
    Sys.c = Qtrans(i,2);
    F(:,i) = Qtrans(i,3)*randn(N,1);
    [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F(:,i)',Sys),Tspan,IC,options); 
    dspl(:,i) = Xsim(:,1);
%    fnum = fir1(10,0.5); fden = 1;
%    dspl(:,i) = filtfilt(fnum,fden,dspl(:,i));
    [z,p,k] = cheby2(20,40,0.5);
    [sos,g] = zp2sos(z,p,k);     
    dspl(:,i) = filtfilt(sos,g,dspl(:,i));
    fdspl(:,i) = dspl(501:2:end-500,i);
    % F(:,i) = filtfilt(fnum,1,F(:,i));
    F(:,i) = filtfilt(sos,g,F(:,i));
    fF(:,i) = F(501:2:end-500,i);
end
N = length(fdspl);

y = zeros(N,K);
x = zeros(N,K);
for i = 1:K
    y(:,i) = fdspl(:,i);
    x(:,i) = fF(:,i);
    figure(i)
    subplot(211),pwelch(fdspl(:,i));
    subplot(212),mscohere(fF(:,i),fdspl(:,i))
end

%% LHS experiments (input gain only)
clear,clc,close all
if ispc
    cd('J:\Documents\MATLABroutines\EarthquakeSimulation\');
    load('J:\Documents\MATLABroutines\EarthquakeSimulation\FittedDist')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/EarthquakeSimulation/');
    load('/media/minas/USB DISK/Documents/MATLABroutines/EarthquakeSimulation/FittedDist/')
end

% Number of experiments
K = 10; 
% Number of samples
N = 5000;
% Number of variables (beta, gamma, n)
M = 2;
% Sampling period
Ts = 1;
%%

% Latin Hypercube Sampling
Q = -1+2*lhsdesign(K,M,'criterion','maximin','iterations',10);
Qtrans(:,1) = 10.^(-1+(0.5*Q(:,1)+0.5));            % sigma
Qtrans(:,2) = 0.5 + (0.5*Q(:,2)+0.5);               % k
% System constant properties
Sys.m = 4;
Sys.knl = 0.2;
Sys.c = 0.05;

if ispc
    cd('J:\Documents\MATLABroutines\PCE\')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/PCE/');
end

% Time index
Tspan = (0:Ts:(N-1)*Ts);
options = odeset('RelTol',1e-3,'AbsTol',1e-6);
IC = zeros(2,1);

for i = 1:K
    % Model properties
    disp(i)
    F(:,i) = Qtrans(i,1)*randn(N,1);
    Sys.k = Qtrans(i,2);     
    [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F(:,i)',Sys),Tspan,IC,options); 
    dspl(:,i) = Xsim(:,1);
    [z,p,k] = cheby2(20,40,0.5);
    [sos,g] = zp2sos(z,p,k);     
    dspl(:,i) = filtfilt(sos,g,dspl(:,i));
    fdspl(:,i) = dspl(501:2:end-500,i);
    F(:,i) = filtfilt(sos,g,F(:,i));
    fF(:,i) = F(501:2:end-500,i);
end

N = length(fdspl);
y = zeros(N,K);
x = zeros(N,K);
for i = 1:K
    y(:,i) = fdspl(:,i);
    x(:,i) = fF(:,i);
    figure(i)
    subplot(211),pwelch(fdspl(:,i));
    subplot(212),mscohere(fF(:,i),fdspl(:,i))
end




%% LHS experiments (modulating function only)
clear,clc,close all
if ispc
    cd('J:\Documents\MATLABroutines\EarthquakeSimulation\');
    load('J:\Documents\MATLABroutines\EarthquakeSimulation\FittedDist')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/EarthquakeSimulation/');
    load('/media/minas/USB DISK/Documents/MATLABroutines/EarthquakeSimulation/FittedDist/')
end

% Number of experiments
K = 10; 
% Number of samples
N = 5000;
% Number of variables (beta, gamma, n)
M = 2;
% Sampling period
Ts = 1;

% Latin Hypercube Sampling
Q = -1+2*lhsdesign(K,M,'criterion','maximin','iterations',10);
Qtrans(:,1) = 10.^(-1+(0.5*Q(:,1)+0.5));            % sigma
Qtrans(:,2) = 0.5 + (0.5*Q(:,2)+0.5);               % k
% System constant properties
Sys.m = 4;
Sys.knl = 0.2;
Sys.c = 0.05;

if ispc
    cd('J:\Documents\MATLABroutines\PCE\')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/PCE/');
end

% Time index
Tspan = (0:Ts:(N-1)*Ts);
options = odeset('RelTol',1e-3,'AbsTol',1e-6);
IC = zeros(2,1);

for i = 1:K
    % Model properties
    disp(i)
    F(:,i) = Qtrans(i,1)*randn(N,1);
    Sys.k = Qtrans(i,2);     
    [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F(:,i)',Sys),Tspan,IC,options); 
    dspl(:,i) = Xsim(:,1);
    [z,p,k] = cheby2(20,40,0.5);
    [sos,g] = zp2sos(z,p,k);     
    dspl(:,i) = filtfilt(sos,g,dspl(:,i));
    fdspl(:,i) = dspl(501:2:end-500,i);
    F(:,i) = filtfilt(sos,g,F(:,i));
    fF(:,i) = F(501:2:end-500,i);
end

N = length(fdspl);
y = zeros(N,K);
x = zeros(N,K);
for i = 1:K
    y(:,i) = fdspl(:,i);
    x(:,i) = fF(:,i);
    figure(i)
    subplot(211),pwelch(fdspl(:,i));
    subplot(212),mscohere(fF(:,i),fdspl(:,i))
end



%% LHS experiments (IRF only)
clear,clc,close all
if ispc
    cd('J:\Documents\MATLABroutines\EarthquakeSimulation\');
    load('J:\Documents\MATLABroutines\EarthquakeSimulation\FittedDist')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/EarthquakeSimulation/');
    load('/media/minas/USB DISK/Documents/MATLABroutines/EarthquakeSimulation/FittedDist/')
end

% Number of experiments
K = 100; 
% Number of samples
N = 5000;
% Number of variables (beta, gamma, n)
M = 3;
% Sampling period
Ts = 0.005;

% Latin Hypercube Sampling
Q = -1+2*lhsdesign(K,M,'criterion','maximin','iterations',10);
Qtrans(:,1) = icdf(PD_wmid,0.5*(1+Q(:,1)))*2*pi;
Qtrans(:,2) = icdf(PD_wdot,0.5*(1+Q(:,2)))*2*pi;
Qtrans(:,3) = icdf(PD_zeta,0.5*(1+Q(:,3)));

% System constant properties
Sys.m = 1;
Sys.knl = 10000*0.2;
Sys.c = 0.05;
Sys.k = 10000*0.5;

if ispc
    cd('J:\Documents\MATLABroutines\PCE\')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/PCE/');
end

% Time index
Tspan = (0:Ts:(N-1)*Ts);
options = odeset('RelTol',1e-3,'AbsTol',1e-6);
IC = zeros(2,1);



for i = 1:K
    cd('J:\Documents\MATLABroutines\EarthquakeSimulation\')
    disp(i)
    % Model properties
    a1 = 1e14;
    a2 = 12;
    a3 = 33;
    q = modulfun([a1 a2 a3]',N,Ts);
    [~,~,tmid(i)] = EQarias(q,Ts);
    wmid = Qtrans(i,1);
    wdot = Qtrans(i,2);
    zeta = Qtrans(i,3);
    % tmid(i) = gaminv(0.45,2*a3,2*a2-1)*N;
    u = randn(N,1);
    F(:,i) = syntheticEQirf([wmid wdot zeta a1 a2 a3],u,Ts,tmid(i),0.25);
    cd('J:\Documents\MATLABroutines\PCE\')
         
    [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessSDOF(t,x,Tspan,F(:,i)',Sys),Tspan,IC,options); 
    dspl(:,i) = Xsim(:,1);
    [z,p,k] = cheby2(20,40,0.5);
    [sos,g] = zp2sos(z,p,k);     
    dspl(:,i) = filtfilt(sos,g,dspl(:,i));
    fdspl(:,i) = dspl(501:2:end-500,i);
    F(:,i) = filtfilt(sos,g,F(:,i));
    fF(:,i) = F(501:2:end-500,i);
end

N = length(fdspl);
y = zeros(N,K);
x = zeros(N,K);
for i = 1:10
    y(:,i) = fdspl(:,i);
    x(:,i) = fF(:,i);
    figure(i)
    subplot(211),pwelch(fdspl(:,i));
    subplot(212),mscohere(fF(:,i),fdspl(:,i))
end



%% Conventional ARX - NARX identification
clc, clear *rss*
disp('Linear    Nonlinear')
disp('===================')
for ii=1:10
    IdData = iddata(fdspl(:,ii),fF(:,ii),Ts);
    na_max = 2;
    for na = 2:na_max
        % disp(['AR/X order: ',num2str(na)]);
        model = arx(IdData,[na na 0]);
        res = resid(model,IdData);
        rss_sss(na) = 100*norm(res.outputdata)^2/norm(fdspl(:,ii))^2;
        sigma(na) = var(res.outputdata);
        bic(na) = log(sigma(na)) + (2*na+1)*log(N)/N;
        NLmodel = idnlarx([na na 0], linear, 'InputName', 'u1', 'OutputName', 'y1','Focus','Prediction',...
            'CustomReg',{'y1(t-1)','y1(t-2)',...
            'u1(t)','u1(t-1)','u1(t-2)',...
            'y1(t-1)^2','y1(t-2)^2',...
            'y1(t-1)^3','y1(t-2)^3',...
            'y1(t-1)^4','y1(t-2)^4',...
            'y1(t-1)*y1(t-2)',...
            'y1(t-1)*u1(t-1)','y1(t-2)*u1(t-2)',...
            'y1(t-1)^2*u1(t-1)^2','y1(t-2)^2*u1(t-2)^2',...
            'y1(t-1)^3*u1(t-1)^3','y1(t-2)^3*u1(t-2)^3'});
        % Estimate the model parameters a1, a2, ... a6.
        nlmodel = nlarx(IdData, NLmodel);
%         nlmodel = nlarx(IdData,[na na 0],'sigmoidnet',...
%                   'NonlinearRegressors','search');
        %     nlmodel = nlarx(IdData,[na na 0],'wavenet',...
        %            'NonlinearRegressors','search');
        % nlmodel = nlarx(IdData,[0 0 0],'linear','CustomReg',{'y1(t-1)','y1(t-2)',...
        %     'u1(t)','u1(t-1)','u1(t-2)'});
%             'y1(t-1)^2','y1(t-2)^2',...
%             'y1(t-1)^3','y1(t-2)^3',...
%             'y1(t-1)^4','y1(t-2)^4',...
%             'y1(t-1)*y1(t-2)',...
%             'y1(t-1)*u1(t-1)','y1(t-2)*u1(t-2)',...
%             'y1(t-1)^2*u1(t-1)^2','y1(t-2)^2*u1(t-2)^2',...
%             'y1(t-1)^3*u1(t-1)^3','y1(t-2)^3*u1(t-2)^3',...
%             'u1(t)','u1(t-1)','u1(t-2)'});
        % NLres = resid(nlmodel,IdData);
        Yhat = predict(nlmodel,IdData,1);
        NLres = acc(:,ii) - Yhat.outputdata;
        NLrss_sss(na) = 100*norm(NLres)^2/norm(acc(:,ii))^2;
        NLsigma(na) = var(NLres);
        NLbic(na) = log(NLsigma(na)) + (2*na+1)*log(N)/N;
        disp([num2str(rss_sss(na)),'   ',num2str(NLrss_sss(na))])
    end
end


%% Conventional ARX - NARX identification (my routine)
if ispc
    cd('J:\Documents\MATLABroutines\PCE\m-files')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/PCE/m-files/');
end


clc, clear *rss*
options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';

P = 3;
q = 1;
maxterms = [];

disp('Linear    Nonlinear')
disp('===================')
for ii=1:100
    Y = fdspl(:,ii);
    X = fF(:,ii);
    na_max = 2;
    for na = 2:na_max
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
    if ii==10
        % Cross-corellation(z,z^2)
        ccf(res-mean(res),(res-mean(res)).^2,100,0.8);
        ccf(NLres-mean(NLres),(NLres-mean(NLres)).^2,100,0.8);
    end
end




%% GA - NARX
cd('J:\Documents\MATLABroutines\PCE\m-files')
clc, clear *rss*
options.focus = 'prediction';
options.GA.PopulationSize = 100 ;
options.GA.PlotFcns = {};
options.GA.Display = 'off';
P = 3;
q = 1;
maxterms = 2;

disp('Linear    Nonlinear')
disp('===================')
for ii=1:10
    Y = fdspl(:,ii);
    X = fF(:,ii);
    na = 2;
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
figure(1)
ccf(res-mean(res),(res-mean(res)).^2,100,0.8);
figure(2)
ccf(NLres-mean(NLres),(NLres-mean(NLres)).^2,100,0.8);


%% PC-ARX modelling
cd('J:\Documents\MATLABroutines\PCE\m-files')
options.maxsize = 1e9;
options.basis = 'legen';
options.focus = 'prediction';
% AR modelling
% Q = Qtrans';
n_max = 6;
q = 1;
for na = 1:n_max
    for p = 0:5
        disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
        indx{1} = combinations(M,p,q);
        indx{2} = indx{1};
        B = length(indx{1});
        [aij,bij,res,criteria] = pcarx(y,x,Q',[na 0 na],indx,[],options);
        rss_sss(na,p+1) = criteria.rss_sss;
        connum(na,p+1) = criteria.connum;
        spp(na,p+1) = criteria.spp;
        sigma = var(res);
        bic(na,p+1) = log(sigma) + ((na+na+1)*B)*log(K*N)/(K*N);
    end
end


%% PC-NARX modelling
if ispc
    cd('J:\Documents\MATLABroutines\PCE\m-files')
else
    cd('/media/minas/USB DISK/Documents/MATLABroutines/PCE/m-files/');
end

options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';
% AR modelling
% Q = Qtrans';
n_max = 4;
Nr = length(regstr);
q = 1;
for na = 1:n_max
    for p = 0:3
        disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
        indx{1} = combinations(M,p,q);
        indx{2} = indx{1};
        B = length(indx{1});
        [~,~,~,InitRegstr] = narx(y(:,1),x(:,1),[na 0 na],[],3,1,2,[],options);
        [aij,res,criteria] = pcnarx(y,x,Q',[na 0 na],indx,InitRegstr,[],options);
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

