%% Transient analysis through ANSYS (sinusoidal excitation) 
clear,clc,close all
cd('J:\Documents\MATLABroutines\PCE\')

% Number of samples
N = 1000;
% Sampling period
Ts = 0.5;
% Time index
Tspan = (0:Ts:(N-1)*Ts);

options = odeset('RelTol',1e-3,'AbsTol',1e-6);
IC = zeros(4,1);
power = -4:1;

figure(1),
clr = colormap(lines(10));
hold on
for i = 1:5
    disp(i)
    % Model properties
    Sys.m1 = 1;         Sys.m2 = 0.05;
    Sys.k1 = 1;         Sys.k2 = 1;
    Sys.knl1 = 0.5;     Sys.knl2 = 0.1;
    U(1,:) = zeros(1,N);
    U(2,:) = (10^(power(i)))*sin(Tspan*4*pi/(N-1));
    [Tsim,Ysim] = ode45(@(t,x) cubicstiffnessTDOF(t,x,Tspan,U,Sys),Tspan,IC,options); 
    figure(1)
    plot(Ysim(:,3),U(2,:),'Color',clr(i,:))
    figure(i+1)
    plot(20*log10(abs(fft(Ysim(:,3)))))
end

A = [0 1 0 0; -(Sys.k1+Sys.k2)/Sys.m1 0 Sys.k2/Sys.m1 0; 0 0 0 1; Sys.k2/Sys.m2 0 -Sys.k2/Sys.m2 0];
wn = abs(eig(A));
zeta = -real(eig(A))./wn;
wn = wn/2/pi;

%% LHS experiments (input gain only)
clear,clc,close all
cd('J:\Documents\MATLABroutines\EarthquakeSimulation\');
load('J:\Documents\MATLABroutines\EarthquakeSimulation\FittedDist')

% Number of experiments
K = 10; 
% Number of samples
N = 1000;
% Number of variables (beta, gamma, n)
M = 3;
% Sampling period
Ts = 0.4;

% Latin Hypercube Sampling
Q = lhsdesign(K,M,'criterion','maximin','iterations',10);
Qtrans(:,1) = 10.^(-5+5*Q(:,1));         % sigma
Qtrans(:,2) = 1 + 0.9*Q(:,2);            % k1
Qtrans(:,3) = 1 + 0.9*Q(:,3);            % knl1

% System constant properties
Sys.m1 = 1;         Sys.m2 = 0.05;
Sys.k2 = 0.5;       Sys.knl2 = 0.1;

cd('J:\Documents\MATLABroutines\PCE\')
% Time index
Tspan = (0:Ts:(N-1)*Ts);
options = odeset('RelTol',1e-3,'AbsTol',1e-6);
IC = zeros(4,1);

for i = 1:K
    % Model properties
    disp(i)
    F1 = zeros(1,N);
    F2(:,i) = Qtrans(i,1)*randn(N,1);
    Sys.k1 = Qtrans(i,2);     Sys.knl1 = Qtrans(i,3);
    [Tsim,Xsim] = ode45(@(t,x) cubicstiffnessTDOF(t,x,Tspan,[F1;F2(:,i)'],Sys),Tspan,IC,options); 
    for t=1:N
        Xsimdot(t,:) = cubicstiffnessTDOF(Tspan(t),Xsim(t,:),Tspan(t),[0 ;F2(t,i)],Sys);
    end
    acc1(:,i) = Xsimdot(:,2);
    acc2(:,i) = Xsimdot(:,4);
end

y = zeros(N,K);
x = zeros(N,K);
for i = 1:K
    y(:,i) = acc1(:,i);
    x(:,i) = F2(:,i);
    figure(i)
    subplot(211),pwelch(acc1(:,i));
    subplot(212),mscohere(F2(:,i),acc1(:,i))
end


%% Conventional ARX - NARX identification
clc, clear *rss*
for ii=1:10
    IdData = iddata(acc1(:,ii),F2(:,ii),Ts);
    na_max = 4;
    AR = zeros(na_max);
    X = zeros(na_max+1);
    for na = 4:na_max
        disp(['AR/X order: ',num2str(na)]);
        model = arx(IdData,[na na 0]);
        AR(na,1:na)= model.a(2:end);
        X(na,1:na)= model.b(1:end);
        res = resid(model,IdData);
        rss_sss(na) = 100*norm(res.outputdata)^2/norm(acc1(:,1))^2
        sigma(na) = var(res.outputdata);
        bic(na) = log(sigma(na)) + (2*na+1)*log(N)/N;
        %     nlmodel = nlarx(IdData,[na na 0],'sigmoidnet',...
        %           'NonlinearRegressors','search');
        %     nlmodel = nlarx(IdData,[na na 0],'wavenet',...
        %            'NonlinearRegressors','search');
        nlmodel = nlarx(IdData,[0 0 0],'linear','CustomReg',{'y1(t-1)','y1(t-2)','y1(t-3)','y1(t-4)',...
            'y1(t-1)^2','y1(t-2)^2','y1(t-3)^2','y1(t-4)^2',...
            'y1(t-1)*u1(t-1)','y1(t-2)*u1(t-2)','y1(t-3)*u1(t-3)','y1(t-4)*u1(t-4)',...
            'y1(t-1)^3','y1(t-2)^3','y1(t-3)^3','y1(t-4)^3',...
            'u1(t)','u1(t-1)','u1(t-2)','u1(t-3)','u1(t-4)'});
        NLres = resid(nlmodel,IdData);
        NLrss_sss(na) = 100*norm(NLres.outputdata)^2/norm(acc1(:,1))^2
        NLsigma(na) = var(NLres.outputdata);
        NLbic(na) = log(NLsigma(na)) + (2*na+1)*log(N)/N;
    end
end

%% PC-ARX modelling
cd('J:\Documents\MATLABroutines\PCE\m-files')

options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';
% AR modelling
Q = Qtrans';
n_max = 8;
q = 1;
for na = 2:n_max
    for p = 1:4
        disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
        indx{1} = combinations(M,p,q);
        indx{2} = indx{1};
        B = length(indx{1});
        [aij,bij,res,criteria] = pcarx(y,x,Q,[na 0 na],indx,[],options);
        rss_sss(na,p) = criteria.rss_sss;
        connum(na,p) = criteria.connum;
        spp(na,p) = criteria.spp;
        sigma = var(res);
        bic(na,p) = log(sigma) + ((na+na+1)*B)*log(K*N)/(K*N);
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
Q = Qtrans';
q = 1;
p = 3;
indx{1} = combinations(M,p,q);
indx{2} = indx{1};
B = length(indx{1});
na = 12;
% [aij,bij,res,criteria] = pcarx(y,x,Q,[na 0 na],indx,[],options);
[ms,criteria] = folsPCARX(y,x,Q,[na 0 na],indx,100,[],options);

% [INDX,aij,bij,res,criteria,GAoutput] = gaPCARX(y,x,Q,[na 0 na],3,1,options);

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

%%
cd('J:\Documents\MATLABroutines\PCE\m-files')
options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';
% AR modelling
Q = Qtrans';
n_max = 4;
q = 1;
for na = 2:n_max
    for p = 1:4
        disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
        indx{1} = combinations(M,p,q);
        indx{2} = indx{1};
        B = length(indx{1});
        [aij,bij,res,criteria] = pcarx(y,x,Q,[na 0 na],indx,[],options);
        rss_sss(na,p) = criteria.rss_sss;
        connum(na,p) = criteria.connum;
        spp(na,p) = criteria.spp;
        sigma = var(res);
        bic(na,p) = log(sigma) + ((na+na+1)*B)*log(K*N)/(K*N);
    end
end

%% Conventional ARX identification
IdData = iddata(acc2(:,1),z(:,1),Ts);
na_max = 4;
AR = zeros(na_max);
X = zeros(na_max+1);
for na = 2:na_max
    disp(['AR/X order: ',num2str(na)]);
    model = arx(IdData,[na na 0]);
    AR(na,1:na)= model.a(2:end); 
    X(na,1:na)= model.b(1:end); 
    res = resid(model,IdData);
    rss_sss(na) = 100*norm(res.outputdata)^2/norm(a2)^2;
    sigma(na) = var(res.outputdata);
    bic(na) = log(sigma(na)) + (2*na+1)*log(N)/N;
    
    nlmodel = nlarx(IdData,[na na 0],'wavenet',...
          'NonlinearRegressors','search');
    NLres = resid(nlmodel,IdData);
    NLrss_sss(na) = 100*norm(NLres.outputdata)^2/norm(a2)^2;
    NLsigma(na) = var(NLres.outputdata);
    NLbic(na) = log(NLsigma(na)) + (2*na+1)*log(N)/N;
end


%% Time histories plot (Node 12 acceleration response)
linemark = [' -';'--';'-.'];

figure(2)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
hold on
clr = colormap(lines(3));

for i = 1:3
    if i==1
        figure(1)
        set(gcf,'PaperOrientation','portrait','papertype','A4',...
            'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
        subplot(5,1,1:3),plot(Tspan,elcentro,'linewidth',1)
        set(gca,'xlim',[1/64 1000/64],'ylim',[-50 50],'FontName','TimesNewRoman','Fontsize',11)
        ylabel('$\alpha $ (m/s$^2$)','Interpreter','Latex','Fontsize',16)
        xlabel('Time (s)','FontName','TimesNewRoman','Fontsize',16)
        grid on
    end
    
    figure(2)
    subplot(7,1,(i-1)*2+1:i*2)
    plot(Tspan,y2dot(:,i),linemark(1,:),'linewidth',1,'color',clr(i,:))
    set(gca,'xlim',[1/64 1000/64],'ylim',[-50 50],'FontName','TimesNewRoman','Fontsize',12)
    if i == 3
        xlabel('Time (s)','FontName','TimesNewRoman')
    else
        set(gca,'xticklabel',[])
    end
    ylabel(['$y_',num2str(i),'[t]$ (m/s$^2$)'],'Interpreter','Latex','Fontsize',12)
    grid on
    box on
    
    P(:,i) = pwelch(y1dot(:,i),512,480,512,1/Ts);
    Txy(:,i) = tfestimate(elcentro,y1dot(:,i),512,480,512,1/Ts);
    [Cxy(:,i),F] = mscohere(elcentro,y1dot(:,i),512,480,512,1/Ts);
end


figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(5,1,1:2),plot(F,20*log10(abs(Txy(:,1))),linemark(1,:),'linewidth',1,'Color',clr(1,:)),hold on
plot(F,20*log10(abs(Txy(:,2))),linemark(2,:),'linewidth',1,'Color',clr(2,:))
plot(F,20*log10(abs(Txy(:,3))),linemark(3,:),'linewidth',1,'Color',clr(3,:))
axis([0 1/(2*Ts) -60 20])
grid on
set(gca,'xticklabel',[],'FontName','TimesNewRoman','Fontsize',11)
ylabel('FRF magnitude (dB)','FontName','TimesNewRoman','Fontsize',12)
legend({'$\xi_1$','$\xi_2$','$\xi_3$'},'Interpreter','Latex','Fontsize',10,'Location','Southwest')
subplot(5,1,3:4),plot(F,Cxy(:,1),linemark(1,:),'linewidth',1,'Color',clr(1,:)),hold on
plot(F,Cxy(:,2),linemark(2,:),'linewidth',1,'Color',clr(2,:))
plot(F,Cxy(:,3),linemark(3,:),'linewidth',1,'Color',clr(3,:))
grid on
axis([0 1/(2*Ts) 0 1])
set(gca,'FontName','TimesNewRoman','Fontsize',11)
xlabel('Frequency (Hz)','FontName','TimesNewRoman','Fontsize',12)
ylabel('Coherence function','FontName','TimesNewRoman','Fontsize',12)
legend({'$\xi_1$','$\xi_2$','$\xi_3$'},'Interpreter','Latex','Fontsize',10,'Location','Southwest')


%% Transient analysis through ANSYS (Pacoima Dam accelerogram) 
clear,clc,close all
cd('C:\Users\sminas\Desktop\PCE\ANSYS\')

% Number of experiments
K = 1; 
% Number of samples
N = 1000;
% Number of variables (Ex, Density, Prxy, Cross-section area)
M = 3;

% Latin Hypercube Sampling
Qtrans = [1.95e11; 0.045; 250e6];
save('BeamProps.txt','Qtrans','-ascii')

% Run Ansys
delete('file.*')
dos(' "C:\Program Files\ANSYS Inc\v130\ansys\bin\winx64\ANSYS130.exe" -b -i "C:\Users\sminas\Desktop\PCE\ANSYS\ShearFrameBilinearSeismic.txt" -o "output.txt"');

% Plot signal realization
load('C:\Users\sminas\Desktop\PCE\ANSYS\ShearFrameSeismic\SFresPACD.txt')
load('C:\Users\sminas\Desktop\PCE\ANSYS\ShearFrameSeismic\pacoima.txt')
figure,
subplot(311),plot(SFresPACD(:,1),pacoima(1:N))
set(gca,'Xticklabel',[])
ylabel('Input acc. (m/s^2)')
subplot(312),plot(SFresPACD(:,1),SFresPACD(:,3:4))
legend({'Node 1','Node 3'})
ylabel('Reaction forces (K)')
set(gca,'Xticklabel',[])
subplot(313),plot(SFresPACD(:,1),SFresPACD(:,5))
ylabel('Top floor displacement (m)')
xlabel('Time (s)')

figure,
plot3(SFresPACD(:,5),-(SFresPACD(:,3)+SFresPACD(:,4)),linspace(0,100,N),'-o')
ylabel('Total force (K)')
xlabel('Top floor displacement (m)')
zlabel('Max input acceleration (m/s^2)')
axis tight
view([-25,5])
grid on


%% Time histories plot (Node 12 acceleration response)
N = 1000;
K = 100;
Ts = 0.01;
linemark = [' -';'--';'-.'];

T0 = 1;
Tf = 250;


figure(2)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
hold on
clr = colormap(lines(3));

for i = 1:3
    out = load(['ShearFrameSeismic\SFres',num2str(i),'.txt']);
    if i==1
        figure(1)
        set(gcf,'PaperOrientation','portrait','papertype','A4',...
            'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
        subplot(5,1,1:3),plot(out(:,1),out(:,2),'linewidth',1)
        set(gca,'xlim',[1/64 1000/64],'ylim',[-20 20],'FontName','TimesNewRoman','Fontsize',11)
        ylabel('$\alpha $ (m/s$^2$)','Interpreter','Latex','Fontsize',16)
        xlabel('Time (s)','FontName','TimesNewRoman','Fontsize',16)
        grid on
        if print_to_eps == 'y'
            print(pictype,resolution,[write_dir,'elcentro']),close;
            result = eps2xxx([write_dir,'elcentro.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
        end
    end
    
    figure(2)
    subplot(7,1,(i-1)*2+1:i*2)
    % subplot(3,3,(i-1)*3+1:i*3-1)
    plot(out(:,1),out(:,15),linemark(1,:),'linewidth',1,'color',clr(i,:))
    set(gca,'xlim',[1/64 1000/64],'ylim',[-50 50],'FontName','TimesNewRoman','Fontsize',12)
    if i == 3
        xlabel('Time (s)','FontName','TimesNewRoman')
    else
        set(gca,'xticklabel',[])
    end
    ylabel(['$y_',num2str(i),'[t]$ (m/s$^2$)'],'Interpreter','Latex','Fontsize',12)
    grid on
    box on
    
    P(:,i) = pwelch(out(:,15),512,480,512,1/Ts);
    Txy(:,i) = tfestimate(out(:,2),out(:,15),512,480,512,1/Ts);
    [Cxy(:,i),F] = mscohere(out(:,2),out(:,15),512,480,512,1/Ts);
end
if print_to_eps == 'y'
    print(pictype,resolution,[write_dir,'signals']),close;
    result = eps2xxx([write_dir,'signals.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end

figure(3)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(5,1,1:2),plot(F,20*log10(abs(Txy(:,1))),linemark(1,:),'linewidth',1,'Color',clr(1,:)),hold on
plot(F,20*log10(abs(Txy(:,2))),linemark(2,:),'linewidth',1,'Color',clr(2,:))
plot(F,20*log10(abs(Txy(:,3))),linemark(3,:),'linewidth',1,'Color',clr(3,:))
axis([0 1/(2*Ts) -60 20])
grid on
set(gca,'xticklabel',[],'FontName','TimesNewRoman','Fontsize',11)
ylabel('FRF magnitude (dB)','FontName','TimesNewRoman','Fontsize',12)
legend({'$\xi_1$','$\xi_2$','$\xi_3$'},'Interpreter','Latex','Fontsize',10,'Location','Southwest')
subplot(5,1,3:4),plot(F,Cxy(:,1),linemark(1,:),'linewidth',1,'Color',clr(1,:)),hold on
plot(F,Cxy(:,2),linemark(2,:),'linewidth',1,'Color',clr(2,:))
plot(F,Cxy(:,3),linemark(3,:),'linewidth',1,'Color',clr(3,:))
grid on
axis([0 1/(2*Ts) 0 1])
set(gca,'FontName','TimesNewRoman','Fontsize',11)
xlabel('Frequency (Hz)','FontName','TimesNewRoman','Fontsize',12)
ylabel('Coherence function','FontName','TimesNewRoman','Fontsize',12)
legend({'$\xi_1$','$\xi_2$','$\xi_3$'},'Interpreter','Latex','Fontsize',10,'Location','Southwest')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'preproc']),close;
     result = eps2xxx([write_dir,'preproc.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end

%% Hysteresis loop
clear,clc,close all
local_dir = 'C:\Users\sminas\Desktop\PCE\';
cd([local_dir,'ANSYS\'])

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

N = 1000;
K = 20;
T0 = 1;
Tf = 250;
out = load('ShearFrameSeismic\SFres1.txt');

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(9,3,[2:3 5:6 8:9]),plot(out(:,1),out(:,2),'linewidth',1)
set(gca,'xlim',[1/64 1000/64],'ylim',[-20 20],'FontName','TimesNewRoman','Fontsize',10)
ylabel('$\alpha $ (m/s$^2$)','Interpreter','Latex','Fontsize',12)
xlabel('Time (s)','FontName','TimesNewRoman','Fontsize',12)
grid on
subplot(9,3,[14:15 17:18 20:21 23:24 26:27]),plot(out(T0:Tf,5),-(out(T0:Tf,3)+out(T0:Tf,4)),'linewidth',1.2)
ylabel('Shear stress (Pa)','Fontname','TimesNewRoman','Fontsize',12)
xlabel('Top floor displacement (m)','Fontname','TimesNewRoman','Fontsize',12)
axis tight
grid on
set(gca,'Fontsize',10,'xlim',[-0.3,0.3],'ylim',[-9e5 9e5],...
    'Xtick',-0.25:0.1:0.25,'Fontname','TimesNewRoman')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'hysteresis']),close;
     result = eps2xxx([write_dir,'hysteresis.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end



%% Material properties plot
clear,clc,close all
local_dir = 'C:\Users\sminas\Desktop\PCE\';
cd([local_dir,'ANSYS\'])

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

% Number of experiments
K = 20; 
% Number of variables (Density, Ex = Ey = Ez, Prxy = Prxz = Prxyz)
M = 3; 
% Qtrans = load('ShearFrameSeismic4g\MProps.txt','-ascii');

load('ShearFrameSeismic\RandomInputs.mat');
Qtrans = [1:20 ;2e11 + 1e10*Q(1,:); 2e11 + 1e10*Q(1,:); ...
    0.065 + 0.025*Q(2,:); 0.065 + 0*Q(2,:); 350e6 + 150e6*Q(3,:) ; 350e6 + 150e6*Q(3,:)]';

figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
subplot(334),plot(Qtrans(:,1),Qtrans(:,2),'linewidth',1.5),
hold on,plot(Qtrans(:,1),Qtrans(:,3),'--g','linewidth',1.5),grid on
xlabel('Experiment number','Interpreter','Latex','Fontsize',9)
ylabel('$E$ (Pa)','Interpreter','Latex','Fontsize',9)
set(gca,'Fontsize',7,'xlim',[1,K],'ylim',[190e9 210e9],'Fontname','TimesNewRoman')
subplot(336),plot(Qtrans(:,1),Qtrans(:,4),'linewidth',1.5)
hold on,plot(Qtrans(:,1),Qtrans(:,5),'--g','linewidth',1.5),grid on
h = legend({'Vertical elements';'Horizontal elements'});
set(h,'Location','NorthOutside','Orientation','Horizontal','Position',[0.325 0.66 0.4 0.05],'Fontsize',9)
set(gca,'Fontsize',7,'xlim',[1,K],'ylim',[0.04 0.09],'Fontname','TimesNewRoman')
xlabel('Experiment number','Interpreter','Latex','Fontsize',9)
ylabel('$A$ (m$^2$)','Interpreter','Latex','Fontsize',9)
subplot(335),plot(Qtrans(:,1),Qtrans(:,6),'linewidth',1.5)
hold on,plot(Qtrans(:,1),Qtrans(:,7),'--g','linewidth',1.5),grid on
xlabel('Experiment number','Interpreter','Latex','Fontsize',9)
ylabel('$\sigma_Y$ (Pa)','Interpreter','Latex','Fontsize',9)
set(gca,'Fontsize',7,'xlim',[1,K],'ylim',[200e6 500e6],'Fontname','TimesNewRoman')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'mprops']),close;
     result = eps2xxx([write_dir,'mprops.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% PCE ARX Identification
clear,clc,close all
local_dir = 'C:\Users\sminas\Desktop\PCE\';
cd([local_dir,'ANSYS\'])

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

if ispc
    local_dir = 'C:\Users\sminas\Desktop\PCE\';
    cd('C:\Users\sminas\Desktop\PCE\m-files')
else
    local_dir = '/home/minas/Desktop/PCE/';
    cd('/home/minas/Desktop/PCE/m-files')
end

K = 20;
M = 3;
N = 1000;
n_max = 20;
q = 1;

load([local_dir,'ANSYS\ShearFrameSeismic\RandomInputs.mat']);
y = zeros(N,K);
Xstack = load([local_dir,'ANSYS\ShearFrameSeismic\elcentro.txt'],'-ascii');
for i = 1:K
    out = load([local_dir,'ANSYS\ShearFrameSeismic\SFres',num2str(i),'.txt']);
    y(:,i) = out(1:N,15);
    x(:,i) = Xstack(1:N);
end

options.maxsize = 1e9;
options.basis = 'legen';
options.focus = 'prediction';
% AR modelling
for na = 1:n_max
    for p = 1:4
        disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
        indx{1} = combinations(M,p,q);
        indx{2} = indx{1};
        B = length(indx{1});
        [aij,bij,res,criteria] = pcarx(y,x,Q,[na 0 na],indx,[],options);
        rss_sss(na,p) = criteria.rss_sss;
        connum(na,p) = criteria.connum;
        spp(na,p) = criteria.spp;
        sigma = var(res);
        bic(na,p) = log(sigma) + ((na+na+1)*B)*log(K*N)/(K*N);
    end
end

%% MSS criteria plot 
figure, 
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
clr = colormap(lines(4));
h(1) = axes('OuterPosition',[0 0.52 1 0.40]);
plot(2:n_max,rss_sss(2:n_max,1),'-o','Color',clr(1,:)),hold on
plot(2:n_max,rss_sss(2:n_max,2),'-d','Color',clr(2,:))
plot(2:n_max,rss_sss(2:n_max,3),'-s','Color',clr(3,:))
plot(2:n_max,rss_sss(2:n_max,4),'-*','Color',clr(4,:))
text(11,6.5,'PC-ARX($n,n$) models','Interpreter','Latex','Fontsize',12,'HorizontalAlignment','center')
text(0,3,'RSS/SSS (\%)','Interpreter','Latex','Fontsize',12,'HorizontalAlignment','center','Rotation',90)
set(gca,'xticklabel',[],'Fontsize',10,'Fontname','Timesnewroman','xlim',[1.5 20.5])
grid on
legend({'$P=1$','$P=2$','$P=3$','$P=4$'},'Interpreter','Latex','Fontsize',10)
annotation('rectangle',[0.375 0.565 0.2 0.08],'Color','k','Linestyle','--')
annotation('line',[0.375 0.5],[0.565 0.69],'Color','k','Linestyle','--')
annotation('line',[0.575 0.75],[0.565 0.69],'Color','k','Linestyle','--')
annotation('line',[0.375 0.5],[0.645 0.86],'Color','k','Linestyle','--')
annotation('line',[0.575 0.61],[0.645 0.69],'Color','k','Linestyle','--')
h(2) = axes('Position',[.5 .69 .25 .17],'Fontsize',8);
plot(8:12,rss_sss(8:12,1),'-o','Color',clr(1,:)),hold on
plot(8:12,rss_sss(8:12,2),'-d','Color',clr(2,:))
plot(8:12,rss_sss(8:12,3),'-s','Color',clr(3,:))
plot(8:12,rss_sss(8:12,4),'-*','Color',clr(4,:))
set(gca,'xtick',8:2:12,'Fontsize',10,'Fontname','Timesnewroman',...
    'xlim',[7.75 12.25],'ylim',[0 1.5])
grid on


h(3) = axes('OuterPosition',[0 0.15 1 0.40]);
plot(2:n_max,connum(2:n_max,1),'-o','Color',clr(1,:)),hold on
plot(2:n_max,connum(2:n_max,2),'-d','Color',clr(2,:))
plot(2:n_max,connum(2:n_max,3),'-s','Color',clr(3,:))
plot(2:n_max,connum(2:n_max,4),'-*','Color',clr(4,:))
text(0,1e10,'$\mathbf \Phi(\mathbf \Xi)$ condition number','Interpreter','Latex','Fontsize',12,'HorizontalAlignment','center','Rotation',90)
set(gca,'yscale','log','Fontsize',10,'Fontname','Timesnewroman','xlim',[1.5 20.5])
grid on
legend({'$P=1$','$P=2$','$P=3$','$P=4$'},'Interpreter','Latex','Fontsize',10)
text(11,6e-5,'AR/X models order ($n$)','Interpreter','Latex','Fontsize',12,'HorizontalAlignment','center')

if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'MSScriteria']),close;
     result = eps2xxx([write_dir,'MSScriteria.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% PC-ARX parameters
clear,clc,close all
local_dir = 'C:\Users\sminas\Desktop\PCE\';
cd([local_dir,'ANSYS\'])

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

if ispc
    local_dir = 'C:\Users\sminas\Desktop\PCE\';
    cd('C:\Users\sminas\Desktop\PCE\m-files')
else
    local_dir = '/home/minas/Desktop/PCE/';
    cd('/home/minas/Desktop/PCE/m-files')
end

K = 20;
M = 3;
N = 1000;
n_max = 20;
t0 = 50;
q = 1;

load([local_dir,'ANSYS\ShearFrameSeismic\RandomInputs.mat']);

y = zeros(N,K);
Xstack = load([local_dir,'ANSYS\ShearFrameSeismic\elcentro.txt'],'-ascii');
for i = 1:K
    out = load([local_dir,'ANSYS\ShearFrameSeismic\SFres',num2str(i),'.txt']);
    y(:,i) = out(1:N,15);
    x(:,i) = Xstack(1:N);
end

% ARX modelling
na = 10;
p = 2;
indx{1} = combinations(M,p,q);
indx{2} = indx{1};
options.focus = 'prediction';
options.maxsize = 1e9;
options.basis = 'legen';
[aij,bij] = pcarx(y,x,Q,[na 0 na],indx,[],options);
% PCE-AR parameters
[phi,an] = PCparam(Q,p,indx{1},aij,options);
% PCE-AR parameters
[phi,bn] = PCparam(Q,p,indx{2},bij,options);

figure(1)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
clr = colormap(lines(11));
symbols = ['-o';'-x';'-+';'-*';'-s';'-d';'-v';'-^';'-<';'->';'-p';'-h'];
subplot(7,1,1:3),hold on
for i = 1:10 
    plot(1:K,an(:,i),symbols(i,:),'Markersize',5,'linewidth',0.8,'Color',clr(i,:))
end
ylabel('AR parameters','Fontname','TimesNewRoman','Fontsize',12)
legend({'$a_1(\xi_k)$','$a_2(\xi_k)$','$a_3(\xi_k)$','$a_4(\xi_k)$','$a_5(\xi_k)$','$a_6(\xi_k)$',...
    '$a_7(\xi_k)$','$a_8(\xi_k)$','$a_9(\xi_k)$','$a_{10}(\xi_k)$'},'Interpreter','Latex',...
    'Fontsize',11,'Location','Eastoutside','Orientation','Vertical')
set(gca,'xticklabel',[],'Fontname','TimesNewRoman')
axis([0.5 20.5 -12 12])
grid on
box on
subplot(7,1,4:6),hold on
for i = 1:11
    plot(1:K,bn(:,i),symbols(i,:),'Markersize',4,'linewidth',0.8,'Color',clr(i,:))
end
ylabel('X parameters','Fontname','TimesNewRoman','Fontsize',12)
xlabel('Experiment number $k$','Interpreter','Latex','Fontname','TimesNewRoman','Fontsize',12)
legend({'$b_0(\xi_k)$','$b_1(\xi_k)$','$b_2(\xi_k)$','$b_3(\xi_k)$','$b_4(\xi_k)$','$b_5(\xi_k)$','$b_6(\xi_k)$',...
    '$b_7(\xi_k)$','$b_8(\xi_k)$','$b_9(\xi_k)$','$b_{10}(\xi_k)$'},'Interpreter','Latex',...
    'Fontsize',11,'Location','Eastoutside','Orientation','Vertical')
set(gca,'Fontname','TimesNewRoman')
axis([0.5 20.5 -0.4 0.8])
grid on
box on
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'ARXparams']),close;
     result = eps2xxx([write_dir,'ARXparams.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% Parameter surfaces
clear,clc,close all
local_dir = 'C:\Users\sminas\Desktop\PCE\';
cd([local_dir,'ANSYS\'])

write_dir = [local_dir,'Figures\'];
print_to_eps = 'y';
pictype = '-depsc2';
resolution = '-r300';

if ispc
    local_dir = 'C:\Users\sminas\Desktop\PCE\';
    cd('C:\Users\sminas\Desktop\PCE\m-files')
else
    local_dir = '/home/minas/Desktop/PCE/';
    cd('/home/minas/Desktop/PCE/m-files')
end

K = 20;
M = 3;
N = 1000;
t0 = 50;
q = 1;

load([local_dir,'ANSYS\ShearFrameSeismic\RandomInputs.mat']);
y = zeros(N,K);

Xstack = load([local_dir,'ANSYS\ShearFrameSeismic\elcentro.txt'],'-ascii');
for i = 1:K
    out = load([local_dir,'ANSYS\ShearFrameSeismic\SFres',num2str(i),'.txt']);
    y(:,i) = out(1:N,15);
    x(:,i) = Xstack(1:N);
end

% AR modelling
na = 10;
p = 2;
options.maxsize = 1e7;
options.basis = 'legen';
options.focus = 'prediction';

indx{1} = combinations(M,p,q);  indx{2} = indx{1};
[aij,bij] = pcarx(y,x,Q,[na 0 na],indx,[],options);
% AR modelling
disp(['AR order: ',num2str(na),' basis degree: ',num2str(p)]);
[phi,an] = PCparam(Q,p,indx{1},aij,options);

clear Qsurf
NoI = 100;
aux = linspace(-1,1,NoI);

for i=1:10
    figure(1)
    set(gcf,'PaperOrientation','portrait','papertype','A4',...
        'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
    Qsurf(1,:) = kron(ones(1,NoI),aux);
    Qsurf(2,:) = kron(aux,ones(1,NoI));
    Qsurf(3,:) = zeros(1,NoI^2);
    [phi,an_interp] = PCparam(Qsurf,p,indx{1},aij,options);
    subplot(2,3,2),surf(2e11 + 1e10*linspace(-1,1,NoI),0.065+0.025*linspace(-1,1,NoI),reshape(an_interp(:,i),NoI,NoI))
    shading interp
    set(gca,'Fontname','TimesNewRoman','Fontsize',6)
    xlabel('$E_x\ $ (Pa)','Interpreter','Latex','Fontsize',7)
    ylabel('$A\ $ (m$^2$)','Interpreter','Latex','Fontsize',7)
    title('$\sigma_Y = $ 350 (MPa)','Interpreter','Latex','Fontsize',7)
    zlabel(['$a_{',num2str(i),'}(\xi)$'],'Interpreter','Latex','Fontsize',7)
    axis tight
    grid on
    box on
    Qsurf(1,:) = kron(ones(1,NoI),aux);
    Qsurf(3,:) = kron(aux,ones(1,NoI));
    Qsurf(2,:) = zeros(1,NoI^2);
    [phi,an_interp] = PCparam(Qsurf,p,indx{1},aij,options);
    subplot(2,3,1),surf(2e11 + 1e10*linspace(-1,1,NoI),350e6+150e6*linspace(-1,1,NoI),reshape(an_interp(:,i),NoI,NoI))
    shading interp
    set(gca,'Fontname','TimesNewRoman','Fontsize',6)
    title('$A\ = $ 0.065 (m$^2$)','Interpreter','Latex','Fontsize',7)
    xlabel('$E_x\ $ (Pa)','Interpreter','Latex','Fontsize',7)
    ylabel('$\sigma_Y\ $ (Pa)','Interpreter','Latex','Fontsize',7)
    zlabel(['$a_{',num2str(i),'}(\xi)$'],'Interpreter','Latex','Fontsize',7)
    axis tight
    grid on
    box on
    Qsurf(2,:) = kron(ones(1,NoI),aux);
    Qsurf(3,:) = kron(aux,ones(1,NoI));
    Qsurf(1,:) = zeros(1,NoI^2);
    [phi,an_interp] = PCparam(Qsurf,p,indx{1},aij,options);
    subplot(2,3,3),surf(0.065 + 0.25*linspace(-1,1,NoI),350e6+150e6*linspace(-1,1,NoI),reshape(an_interp(:,i),NoI,NoI))
    shading interp
    set(gca,'Fontname','TimesNewRoman','Fontsize',6)
    title('$E_x\ = 200 $ (GPa)','Interpreter','Latex','Fontsize',7)
    xlabel('$A\ $ (m$^2$)','Interpreter','Latex','Fontsize',7)
    ylabel('$\sigma_Y\ $ (Pa)','Interpreter','Latex','Fontsize',7)
    zlabel(['$a_{',num2str(i),'}(\xi)$'],'Interpreter','Latex','Fontsize',7)
    grid on
    box on
    axis tight
    if print_to_eps=='y';
        print(pictype,resolution,'-zbuffer',[write_dir,'ARsurfs',num2str(i)]),close;
        result = eps2xxx([write_dir,'ARsurfs',num2str(i),'.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    end
end

for i=1:11
    figure(1)
    set(gcf,'PaperOrientation','portrait','papertype','A4',...
        'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
    Qsurf(1,:) = kron(ones(1,NoI),aux);
    Qsurf(2,:) = kron(aux,ones(1,NoI));
    Qsurf(3,:) = zeros(1,NoI^2);
    [phi,bn_interp] = PCparam(Qsurf,p,indx{2},bij,options);
    subplot(2,3,2),surf(2e11 + 1e10*linspace(-1,1,NoI),0.065+0.025*linspace(-1,1,NoI),reshape(bn_interp(:,i),NoI,NoI))
    shading interp
    set(gca,'Fontname','TimesNewRoman','Fontsize',6)
    xlabel('$E_x\ $ (Pa)','Interpreter','Latex','Fontsize',7)
    ylabel('$A\ $ (m$^2$)','Interpreter','Latex','Fontsize',7)
    title('$\sigma_Y = $ 350 (MPa)','Interpreter','Latex','Fontsize',7)
    zlabel(['$b_{',num2str(i-1),'}(\xi)$'],'Interpreter','Latex','Fontsize',7)
    axis tight
    grid on
    box on
    Qsurf(1,:) = kron(ones(1,NoI),aux);
    Qsurf(3,:) = kron(aux,ones(1,NoI));
    Qsurf(2,:) = zeros(1,NoI^2);
    [phi,bn_interp] = PCparam(Qsurf,p,indx{2},bij,options);
    subplot(2,3,1),surf(2e11 + 1e10*linspace(-1,1,NoI),350e6+150e6*linspace(-1,1,NoI),reshape(bn_interp(:,i),NoI,NoI))
    shading interp
    set(gca,'Fontname','TimesNewRoman','Fontsize',6)
    title('$A\ = $ 0.065 (m$^2$)','Interpreter','Latex','Fontsize',7)
    xlabel('$E_x\ $ (Pa)','Interpreter','Latex','Fontsize',7)
    ylabel('$\sigma_Y\ $ (Pa)','Interpreter','Latex','Fontsize',7)
    zlabel(['$b_{',num2str(i-1),'}(\xi)$'],'Interpreter','Latex','Fontsize',7)
    axis tight
    grid on
    box on
    Qsurf(2,:) = kron(ones(1,NoI),aux);
    Qsurf(3,:) = kron(aux,ones(1,NoI));
    Qsurf(1,:) = zeros(1,NoI^2);
    [phi,bn_interp] = PCparam(Qsurf,p,indx{2},bij,options);
    subplot(2,3,3),surf(0.065 + 0.25*linspace(-1,1,NoI),350e6+150e6*linspace(-1,1,NoI),reshape(bn_interp(:,i),NoI,NoI))
    shading interp
    set(gca,'Fontname','TimesNewRoman','Fontsize',6)
    title('$E_x\ = 200 $ (GPa)','Interpreter','Latex','Fontsize',7)
    xlabel('$A\ $ (m$^2$)','Interpreter','Latex','Fontsize',7)
    ylabel('$\sigma_Y\ $ (Pa)','Interpreter','Latex','Fontsize',7)
    zlabel(['$b_{',num2str(i-1),'}(\xi)$'],'Interpreter','Latex','Fontsize',7)
    grid on
    box on
    axis tight
    if print_to_eps=='y';
        print(pictype,resolution,'-zbuffer',[write_dir,'Xsurfs',num2str(i)]),close;
        result = eps2xxx([write_dir,'Xsurfs',num2str(i),'.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
    end
end



%% Predictions - Simulations
clear,clc,close all
local_dir = 'C:\Users\sminas\Desktop\PCE\';
cd([local_dir,'ANSYS\'])

write_dir = [local_dir,'Figures\'];
print_to_eps = 'n';
pictype = '-depsc2';
resolution = '-r300';

if ispc
    local_dir = 'C:\Users\sminas\Desktop\PCE\';
    cd('C:\Users\sminas\Desktop\PCE\m-files')
else
    local_dir = '/home/minas/Desktop/PCE/';
    cd('/home/minas/Desktop/PCE/m-files')
end

K = 20;
M = 3;
N = 1000;
q = 1;
Ts = 1/64;

load([local_dir,'ANSYS\ShearFrameSeismic\RandomInputs.mat']);
y = zeros(N,K);

Xstack = load([local_dir,'ANSYS\ShearFrameSeismic\elcentro.txt'],'-ascii');
for i = 1:K
    out = load([local_dir,'ANSYS\ShearFrameSeismic\SFres',num2str(i),'.txt']);
    y(:,i) = out(1:N,15);
    x(:,i) = Xstack(1:N);
end

% ARX modelling
na = 10;
p = 2;
indx{1} = combinations(M,p,q);  indx{2} = indx{1};
options.focus = 'prediction';
options.maxsize = 1e9;
options.basis = 'legen';
[aij0,bij0] = pcarx(y,x,Q,[na 0 na],indx,[],options);
options.focus = 'simulation';
[aij,bij] = pcarx(y,x,Q,[na 0 na],indx,[aij0;bij0],options);

T0 = 160;
Tf = 350;

Qtrans = [1.95e11; 0.045; 250e6];
Q = [(Qtrans(1)-2e11)/1e10 ,(Qtrans(2)-0.065)/0.025 ,(Qtrans(3)-350e6)/150e6]';
pacoima = load([local_dir,'ANSYS\ShearFrameSeismic\pacoima.txt']);
yr = load([local_dir,'ANSYS\ShearFrameSeismic\SFresPACD.txt']);

% PCE-AR parameters
[~,an0] = PCparam(Q,p,indx{1},aij0,options);
% PCE-AR parameters
[~,bn0] = PCparam(Q,p,indx{2},bij0,options);
% PCE-AR parameters
[~,an] = PCparam(Q,p,indx{1},aij,options);
% PCE-AR parameters
[~,bn] = PCparam(Q,p,indx{2},bij,options);

% Dynamic response one-step-ahead predictions
ResDat = iddata(yr(:,15),pacoima,Ts);
model = idpoly([1 an],bn);
yp = predict(model,ResDat,1);

% Dynamic response simulations
ys = filter(bn,[1 an],pacoima,zeros(na,1));


figure(2)
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
clr = colormap(lines(3));
subplot(3,1,1),plot(Ts:Ts:Ts*N,pacoima,'linewidth',1)
set(gca,'xlim',[Ts Ts*1000],'ylim',[-20 20],'FontName','TimesNewRoman','Fontsize',11)
ylabel('$\alpha$ (m/s$^2$)','Interpreter','Latex','Fontsize',14)
xlabel('Time (s)','FontName','TimesNewRoman')
grid on
subplot(3,1,2:3)
plot(T0*Ts:Ts:Tf*Ts,yr(T0:Tf,15),'-bo','linewidth',0.8,'Markersize',5),hold on
plot(T0*Ts:Ts:Tf*Ts,yp.outputdata(T0:Tf),'rx','linewidth',1,'Markersize',5),hold on
plot(T0*Ts:Ts:Tf*Ts,ys(T0:Tf),'g+','linewidth',1,'Markersize',5),
legend({['$y[t]$'];'$\hat{y}[t|t-1]$';'$\bar{y}[t]$'},'Interpreter','Latex','Fontsize',12,'Orientation','Horizontal')
set(gca,'xlim',[T0/64 Tf/64],'ylim',[-30 20],'FontName','TimesNewRoman','Fontsize',11)
ylabel('Acceleration (m/s$^2$)','Interpreter','Latex','Fontsize',14)
xlabel('Time (s)','FontName','TimesNewRoman','Fontsize',12)
grid on
box on
if print_to_eps == 'y'
    print(pictype,resolution,[write_dir,'predsPACDsim']),close;
    result = eps2xxx([write_dir,'predsPACDsim.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end

pred_err = 100*norm(yr(:,15)-yp.outputdata)^2/norm(yr(:,15))^2
sim_err = 100*norm(yr(:,15)-ys)^2/norm(yr(:,15))^2