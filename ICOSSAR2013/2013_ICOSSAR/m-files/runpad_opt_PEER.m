%% Save earthquake accelerograms in one file
clear,clc,pack
local_dir = [pwd,'\'];

cd([local_dir,'PEER'])
listing = ls('*.AT2');

p = 0;
for j = 1:size(listing,1)
    if isempty(strfind(listing(j,:),'DWN')) && isempty(strfind(listing(j,:),'UP'))
        A = importdata(listing(j,:),' ',4);
        A.data = A.data';
        % if length(A.data(:)) < 10000
            p = p+1;
            N(p) = length(A.data(:));
            EQ{p} = A.data(:);
            if isfield(A,'colheaders')
                disp('Y')
                DT(p) = str2double(A.colheaders{4});
            else
                k = strfind(A.textdata{4},' ');
                NoS(p) = str2double(A.textdata{4}(1:k(1)-1));
                q = 1;
                while (k(2)-k(1))==1,k = k(2:end);end
                DT(p) = str2double(A.textdata{4}(k(1)+1:k(2)));
            end
        % end
    end
end
clear A j
save('PEER_ALL','EQ','N','DT');


%% PEER Database criteria computation
clear,clc,pack
local_dir='C:\Users\echatzi\Dropbox\Conferences\ICOSSAR2013\2013_ICOSSAR\m-files\'
load([local_dir,'PEER\PEER_ALL'],'EQ','N','DT');

for j = 1:length(EQ)
    disp(j)
    Y = EQ{j};
    Y = Y(~isnan(Y));
    N = length(Y);
    Ts = DT(j);
    % a) Peak ground acceleration
    PGA(j) = max(abs(Y));
    % b) Root-mean-square acceleration
    RMSA(j) = sqrt((1/N)*sum(Y.^2));
end
figure(j)
subplot(211),plot(sort(PGA))
subplot(212),plot(sort(RMSA))


%% Modulating function parameters optimization (Arias Intensity)
clear,clc,pack
local_dir = 'C:\Users\echatzi\Dropbox\Conferences\ICOSSAR2013\2013_ICOSSAR\m-files\';
load([local_dir,'PEER\PEER_ALL'],'EQ','N','DT');

PSOoptions.StallGenLimit = 20; 
PSOoptions.PopulationSize = 100;
PSOoptions.InitialPopulationSize = 100;
PSOoptions.Generations = 500;
PSOoptions.TolFun = 1e-6;
PSOoptions.SocialAttraction = 1.45;
PSOoptions.CognitiveAttraction = 1.45;
PSOoptions.Display = 'off';
PSOoptions.Vectorized = 'on';
PSOoptions.PlotFcns = {};%{@psoplotbestf,@psoplotswarm,@psoplotscorediversity}; 
PSOoptions.ConstrBoundary = 'soft';

NLoptions = optimset('Algorithm','interior-point',...
    'MaxFunEvals',20000,'MaxIter',20000,'Display','off',...
    'TolFun',1e-4,'TolX',1e-8,'TolCon',1e-6);

nvar = 2;   
theta = zeros(length(EQ),4);
th0 = [1.01 0.5];           % a2 and a3 initial values
lb = [1+eps eps];           % a2 and a3 lower bound
ub = [inf 1];               % a2 and a3 upper bound
[RSS, extflg, Ia, D5_45, D5_95, IAs, D545s, D595s] = deal(zeros(length(EQ),1)); 
plotfigs = 'n';

for j = 1:2  %length(EQ)
    disp(j)
    Y = EQ{j};
    Y = Y(~isnan(Y));
    N = length(Y);
    Ts = DT(j);
    % PSO optimization
    thPSO = pso(@(theta0) Qoptimization(theta0,Y,Ts),nvar,[],[],[],[],lb,ub,[],PSOoptions);
    % Nonlinear optimization
    [thNL,~,extflg(j)] = fmincon(@(theta0) Qoptimization(theta0,Y,Ts),thPSO,[],[],[],[],lb,ub,[],NLoptions);
    % Arias intensity calculation
    [Ia(j),t5,t45,t95,Iat,D5_45(j),D5_95(j)] = EQarias(Y,Ts);
    % Normalized Arias intensity as a function of time
    Iaindx = Iat./Ia(j);
    [RSS(j),theta(j,:),IAs(j),t5s,t45s,t95s,IATs,D545s(j),D595s(j)] = Qoptimization(thNL,Y,Ts);
    T = 0:Ts:Ts*(N-1);
    T0 = theta(j,4);
    if strcmp(plotfigs,'y');
        figure(j)
        plot(T/(Ts*(N-1)),Iaindx),hold on
        plot([t5 t5 t5],[0 0.5 1],'-or',[t45 t45 t45],[0 0.5 1],'-^r',[t95 t95 t95],[0 0.5 1],'-sr')
        plot(T0+T/(Ts*(N-1)),IATs,'-c')
        plot(T0+[t5s t5s t5s],[0 0.5 1],'-om',T0+[t45s t45s t45s],[0 0.5 1],'-^m')
        plot(T0+[t95s t95s t95s],[0 0.5 1],'-sm')
        hold off
        xlabel('Sample')
        ylabel('Normalized Arias intensity')
    end
end
%save('PEER_ALL_MODopt_Arias.mat','theta*','RSS*','ext*','IAs','Ia','D5*')

%% Frequency-domain parameters optimization (omega - interpolation)
clear,clc,pack
load('D:\MATLAB\PEER\PEER_ALL.mat','EQ','N','DT');
cd('C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\');

options.points = 9;
for j = 1:length(EQ)
    disp(j)
    Y = EQ{j};
    Y = Y(~isnan(Y));
    N(j) = length(Y);
    Ts = DT(j);
    [theta(j,:),tmid(j)] = spectralINTERP(Y,Ts,options);
    theta(j,:) = theta(j,:);
end
save('PEER_ALL_IRFinterp.mat','theta*','tmid')


%% Frequency-domain parameters optimization (zeta - interpolation)
clear,clc,pack
load('D:\PEER\PEER_ALL','EQ','N','DT');
cd('C:\Users\sminas\Dropbox\MATLAB\ICOSSAR2013\');

IRFth = load('PEER_ALL_IRFinterp.mat','theta*');
load('PEER_ALL_IRFinterp.mat','tmid');
qth = load('PEER_ALL_MODopt.mat','theta*');
options.iter = 10;

Zinit_indx = 0.1:0.1:0.9;
Zindx = zeros(length(EQ),19);
RSSref = zeros(length(EQ),19);

for j = 1:length(EQ)
    disp(j)
    Y = EQ{j};
    Y = Y(~isnan(Y));
    Ts = DT(j);
    p=0;
    if N(j)<10000
        RSS = zeros(length(Zinit_indx),1);
        for i = Zinit_indx
            p = p+1;
            zeta = i;
            TH = [IRFth.theta(j,:) zeta, qth.theta(j,1:3)];
            [RSS(p),PN,PNs] = bandwidthERRinterp(TH,Y,Ts,tmid(j),options);
            if p>2 && (RSS(p)>RSS(p-1)) && (RSS(p-1)>RSS(p-2))
                break
            end
        end
        RSS(RSS == 0) = NaN;
        [minRSS,J] = min(RSS);
        Zindx(j,:) = Zinit_indx(J)-0.09:0.01:Zinit_indx(J)+0.09;
        q = 0;
        for i = Zindx(j,:)
            q = q+1;
            zeta = i;
            TH = [IRFth.theta(j,:) zeta, qth.theta(j,1:3)];
            [RSSref(j,q),PN,PNs] = bandwidthERRinterp(TH,Y,Ts,tmid(j),options);
        end
        [minRSSref,K] = min(RSSref(j,:));
        zeta_ref(j) = Zindx(j,K);
    end
    save('PEER_Nleq10000_ZETAinterp.mat','RSSref','zeta_ref')
end