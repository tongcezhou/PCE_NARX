%%  *****  INDEX OF DATA SETS ON RESULTS FILE  *****
% =======================================================
%    SET   TIME/FREQ    LOAD STEP   SUBSTEP  CUMULATIVE
%      3 0.73404             1         3         3
%      6  2.1501             1         6         6
%      8  3.4322             1         8         8
%     10  4.4252             1        10        10
%     12  5.0605             1        12        12
    
%% Section integration points
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd(local_dir)
ANSYSplots('preliminary/SF2Dtests3.out');




%% Sinusodial acceleration (ANSYS run)
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd(local_dir)

Ts = 0.05;
time = Ts*[1:1000];
Acc = 25*sin(time*pi/12.5)';
Acc = (linspace(0,1,1000)').*Acc;
plot(time,Acc)
save('preliminary\testINincr.txt','Acc','-ascii')


%% Comparison wrt to time substeps
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'preliminary/'])

stress10 = load(['SF2Dtest3_10_Sx_1.txt']);
stress100 = load(['SF2Dtest3_Sx_1.txt']);
stress1000 = load(['SF2Dtest3_1000_Sx_1.txt']);

figure, 
clr = colormap(lines(3)); 
hold on
plot(stress10(1:1000,3),'-','Color',clr(1,:))
plot(stress100(1:1000,3),'--','Color',clr(2,:))
plot(stress1000(1:1000,3),'-.','Color',clr(3,:))
legend({'10 substeps','100 substeps','1000 substeps'})


%% Comparison wrt to number of column elements
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'preliminary/'])

stress1 = load(['SF2Dtest1_Sx_1.txt']);
stress3 = load(['SF2Dtest3_Sx_1.txt']);
stress5 = load(['SF2Dtest5_Sx_1.txt']);
stress10 = load(['SF2Dtest10_Sx_1.txt']);

figure, 
clr = colormap(lines(4)); 
hold on
plot(stress1(:,3),'-','Color',clr(1,:))
plot(stress3(:,3),'--','Color',clr(2,:))
plot(stress5(:,3),'-.','Color',clr(3,:))
plot(stress10(:,3),':','Color',clr(4,:))
title('BEAM188')
legend({'1 element','3 elements','5 elements','10 elements'})


%% BEAM188 - increasing load
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'preliminary/'])

stress = load(['SF2Dtest1_incr_Sx_1.txt']);
strain = load(['SF2Dtest1_incr_Ex_1.txt']);

figure(1)
plot(stress(:,3:10))

figure(2)
plot(strain(:,3:10),stress(:,3:10))



%% Comparison between BEAM23 and BEAM188
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'preliminary/'])

stress188 = load(['SF2Dtest1_Sx_1.txt']);
strain188 = load(['SF2Dtest1_Ex_1.txt']);
stress23 = load(['SF2Dtest1legacy_Sx_1.txt']);
% strainP23 = load(['SF2Dtest1legacy_Epx_1.txt']);
% strainE23 = load(['SF2Dtest1legacy_Eex_1.txt']);
% dspl23 = load(['SF2Dtest1legacy_dspl.txt']);


figure(1) 
clr = colormap(lines(3)); 
hold on
plot(stress188(1:end,3:end),'Color',clr(1,:))
plot(stress23(1:end,3:end),'Color',clr(2,:))


%%
figure(2)
clr = colormap(lines(2)); 
plot(strain188(:,3),stress188(:,3),'Color',clr(1,:))
hold on
plot(strainE23(:,3)+strainP23(:,3),stress23(:,3),'Color',clr(2,:))



%% Comparison between BEAM23 and BEAM188
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'preliminary/'])

stress188 = load(['SF2Dtest1_Sx_1.txt']);
strain188 = load(['SF2Dtest1_Ex_1.txt']);
stress23 = load(['SF2Dtest1legacy_Sx_1.txt']);
% strainP23 = load(['SF2Dtest1legacy_Epx_1.txt']);
% strainE23 = load(['SF2Dtest1legacy_Eex_1.txt']);
% dspl23 = load(['SF2Dtest1legacy_dspl.txt']);
% dspl23new = load(['SF2Dtest1legacy_dspl_new.txt']);

figure(1) 
clr = colormap(lines(2)); 
hold on
plot(stress188(1:end,3:end),'Color',clr(1,:))
plot(stress23(1:end,3:end),'Color',clr(2,:))


figure(2)
clr = colormap(lines(2)); 
plot(strain188(:,3),stress188(:,3),'Color',clr(1,:))
hold on
plot(strainE23(:,3)+strainP23(:,3),stress23(:,3),'Color',clr(2,:))


%% BEAM23 results
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'preliminary/'])

for i = 1:15
    stress = load(['SF2Dtest1legacy_Sx_',num2str(i),'.txt']);
    % strainP = load(['SF2Dtest1legacy_Epx_1.txt']);
    % strainE = load(['SF2Dtest1legacy_Eex_1.txt']);
    figure
    plot(stress(:,3:end))
end


%%
IP = 3;
figure(2)
clr = colormap(lines(2)); 
hold on
plot(strainE(:,3:end)+strainP(:,3:end),stress(:,3:end),'Color',clr(2,:))



%% Section integration points
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'preliminary/'])

figure, 
clr = colormap(lines(35)); 
hold on
for elem = 1:30
    stress = load(['SF2Dtest3_Sx_',num2str(elem),'.txt']);
    strain = load(['SF2Dtest3_Ex_',num2str(elem),'.txt']);
    % plot(stress(:,3:end))
    % plot(strain(:,3:end),stress(:,3:end),'-','Color',clr(elem,:))
    % plot(stress(:,3),'-o','Color',clr(elem,:))    
    % plot(stress(:,3),'-o','Color',clr(elem,:))
    % figure(elem)
    % plot(stress(:,3:10),'-b'),hold on,plot(stress(:,11:18),'-r')
    plot(stress(:,3:end),'Color',clr(elem,:))
end

dspl = load('SF2Dtest3_dspl.txt');
vel = load('SF2Dtest3_vel.txt');
acc = load('SF2Dtest3_acc.txt');


%%
Sx = load(['SF2Dtest3_Sx_1.txt']);
Sy = load(['SF2Dtest3_Sy_1.txt']);
Sz = load(['SF2Dtest3_Sz_1.txt']);
Sxy = load(['SF2Dtest3_Sxy_1.txt']);
Syz = load(['SF2Dtest3_Syz_1.txt']);
Sxz = load(['SF2Dtest3_Sxz_1.txt']);

IP = 5;
figure, 
clr = colormap(lines(6)); 
hold on
plot(Sx(:,IP),'-o','Color',clr(1,:))
plot(Sy(:,IP),'Color',clr(2,:))
plot(Sz(:,IP),'Color',clr(3,:))
plot(Sxy(:,IP),'Color',clr(4,:))
plot(Syz(:,IP),'Color',clr(5,:))
plot(Sxz(:,IP),'Color',clr(6,:))


%%
N = 1000;

figure(1),
subplot(311),plot(0.025*(1:N),SF2D_test_strprin1(:,2))
xlabel('Time (s)')
ylabel('Input acc (m/s^2)')
subplot(3,1,2:3),plot(0.025*(1:N),SF2D_test_strprin1(:,4:2:13)/1e6,'-')
hold on,plot(0.025*(1:N),SF2D_test_strprin1(:,5:2:13)/1e6,'-.')
xlabel('Time (s)')
ylabel('Stress (MPa)')
legend({'Elem. 1','Elem. 3','Elem. 5','Elem. 7','Elem. 9','Elem. 2','Elem. 4','Elem. 6','Elem. 8','Elem. 10',})
print('-dpng','-r300','stresses_sin_zoom.png')

figure(2),
plot(SF2D_test_strain1(:,13),SF2D_test_strprin1(:,13),'-o')
xlabel('')
ylabel('')


%% Space frame -- bilinear isotropic material (increasing force) 
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\EqEngAppSHM\PreliminaryAnalysis\')

NoS = 2500;
NoExp = 1;
% Generation of the input force
Finput = linspace(5e2,1e5,NoS)';                                        % Linear
save('InputForce.txt','Finput','-ascii')

% Run Ansys
delete('file.*')
% dos(' "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ANSYS140.exe" -b -i "C:\Users\sminas\Dropbox\MATLAB\EqEngAppSHM\PreliminaryAnalysis\SpaceFrame2DForce.txt" -o "output.txt"');


%% Space frame -- bilinear isotropic material (increasing force) -- Figure 
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\EqEngAppSHM\PreliminaryAnalysis\')

load('SF_ForceResponse.txt');
time = SF_ForceResponse(:,1);
Fi = SF_ForceResponse(:,2)*36;
sumreact = SF_ForceResponse(:,3);
Sx = SF_ForceResponse(:,1:5);
Sprin = SF_ForceResponse(:,6:10);

figure,
subplot(4,1,1),plot(time,sumreact/1e3)
xlabel('Time (s)')
ylabel('F_{X_{reac}} (kN)')
subplot(4,2,3:2:7),plot(time,Sprin/1e6,'-o','Markersize',2)
legend({'Element 2','Element 4','Element 6','Element 8','Element 10'},'Location','NorthWest')
xlabel('Time (s)')
ylabel('\sigma_1 (MPa)')
axis tight
subplot(4,2,4:2:8),plot(time,Sx/1e6,'-') 
legend({'Element 1','Element 2','Element 3','Element 4','Element 5',...
    'Element 6','Element 7','Element 8','Element 9','Element 10'},'Location','NorthWest')
xlabel('Time (s)')
ylabel('\sigma_X (MPa)')
axis tight
print('-dpng','-r300','stresses.png')


%% Steel Properties
% Properties                    Carbon     Alloy        Stainless     Tool
% ==========================================================================
% Density (1000 kg/m3)        7.85          7.85        7.75-8.1    7.72-8.0
% Elastic Modulus (GPa)       190-210       190-210     190-210     190-210
% Poisson's Ratio             0.27-0.3      0.27-0.3    0.27-0.3    0.27-0.3
% Tensile Strength (MPa)      276-1882      758-1882    515-827     640-2000
% Yield Strength (MPa)        186-758       366-1793    207-552     380-440
% ==========================================================================