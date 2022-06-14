%%  *****  INDEX OF DATA SETS ON RESULTS FILE  *****
% =======================================================
%  SET   TIME/FREQ
%   1   2.3689
%   2   5.0725
%   3   8.1379
%   4   12.304
%   5   18.816


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
time = Ts*[1:2000];
Acc = 50*sin(time*pi/12)';
Acc = (linspace(0.2,1,2000)').*Acc;
plot(time,Acc)
save('final\testIN.txt','Acc','-ascii')



%% Space frame -- bilinear isotropic material (increasing force) 
% =========================================================================
clear,clc,close all
cd('C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\')

NoS = 2000;
NoExp = 1;
% Generation of the input force
Finput = linspace(1e2,1e4,NoS)';                                        % Linear
save('final\InputForce.txt','Finput','-ascii')

% Run Ansys
delete('file.*')


%% BEAM23 results (increasing sinusoidal acceleration)
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'final/'])

for i = 1:15
    stress = load(['SF2Dtest1legacy_Sx_',num2str(i),'.txt']);
    time = stress(:,1);
    figure
    plot(time,stress(:,3:4),'-b'),hold on
    plot(time,stress(:,5:6),'--r')
    plot(time,stress(:,7:8),'-.k')
    title(['Element ',num2str(i)])
    ylabel('Stress (Pa)')
    xlabel('Time (s)')
    legend({'Max (bottom)','Min (bottom)','Max (middle)','Min (middle)','Max (top)','Min (top)'})
end




%% BEAM23 results (increasing sinusoidal acceleration)
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'final/'])

for i = 1:15
    stress = load(['SF2Dtest1etable_Sx_',num2str(i),'.txt']);
    time = stress(:,1);
    figure
    plot(time,stress(:,3:4),'-b'),hold on
    plot(time,stress(:,5:6),'--r')
    plot(time,stress(:,7:8),'-.k')
    title(['Element ',num2str(i)])
    ylabel('Stress (Pa)')
    xlabel('Time (s)')
    legend({'Max (bottom)','Min (bottom)','Max (middle)','Min (middle)','Max (top)','Min (top)'})
end




%% BEAM23 results (increasing force)
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'final/'])

for i = 1:15
    stress = load(['SF2Dtest1force_Sx_',num2str(i),'.txt']);
    time = stress(:,1);
    figure
    plot(time,stress(:,3:4),'-b'),hold on
    plot(time,stress(:,5:6),'--r')
    plot(time,stress(:,7:8),'-.k')
    title(['Element ',num2str(i)])
    ylabel('Stress (Pa)')
    xlabel('Time (s)')
    legend({'Max (bottom)','Min (bottom)','Max (middle)','Min (middle)','Max (top)','Min (top)'})
end



%% BEAM188 - increasing load
%==========================================================================
clear,clc,close all
if ispc
    local_dir = 'C:\Users\sminas\Dropbox\MATLAB\SF2Dtests\';
else
    local_dir = '/home/minas/Dropbox/MATLAB/SF2Dtests/';
end
cd([local_dir,'final/'])


for k=1:15
    stress = load(['SF2Dtest1_Sx_',num2str(k),'.txt']);
    strain = load(['SF2Dtest1_Ex_',num2str(k),'.txt']);
    figure(k)
    subplot(211),plot(stress(:,3:10))
    subplot(212),plot(strain(:,3:10),stress(:,3:10))
end



%% Steel Properties
% Properties                    Carbon     Alloy        Stainless     Tool
% ==========================================================================
% Density (1000 kg/m3)        7.85          7.85        7.75-8.1    7.72-8.0
% Elastic Modulus (GPa)       190-210       190-210     190-210     190-210
% Poisson's Ratio             0.27-0.3      0.27-0.3    0.27-0.3    0.27-0.3
% Tensile Strength (MPa)      276-1882      758-1882    515-827     640-2000
% Yield Strength (MPa)        186-758       366-1793    207-552     380-440
% ==========================================================================