%% Hyperbolic PCE (Fig. 2 [1])
clear,clc,close all
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';
M = 2;  	  % number of variables
p = 3;	      % total maximum degree
[indx,P,J] = combinations(M,p,1);


figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
% NL regressors
subplot(4,3,[1 2 4 5 7 8])
degree = 3;
maxterms = 1;
na = 4;		nd=0;	nb=4;
plot(0,0)
str = {};
for pow = 1:degree
    for j=1:na+nb-nd+1
        if j<=na
            if pow >1
                str = [str, ['$y[t-',num2str(j),']^',num2str(pow),'$']]
            else
                str = [str, ['$y[t-',num2str(j),']$']]
            end
        elseif j==na+1
            if pow >1
                str = [str, ['$x[t]^',num2str(pow),'$']]
            else
                str = [str, ['$x[t]$']]
            end
        else
            if pow >1
                str = [str, ['$x[t-',num2str(j-na-1),']^',num2str(pow),'$']]
            else
                str = [str, ['$x[t-',num2str(j-na-1),']$']]
            end
        end
    end
end
axis([0 degree*(na+nb-nd+1) 0 degree*(na+nb-nd+1)]);
set(gca,'xtick',[0:length(str)],'xticklabel',[],'ytick',[0:length(str)],'yticklabel',[])
text(-1,-1,'1','Interpreter','Latex','HorizontalAlignment','Center')
for k = 1:length(str)
    text(k,-1.5,str{k},'Interpreter','Latex','Rotation',90,'HorizontalAlignment','Center','Fontsize',7)
end
for k = 1:(degree-1)*(na+nb-nd+1)
    text(k,-1.5,str{k},'Interpreter','Latex','Rotation',90,'HorizontalAlignment','Center','Fontsize',7)
    text(-1.5,k,str{k},'Interpreter','Latex','HorizontalAlignment','Center','Fontsize',7)
end

hold on
grid on
aux = (1:length(str))';
aux = kron(ones(1,length(str)),aux);
plot(aux,aux','rx')
plot(0,0,'sb','MarkerFacecolor','b')
% plot(0,1:length(str),'sb','MarkerFacecolor','b')
plot(1:length(str),zeros(1,length(str)),'sb','MarkerFacecolor','b')
% plot(zeros(1,length(str)),1:length(str),'sb','MarkerFacecolor','b')


for ii = 1:length(str)
    powIindx = findstr(str{ii},'^');
    if isempty(powIindx), powI = 1;     else powI = str2num(str{ii}(powIindx+1)); end
    for jj = 1:length(str)
        powJindx = findstr(str{jj},'^');
        if isempty(powJindx), powJ = 1; else powJ = str2num(str{jj}(powJindx+1)); end
        if (powI + powJ) > 1 && (powI + powJ) <= degree
            plot(ii,jj,'sb','MarkerFacecolor','b')
            plot(jj,ii,'sb','MarkerFacecolor','b')
        end
    end
end
title('Nonlinear regressors','Fontname','TimesNewRoman')
ylim([0 (degree-1)*(na+nb-nd+1)])
% PC basis
subplot(4,3,6),hold on, box on, grid on
for i=0:3
    for j = 0:3
        plot(i,j,'rx','Markersize',5)
    end
end
for i = 1:size(indx,1)
    plot(indx(i,1),indx(i,2),'s','MarkerFacecolor','b')
end
axis([0 3 0 3])
set(gca,'xtick',0:3,'ytick',0:3)
xlabel('$\xi_1$','Interpreter','Latex','HorizontalAlignment','Center')
ylabel('$\xi_2$','Interpreter','Latex','HorizontalAlignment','Center')
title('PC basis','Fontname','TimesNewRoman')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'searchspace']),close;
     result = eps2xxx([write_dir,'searchspace.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end


%% Hyperbolic PCE (Fig. 2 [1])
clear,clc,close all
write_dir = 'C:\Users\sminas\Dropbox\MATLAB\PCNARXjournal\Figures\';
pictype = '-depsc2';
resolution = '-r300';
print_to_eps = 'y';
M = 2;  	  % number of variables
p = 3;	      % total maximum degree
[indx,P,J] = combinations(M,p,1);


figure(1),
set(gcf,'PaperOrientation','portrait','papertype','A4',...
    'paperunits','centimeters','paperPosition',[0.63, 7.75, 19.72, 14.17])
% NL regressors
subplot(2,3,[1 2])
degree = 3;
maxterms = 1;
na = 4;		nd=0;	nb=4;
plot(0,0)
str = {'0'};
for j=1:na+nb-nd+1
    if j<=na
        str = [str, ['$y[t-',num2str(j),']$']]
    elseif j==na+1
        str = [str, ['$x[t]$']]
    else
        str = [str, ['$x[t-',num2str(j-na-1),']$']]
    end
end
axis([0 (na+nb-nd+1) 0 degree]);
set(gca,'xtick',[0:length(str)],'xticklabel',[],'ytick',[0:degree],'yticklabel',[0:degree])
for k = 0:(na+nb-nd+1)
    text(k,-0.2,str{k+1},'Interpreter','Latex','Rotation',0,'HorizontalAlignment','Center','Fontsize',7)
end

hold on
grid on
aux = (1:length(str))';
aux = kron(ones(1,degree+1),aux);
plot(aux,0:degree,'rx')
plot(0,0,'sb','MarkerFacecolor','b')
plot(1:length(str),zeros(1,length(str)),'sb','MarkerFacecolor','b')
plot(aux,0:degree,'sb','MarkerFacecolor','b')
text(4.5,-0.75,'Regressors','Fontname','TimesNewRoman','HorizontalAlignment','Center','Fontsize',9)
ylabel('Exponential','Fontname','TimesNewRoman','Fontsize',9)
title('Nonlinear regressors','Fontname','TimesNewRoman','Fontsize',9)
ylim([0 degree])


% PC basis
subplot(2,3,3),hold on, box on, grid on
for i=0:3
    for j = 0:3
        plot(i,j,'rx','Markersize',5)
    end
end
for i = 1:size(indx,1)
    plot(indx(i,1),indx(i,2),'s','MarkerFacecolor','b')
end
axis([0 3 0 3])
set(gca,'xtick',0:3,'ytick',0:3)
xlabel('$\xi_1$','Interpreter','Latex','HorizontalAlignment','Center')
ylabel('$\xi_2$','Interpreter','Latex','HorizontalAlignment','Center')
title('PC basis','Fontname','TimesNewRoman')
if print_to_eps=='y';
     print(pictype,resolution,[write_dir,'searchspace']),close;
     result = eps2xxx([write_dir,'searchspace.eps'],{'pdf'},'C:\Program Files\gs\gs9.05\bin\gswin64c.exe');
end