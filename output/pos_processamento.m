% _______________________________________________________________________________________________________
% 
%
% Algoritmo para pós-processamento dos módulos hidrodinâmico e de qualidade de água do programa SihQual 
% Autora: Danieli Mara Ferreira | contato: danielimaraferreira@gmail.com
% Versão: Tutorial simplificado (10/03/2021) | Projeto Enquadramento Rio Paranapanema - ANA/UFPR
%
% _______________________________________________________________________________________________________

clear all; close all; clc

abs='C:\SIHQUAL\input\Dados_Entrada.xlsx';
data1 = xlsread(abs,1);    % vazões diárias (2012) em estações fluviométricas (informação hidroweb)
SE2=data1(:,2);            % 64278080	
data3 = xlsread(abs,3);    % concentrações observadas em período definido na seção de referência (informação hidroweb/cetesb)
conc_obs=data3(:,1);    

figure(1)
fileID = fopen('C:\SIHQUAL\simulacao\res_hidrodinamico.txt','r');
formatSpec = '%f';
At1 = fscanf(fileID,formatSpec);
fclose(fileID);
cp3=At1(1:(length(At1)-1),1);
cp3r=reshape(cp3,[86400/20,365]);
At1=mean(cp3r);
At1(1)=At1(2);
At1(2:366)=At1(1:365);
plot(At1)
title('64278080')
hold on; plot(SE2)
legend('Simulado','Observado','Orientation','horizontal','location','southeast')
startDate = datenum('01/01/2012','dd/mm/yyyy');
endDate = datenum('31/12/2012','dd/mm/yyyy');
xdata = linspace(startDate,endDate,5);
set(gca,'XTick',xdata);
datetick('x','m'); 
ylabel('\itQ (m^3/s)')
xlabel('{\itData}')
set(gcf, 'PaperUnits', 'inches');
x_width=4 ;y_width=1.8;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
set(gca ,'FontWeight','normal','FontSize',7,'FontName','Times')
print('hidrograma-64278080', '-dpng', '-r800');

fileID = fopen('C:\SIHQUAL\simulacao\res_qualidade.txt','r');
formatSpec = '%f';
At1 = fscanf(fileID,formatSpec);
fclose(fileID);

figure(2)
group = [ones(size(At1));ones(size(conc_obs))+1];
x=cat(1,At1,conc_obs);
boxplot(conc_obs)
h=boxplot(x,group,'symbol','.','labels',{'Simulado','Observado'},'color',[0 0 0]);%
title('64326000/PARP 02500')
ylabel('\itC_P_T (mg-P/L)')
box('on')
ylim([0 0.15])
set(gcf, 'PaperUnits', 'inches');
x_width=4 ;y_width=1.8;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
title('64326000/PARP 02500')
set(gca ,'FontWeight','normal','FontSize',7,'FontName','Times')
print('boxplot', '-dpng', '-r800'); 

figure(3)
D=sort(At1,'descend');    
cc=1:1:length(At1); 
N=length(cc);
fre=(cc/(N+1))*100;   
plot(fre,D); 
D=sort(conc_obs,'descend');    
cc=1:1:length(conc_obs); 
N=length(cc);
fre=(cc/(N+1))*100;   
hold on; plot(fre,D);
ylabel('\itC_P_T (mg-P/L)')
xlabel('\itFrequência (%)\rm')
box('on')
ylim([0 0.15])
x_width=4 ;y_width=1.8;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
legend('Simulado','Observado','Orientation','horizontal')
title('64326000/PARP 02500')
set(gca ,'FontWeight','normal','FontSize',7,'FontName','Times')
print('curva_permanencia', '-dpng', '-r800'); 
