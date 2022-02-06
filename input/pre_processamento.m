% _______________________________________________________________________________________________________
% 
%
% Algoritmo para pré-processamento dos módulos hidrodinâmico e de qualidade de água do programa SihQual 
% Autora: Danieli Mara Ferreira | contato: danielimaraferreira@gmail.com
% Versão: Tutorial simplificado (10/03/2021) | Projeto Enquadramento Rio Paranapanema - ANA/UFPR
%
% _______________________________________________________________________________________________________

clear all; close all; clc

dx=500;          	   % -> passo espacial
dt=20;           	   % -> passo temporal
LL=53556.6475;   	   % -> comprimento do trecho simulado  
xx=0:dx:LL;      	   % -> vetor de discretização espacial dx
J=length(xx);
t=0:dt:31536000; 	   % -> vetor de discretização temporal dt
ta=0:86400:31536000;       % -> vetor de discretização temporal diária   
N=length(t);
alfa=0.1;                  % -> parâmetro do método numérico
g=9.81;                    % -> aceleração da gravidade

% distância em relação ao início do trecho entre as estações fluviométricas de interesse (determinada a partir do shapefile do rio):
dist_fluv=[0	17592.9822	53556.6475];

% posições no vetor de discretização onde estão as estações:
TMP = bsxfun(@(x,y) abs(x-y), dist_fluv(:), reshape(xx,1,[]));
[D, idxB] = min(TMP,[],2) ;
dist_fluv_discretizacao=xx(idxB); % -> dist_fluv na discretização do modelo
for jj=1:length(dist_fluv_discretizacao)
    posicao_fluv(jj)=find(xx==dist_fluv_discretizacao(jj));
end

% dados de interesse:
abs='Dados_Entrada.xlsx';
data1 = xlsread(abs,1);    % vazões diárias (2012) em estações fluviométricas (informação monitoramento)
SE1=data1(:,1);            % estação 64270080		
SE2=data1(:,2);            % estação 64278080	
SE3=data1(:,3);            % estação 64332080	
data2 = xlsread(abs,2);    % concentrações diárias (2012) no contorno de montante (informação gerada após pré-processamento)
conc_montante=data2(:,1);            	
data4 = xlsread(abs,4);    % cargas aportadas em cada nó computacional (informação gerada por processamento dos resultados do Modelo de Bacias)
carga_bacia=data4(:,2); 

% declividade do fundo (calibrada a partir de informações de altitude e distância entre as estações fluviométricas)
So(1:J)=0.000470;
fid = fopen('so.txt','wt');
fprintf(fid, '%f\n', [So']); 
fclose(fid);
movefile('so.txt','C:\SIHQUAL\simulacao');

% largura de fundo (calibrada a partir de informações de monitoramento)
b1(1:J)=60;
fid = fopen('bb.txt','wt');
fprintf(fid, '%f\n', [b1']); 
fclose(fid);
movefile('bb.txt','C:\SIHQUAL\simulacao');

% inclinação do talude do canal (vem de informação monitoramento + pré-processamento)
mm(1:J)=interp1([0 LL],[4  4.5],[xx(1):dx:LL],'pchip');  
fid = fopen('mm.txt','wt');
fprintf(fid, '%f\n', [mm']); 
fclose(fid);
movefile('mm.txt','C:\SIHQUAL\simulacao');

n(1:J)=0.02;
fid = fopen('manning.txt','wt');
fprintf(fid, '%f\n', [n']); 
fclose(fid);
movefile('manning.txt','C:\SIHQUAL\simulacao');         	

% cálculo de contribuição lateral:
qq=zeros(length(ta),J);

q1=(SE2-SE1)/(dist_fluv(2)-dist_fluv(1));
q2=(SE3-SE2)/(dist_fluv(3)-dist_fluv(2));

qq(:, posicao_fluv(1):posicao_fluv(2))=repmat(q1,1,length(posicao_fluv(1):posicao_fluv(2)));
qq(:, posicao_fluv(2)+1:posicao_fluv(3))=repmat(q2,1,length(posicao_fluv(2)+1:posicao_fluv(3)));

% interpolações em dt 
X = qq';
    for i = 1:size(X,1)  % linha
        new_X(i,:) = interp1(ta,X(i,:)',t, 'pchip');
    end
q=new_X';
fid = fopen('hm_lateralcontrib.txt','wt');
for ii = 1:size(q,1)
    fprintf(fid,'%g\t',q(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
movefile('hm_lateralcontrib.txt','C:\SIHQUAL\simulacao');

% vazao - contorno de montante:
Qaam=interp1(ta,SE1,t,'pchip');   
fid = fopen('vazao-m.txt','wt');
fprintf(fid, '%f\n', [Qaam]); 
fclose(fid);
movefile('vazao-m.txt','C:\SIHQUAL\simulacao');

% profundidade - contorno de montante:
cota_m=(-0.4643+(Qaam./20.0).^(1/2.5229)); % eq. ajustada da curva de descarga da estação 64270050 (selecionada devido à limitação de dados no contorno)
fid = fopen('cota-m.txt','wt');
fprintf(fid, '%f\n', [cota_m]); 
fclose(fid);
movefile('cota-m.txt','C:\SIHQUAL\simulacao');

% condições iniciais:
y1(1)=cota_m(1);      
y1(J)=1.5*y1(1);                           % hipótese assumida
y1(2:(J-1))=interp1([1  J],[y1(1)  y1(J)],(2:(J-1)),'pchip');
fid = fopen('prof.txt','wt');
fprintf(fid, '%f\n', [y1(1:J)]); 
fclose(fid);
movefile('prof.txt','C:\SIHQUAL\simulacao');

A1(1:J)=b1(1:J).*y1(1:J)+mm(1:J).*y1(1:J).^2;
v1(1)=Qaam(1)/A1(1);
v1(J)=1.5*v1(1);                          % hipótese assumida
v1(2:(J-1))=interp1([1  J],[v1(1)  v1(J)],(2:(J-1)),'pchip');
fid = fopen('vel.txt','wt');
fprintf(fid, '%f\n', [v1(1:J)]); 
fclose(fid);
movefile('vel.txt','C:\SIHQUAL\simulacao');

% carga difusa e pontual (vem do Modelo de Bacias + processamento para adequar escala):
fid = fopen('carga.txt','wt');
fprintf(fid, '%f\n', [carga_bacia]); 
fclose(fid);
movefile('carga.txt','C:\SIHQUAL\simulacao');

% concentração - contorno de montante:
centr=interp1(ta,conc_montante,t,'pchip');  
fid = fopen('upboundary.txt','wt');
fprintf(fid, '%f\n', [centr']); 
fclose(fid);
movefile('upboundary.txt','C:\SIHQUAL\simulacao');

% concentração inicial (assumida):
fid = fopen('initial_c.txt','wt');
fprintf(fid, '%f\n', [0.1/1000*ones(length(1:J),1)]); 
fclose(fid);
movefile('initial_c.txt','C:\SIHQUAL\simulacao');


