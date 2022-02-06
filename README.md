# SihQual-Tutorial
 
#### Passos 

Fontes dos dados:

1. criar uma pasta chamada 'SIHQUAL' (em letras maiúsculas) no disco  'C', de modo que o caminho de trabalho seja 'C:\SIHQUAL’
2. fazer download da pasta "input" e salvar no diretório criado em 1.  ('C:\SIHQUAL\input’)
3. fazer download da pasta "simulacao" e salvar no diretório criado em 1. ('C:\SIHQUAL\simulacao’).
4. fazer download da pasta "output" e salvar no diretório criado em 1. ('C:\SIHQUAL\output’). 
5. No programa Matlab (ou Octave), executar o script “pre_processamento.m”, que está na pasta 'C:\SIHQUAL\input.
6. Executar o programa principal abrindo o arquivo “sihqual.exe”, que está na pasta 'C:\SIHQUAL\simulação’ (uma janela será aberta e fechará automaticamente, indicando o fim da simulação)  
7. No programa Matlab (ou Octave), executar o script “pos_processamento.m”, que está na pasta 'C:\SIHQUAL\output. Nessa mesma pasta serão salvas as figuras com resultados na estação de referência (os resultados são gerados a cada passo temporal e espacial do modelo, porém impressos somente em estações definidas previamente, onde há dados de monitoramento para comparação: para o caso de interesse, são a estação fluviométrica de código 64278080 e a seção de qualidade de água 64326000/PARP 02500)
