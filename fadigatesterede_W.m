function []=fadigatesterede_W()

% Copia o diretorio corrente
diretorio_corrente=cd;

% Definindo os valores de Tmed Normalizado
tmed_normalizada=0:0.01:1;
tmed_normalizada=[tmed_normalizada]';

% Obtendo os valores de tensao ultima a tracao e a compressao
% Entrando com o valores
mensagem={' Digite o numero de ciclos a ser avaliado log10(N): '...
        ' Digite o Valor Maximo de Numero de ciclos log10(N): '...
        'Digite o valor maximo de W'};
titulo=' Limite de Resistencia Estatico dos Cabos';
num_de_linhas=[1 1 1]';
valor_predefinido={'6' '8' '1500'};
valores=inputdlg(mensagem,titulo,num_de_linhas,valor_predefinido);
if (isempty(valores))
    warndlg('Nao se forneceu nenhum valor.');
    
    % Volta para o diretorio corrente
    cd (diretorio_corrente)
    return;
end
tult=300; % O Valor aqui definido não vai influenciar os resultados obtidos, ja que se deseja valores normalizados
ciclo_aval=str2num(char(valores(1))); % Numero de Ciclos Avaliado
maxciclos=str2num(char(valores(2)));  % Valor maximo de numero de ciclos 'lembrar que ciclos=log10(N)'
W_max=str2num(char(valores(3)));     % EI maximo

% Obtendo a matriz de W em um arquivo
[arq, Caminho] = uigetfile({'*.*'}, 'Escolha o arquivo da matriz de pesos W' );
if ~ischar(arq)  % Verifica se algum dado foi fornecido
    warndlg('Nome de arquivo nao fornecido.');
    return;
end
Arquivo=fullfile(Caminho,arq); % Compoe nome do arquivo
cd (Caminho)
W=dlmread(arq,'\t');
tamanho_de_W=size(W);

% Obtendo a matriz de V em um arquivo
arq=['V' arq(2:length(arq))];
V=dlmread(arq,'\t');
V=V';

% Volta para o diretorio corrente
cd (diretorio_corrente)

% PESOS DA REDE
Ni=4;                 % Neuronios de entrada
Ns=1;                 % Neuronios de saida
Nh=tamanho_de_W(1,1); % Neuronios ocultos

% Criando Matriz tridimensional
W_normalizado=0.4:0.001:1;
W_normalizado=W_normalizado';
for i=1:length(W_normalizado)
    [ta_normalizado(:,i), ta(:,i), tmed(:,i)]=curva(ciclo_aval, tmed_normalizada,...
        tult, maxciclos, W_normalizado(i,1), W, V, tult);
end

% Para 700 Kg/Km
W_normalizada=700/W_max;
[ta_normalizado8, ta3, tmed3] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, W_normalizada, W, V, tult);

% Para 900 Kg/Km
W_normalizada=900/W_max;
[ta_normalizado4, ta4, tmed4] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, W_normalizada, W, V, tult);

% Para 1100 Kg/Km
W_normalizada=1100/W_max;
[ta_normalizado5, ta5, tmed5] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, W_normalizada, W, V, tult);

% Para 1300 Kg/Km
W_normalizada=1300/W_max;
[ta_normalizado6, ta6, tmed6] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, W_normalizada, W, V, tult);

% Para 1500 Kg/Km
W_normalizada=1500/W_max;
[ta_normalizado7, ta7, tmed7] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, W_normalizada, W, V, tult);

% Para 889 Kg/Km
W_normalizada=889/W_max;
[ta_normalizado8, ta8, tmed8] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, W_normalizada, W, V, tult);

% Para 1046 Kg/Km
W_normalizada=1046/W_max;
[ta_normalizado9, ta9, tmed9] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, W_normalizada, W, V, tult);

% Para 1252 Kg/Km
W_normalizada=1252/W_max;
[ta_normalizado10, ta10, tmed10] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, W_normalizada, W, V, tult);

% Para 1466 Kg/Km
W_normalizada=1466/W_max;
[ta_normalizado11, ta11, tmed11] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, W_normalizada, W, V, tult);


% Obtendo os dados de Wfixo
[arq, Caminho] = uigetfile({'*.*'}, 'Escolha o arquivo com todos os dados' );
Arquivo=fullfile(Caminho,arq); % Compoe nome do arquivo
cd (Caminho)
dad=dlmread(arq);

% Retirando os dados do arquivo todos
local=find(dad(:,1)==ciclo_aval);
dad3=[dad(local,2)./dad(local,4) dad(local,3)./dad(local,4) dad(local,8)];

% Organizando os dados de comparacao dos cabos
local=find(dad3(:,3)==889);
orchid=dad3(local,1:2);
local=find(dad3(:,3)==1046);
acar=dad3(local,1:2);
local=find(dad3(:,3)==1466);
cal=dad3(local,1:2);
local=find(dad3(:,3)==1252);
aaac=dad3(local,1:2);


% Saida Grafica
figure
plot(tmed_normalizada,ta_normalizado8,'-k',tmed_normalizada,ta_normalizado4,'-r',...
    tmed_normalizada,ta_normalizado5,'-c',tmed_normalizada,ta_normalizado6,'-m',tmed_normalizada,ta_normalizado7)
legend('W=700 Kg/Km','W=900 Kg/Km', 'W=1100 Kg/Km','W=1300 Kg/Km','W=1500 Kg/Km')
xlabel('Tensao Media/Tensao Ultima')
ylabel('Amplitude de Tensao/Tensao Ultima')

figure
plot(tmed_normalizada,ta_normalizado8,'-k',tmed_normalizada,ta_normalizado9,'-r',...
    tmed_normalizada,ta_normalizado10,'-c',tmed_normalizada,ta_normalizado11,'-m',orchid(:,1), orchid(:,2),'ok',...
    acar(:,1), acar(:,2),'xr',aaac(:,1), aaac(:,2),'oc',cal(:,1), cal(:,2),'xm')
legend('W=889 Kg/Km','W=1046 Kg/Km', 'W=1252 Kg/Km','W=1466 Kg/Km')
xlabel('Tensao Media/Tensao Ultima')
ylabel('Amplitude de Tensao/Tensao Ultima')

figure
[W_NOR, TMED]=meshgrid(W_normalizado,tmed_normalizada);
mesh(TMED,W_NOR,ta_normalizado)
xlabel('Tensao Media/Tensao Ultima')
ylabel('Peso Linear/Peso Linear Maximo')
zlabel('Amplitude de Tensao/Tensao Ultima')

% Salvando dados em arquivo
tipo={'*.dat'; '*.txt'};
titulo=' Arquivando os dados da Rede para a Curva de vida constante';
[nome,Caminho]=uiputfile(tipo,titulo);
if ~ischar(nome)
    warndlg('Nome de arquivo nao fornecido.');
    
else
    % Segundo arquivo
    nome_comp=['Comparacao_' nome];
    
    % Criando o arquivo
    Curva_VC=[tmed_normalizada ta_normalizado3 ta_normalizado4 ta_normalizado5 ta_normalizado6 ta_normalizado7];
    Curva_VC2=[tmed_normalizada ta_normalizado8 ta_normalizado9 ta_normalizado10 ta_normalizado11];
    cd (Caminho)
    dlmwrite(nome,Curva_VC,'\t')
    dlmwrite(nome_comp,Curva_VC2,'\t')
end

% Volta para o diretorio corrente
cd (diretorio_corrente)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ta_normalizado, ta, tmed] = curva (ciclo, tmed_normalizada,...
    maxta, maxciclos, peso_linear_normalizado, W, V, maxtmed)

% Definindo o numero de ciclos 'lembrar que ciclos=log10(N)' 
ciclos=ciclo.*ones(size(tmed_normalizada));

pesos_nor=peso_linear_normalizado.*ones(size(tmed_normalizada));

% Dados para serem utilizados na rede
xteste=[(ciclos./maxciclos) tmed_normalizada pesos_nor]; 

% TESTE DA REDE MLP

% COMPUTACAO NO SENTIDO DIRETO
for n=1:length(xteste)
    % Camada de entrada
    x=[-1 xteste(n,:)]';
    nety=W*x;
    
    % CAMADA OCULTA
    y=(1)./(1+exp(-nety'));
    netz=V*[-1 y]';
    
    % CAMADA DE SAIDA
    zteste(n,:)=netz';
end

% Saida dos resultados de Amplitude de tensao normalizado
ta_normalizado=zteste;

% Desnormalizando os valores da amplitude de tensao
ta=(ta_normalizado).*maxta;
tmed=tmed_normalizada.*maxtmed;

% Corrigindo erros
exemplos=length(tmed_normalizada);
i=1;
while i<=exemplos
    % Se a Amplitude de Tensao for menor que zero (valor absurdo)
    if (ta_normalizado(i,1)<0)
        ta_normalizado(i,1)=0;
        tmed(i,1)=maxtmed;
        ta(i,1)=0;
    end
    
    % Contador
    i=i+1;
end
