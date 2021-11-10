function []=fadigatesterede_EIW()

% Copia o diretorio corrente
diretorio_corrente=cd;

% Definindo os valores de Tmed Normalizado
tmed_normalizada=0:0.01:1;
tmed_normalizada=[tmed_normalizada]';

% Obtendo os valores de tensao ultima a tracao e a compressao
% Entrando com o valores
mensagem={' Digite o numero de ciclos a ser avaliado log10(N): '...
        ' Digite o peso linear analisado'
        ' Digite o Valor Maximo de Numero de ciclos log10(N): '...
        'Digite o valor maximo de EI'};
titulo=' Limite de Resistencia Estatico dos Cabos';
num_de_linhas=[1 1 1 1]';
valor_predefinido={'7' '800' '8' '31000000'};
valores=inputdlg(mensagem,titulo,num_de_linhas,valor_predefinido);
if (isempty(valores))
    warndlg('Nao se forneceu nenhum valor.');
    
    % Volta para o diretorio corrente
    cd (diretorio_corrente)
    return;
end
tult=300; % O Valor aqui definido não vai influenciar os resultados obtidos, ja que se deseja valores normalizados
ciclo_aval=str2num(char(valores(1))); % Numero de Ciclos Avaliado
pl_aval=str2num(char(valores(2)));
maxciclos=str2num(char(valores(3)));  % Valor maximo de numero de ciclos 'lembrar que ciclos=log10(N)'
EI_max=str2num(char(valores(4)));     % EI maximo
pl_normalizado=pl_aval/1500;

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
EI_normalizado=0.4:0.001:1;
EI_normalizado=EI_normalizado';
for i=1:length(EI_normalizado)
    [ta_normalizado(:,i), ta(:,i), tmed(:,i)]=curva(ciclo_aval, tmed_normalizada,...
        tult, maxciclos, EI_normalizado(i,1), W, V, tult, pl_normalizado);
end

% Para 10000 GPa.mm4
EI_normalizada=15000000/EI_max;
[ta_normalizado3, ta3, tmed3] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, EI_normalizada, W, V, tult, pl_normalizado);

% Para 15000 GPa.mm4
EI_normalizada=20000000/EI_max;
[ta_normalizado4, ta4, tmed4] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, EI_normalizada, W, V, tult, pl_normalizado);

% Para 20000 GPa.mm4
EI_normalizada=25000000/EI_max;
[ta_normalizado5, ta5, tmed5] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, EI_normalizada, W, V, tult, pl_normalizado);

% Para 25000 GPa.mm4
EI_normalizada=30000000/EI_max;
[ta_normalizado6, ta6, tmed6] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, EI_normalizada, W, V, tult, pl_normalizado);

% Para 30000 GPa.mm4
EI_normalizada=35000000/EI_max;
[ta_normalizado7, ta7, tmed7] = curva (ciclo_aval, tmed_normalizada,...
    tult, maxciclos, EI_normalizada, W, V, tult, pl_normalizado);

% Saida Grafica
figure
plot(tmed_normalizada,ta_normalizado3,'-k',tmed_normalizada,ta_normalizado4,'-r',...
    tmed_normalizada,ta_normalizado5,'-c',tmed_normalizada,ta_normalizado6,'-m',tmed_normalizada,ta_normalizado7)
legend('EI=15000 GPa.mm4','EI=20000 GPa.mm4', 'EI=25000 GPa.mm4','EI=30000 GPa.mm4','EI=35000 GPa.mm4')
xlabel('Tensao Media/Tensao Ultima')
ylabel('Amplitude de Tensao/Tensao Ultima')

figure
[EI_NOR, TMED]=meshgrid(EI_normalizado,tmed_normalizada);
mesh(TMED,EI_NOR.*EI_max./(10.^6),ta_normalizado)
xlabel('Tensao Media (MPa)')
ylabel('EI (N.m^2)')
zlabel('Amplitude de Tensao (MPa)')

% Salvando dados em arquivo
tipo={'*.dat'; '*.txt'};
titulo=' Arquivando os dados da Rede para a Curva de vida constante';
[nome,Caminho]=uiputfile(tipo,titulo);
if ~ischar(nome)
    warndlg('Nome de arquivo nao fornecido.');
    
else
    % Criando o arquivo
    Curva_VC=[tmed_normalizada ta_normalizado3 ta_normalizado4 ta_normalizado5 ta_normalizado6 ta_normalizado7];
    cd (Caminho)
    dlmwrite(nome,Curva_VC,'\t')
end

% Volta para o diretorio corrente
cd (diretorio_corrente)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ta_normalizado, ta, tmed] = curva (ciclo, tmed_normalizada,...
    maxta, maxciclos, peso_linear_normalizado, W, V, maxtmed, pl_normalizado)

% Definindo o numero de ciclos 'lembrar que ciclos=log10(N)' 
ciclos=ciclo.*ones(size(tmed_normalizada));
pl_normalizados=pl_normalizado.*ones(size(tmed_normalizada));
pesos_nor=peso_linear_normalizado.*ones(size(tmed_normalizada));

% Dados para serem utilizados na rede
xteste=[(ciclos./maxciclos) tmed_normalizada pesos_nor pl_normalizados]; 

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
