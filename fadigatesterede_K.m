function []=fadigatesterede_K()

% Copia o diretorio corrente
diretorio_corrente=cd;

% Definindo os valores de Tmed Normalizado
K_normalizada=0.01:0.01:1;
K_normalizada=[K_normalizada]';

maxciclos=8;  % Valor maximo de numero de ciclos 'lembrar que ciclos=log10(N)'

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
ciclos=5:0.1:8;
ciclos=ciclos';
for i=1:length(ciclos)
    ta_normalizado(:,i)=curva(ciclos(i,1), K_normalizada, maxciclos, W, V);
end

% Para 10^5
ciclo=5;
ta_normalizado3=curva(ciclo, K_normalizada, maxciclos, W, V);

% Para 10^6
ciclo=6;
ta_normalizado4= curva(ciclo, K_normalizada, maxciclos, W, V);

% Para 10^7
ciclo=7;
ta_normalizado5=curva(ciclo, K_normalizada, maxciclos, W, V);

% Para 10^8
ciclo=8;
ta_normalizado6=curva(ciclo, K_normalizada, maxciclos, W, V);

% Obtendo os dados de Wfixo
[arq, Caminho] = uigetfile({'*.*'}, 'Escolha o arquivo com todos os dados' );
if ~ischar(arq)  % Verifica se algum dado foi fornecido
    warndlg('Nome de arquivo nao fornecido.');
    return;
end
Arquivo=fullfile(Caminho,arq); % Compoe nome do arquivo
cd (Caminho)
dad=dlmread(arq);

% Retirando os dados do arquivo todos
local=find(dad(:,1)==5);
dad3=[dad(local,5)./dad(local,6) dad(local,3)./dad(local,4)];
local=find(dad(:,1)==6);
dad4=[dad(local,5)./dad(local,6) dad(local,3)./dad(local,4)];
local=find(dad(:,1)==7);
dad5=[dad(local,5)./dad(local,6) dad(local,3)./dad(local,4)];
local=find(dad(:,1)==8);
dad6=[dad(local,5)./dad(local,6) dad(local,3)./dad(local,4)];

% Saida Grafica
figure
plot(K_normalizada,ta_normalizado3,'-k',K_normalizada,ta_normalizado4,'-r',...
    K_normalizada,ta_normalizado5,'-c',K_normalizada,ta_normalizado6,'-m',dad3(:,1),dad3(:,2),'k<',...
    dad4(:,1),dad4(:,2),'rd',dad5(:,1),dad5(:,2),'cp',dad6(:,1),dad6(:,2),'mh')
legend('N=10^5','N=10^6','N=10^7','N=10^8')
xlabel('K/K_{max}')
ylabel('Amplitude de Tensao/Tensao Ultima')

figure
[N_NOR, TMED]=meshgrid(ciclos,K_normalizada);
mesh(TMED,N_NOR,ta_normalizado)
xlabel('K/K_{max}')
ylabel('Numero de ciclos (log10(N))')
zlabel('Amplitude de Tensao/Tensao Ultima')

% Salvando dados em arquivo
tipo={'*.dat'; '*.txt'};
titulo=' Arquivando os dados da Rede para a Curva de vida constante';
[nome,Caminho]=uiputfile(tipo,titulo);
if ~ischar(nome)
    warndlg('Nome de arquivo nao fornecido.');
    
else
    % Criando o arquivo
    Curva_VC=[K_normalizada ta_normalizado3 ta_normalizado4 ta_normalizado5 ta_normalizado6];
    
    cd (Caminho)
    dlmwrite(nome,Curva_VC,'\t')
end

% Volta para o diretorio corrente
cd (diretorio_corrente)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ta_normalizado] = curva(ciclo, K_normalizado,...
    maxciclos, W, V)

% Definindo o numero de ciclos 'lembrar que ciclos=log10(N)' 
ciclos=ciclo.*ones(size(K_normalizado));

% Dados para serem utilizados na rede
xteste=[(ciclos./maxciclos) K_normalizado]; 

% TESTE DA REDE MLP
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

% Corrigindo erros
exemplos=length(ta_normalizado);
i=1;
while i<=exemplos
    % Se a Amplitude de Tensao for menor que zero (valor absurdo)
    if (ta_normalizado(i,1)<0)
        ta_normalizado(i,1)=0;
    end
    
    % Contador
    i=i+1;
end
