function [ta, K, Kmax, EI] = Sa_cabo_AL_KW(log_N, tmed, tult, pl, A, d, n)
% This Function Calculate Sa from an aluminum cable used KW Model
%
% Na entrada, log_N é o número de ciclos em escala log10 desejado, tmed é 
% a tensão média aplicada ao cabo (MPa), tult é a tensão ultima suportada 
% pelo cabo (MPa), pl é o seu peso linear (kg/m), E é o modulo de 
% elasticidade (GPa), A é a área total do cabo (mm^2), d é o diamentro do fio
% (mm), e n é o número de fios do cabo e as saídas são ta o valor da amplitude
% de tensão (MPa), K, Kmax e EI com o valor de rigidez a flexão (N*m^2)

% Valores fixos
maxciclos=8; % Valor maximo de numero de ciclos 'lembrar que ciclos=log10(N)'
x= 89;       % ultimo ponto de contato (mm)
pl_max=1500; % peso linear maximo (kg/km)
E=69000;     % Modulo de Elasticidade do Al (MPa)

% Calculo de EI e pl
EI=(n.*E.*pi.*d.^4)/64;
pl=pl.*1000;

% Calculo de K
p=((tmed.*A)/EI).^(1/2);
px=p.*x;
e_px=exp(-p.*x);         
K=(E.*d.*(p.^2))./(4.*(e_px-1+px));

% Calculo de Kmax
p_max=((tult.*A)/EI).^(1/2);
px_max=p_max.*x;
e_px_max=exp(-p_max.*x);         
Kmax=(E.*d.*(p_max.^2))./(4.*(e_px_max-1+px_max));

% Normalizando os valores
K_nor=K./Kmax;
pl_nor=pl./pl_max;
log_N_nor=log_N./maxciclos;

% Verificacao de erros
if K>Kmax
    warndlg('O valor de tmed deve ser inferior a tult','Erro em tmed');
    return;
end
if pl>pl_max
    warndlg('O valor de pl deve ser inferior a 1500 Kg/m','Erro em pl');
end
if log_N>maxciclos
    warndlg('O valor de log(N) deve ser inferior a 8','Erro no log(N)');
end


% ANN Weights
W=[1.1967 -0.85809 -0.86429 -0.90482;-1.1915 -6.5330 2.2779 -2.6607;0.81309 -1.4217 -1.2068 -1.1556;1.0539 -1.0883 -1.1251 -0.90488;0.99379 -0.95201 -0.85212 -1.3963;1.2380 -0.63754 -0.94961 -0.92551;1.1419 -0.66303 -0.83848 -1.2504;1.2429 -0.61019 -0.80194 -1.0303;1.3201 -0.60186 -0.69994 -1.0311;1.2652 -0.83176 -0.70433 -0.93348;1.2138 -0.56096 -0.89864 -0.98467;0.73997 -1.5253 -0.49134 -2.4110;1.1940 -0.76818 -1.0123 -0.89452;1.3276 -0.69328 -0.74764 -0.91612;1.0242 -0.97196 -1.0513 -1.0430;1.2908 -0.72147 -0.82335 -0.89874;-1.0618 -3.3374 -2.7149 -1.1603;1.2732 -0.63704 -0.88695 -0.91948;0.70619 -1.2231 -0.96024 -1.5587;0.80471 -1.1810 -0.76842 -1.8407;1.2880 -0.48282 -0.82233 -0.99193;1.1618 -0.82056 -0.66815 -1.1208;1.1729 -0.85986 -0.98793 -0.87397;1.1574 -0.89129 -0.68408 -1.4333;1.2979 -0.96758 -0.95386 -0.81295;1.1757 -0.79961 -1.0042 -0.88511;1.2154 -0.89246 -0.91559 -0.86867;1.0622 -0.91046 -0.79631 -1.4551;1.2145 -0.81131 -0.67242 -1.3373;1.3504 -1.1459 -1.0580 -0.75424;1.1275 -0.78617 -0.96429 -0.96393;0.99894 -0.78925 -0.97050 -1.1698;1.3172 -0.82859 -0.63108 -0.93290;1.3624 -1.0671 -0.84647 -0.78979;1.1903 -1.2135 -0.95481 -0.82623;-6.3554 -5.1468 -1.6639 -5.0719;1.1165 -0.69857 -0.88788 -1.1026;1.3101 -0.22654 -0.92099 -1.0867;1.0205 -1.4855 -1.3024 -0.84939;1.3484 -0.52063 -0.79672 -0.94485;1.0239 -0.92076 -0.79947 -1.4088;1.0881 -0.70127 -0.93525 -1.0367;1.1469 -0.77339 -0.94626 -0.94405;1.2027 -0.70532 -0.86500 -0.95574;1.0748 -0.85685 -0.95445 -1.0056;1.1591 -0.76336 -0.77232 -1.0521;1.1516 -0.71628 -0.85289 -1.0109;1.3359 -0.92285 -1.0672 -0.79870;1.3981 -0.67289 -1.0253 -0.85832;1.2705 -0.62325 -0.87246 -0.92891;1.9170 -1.6804 3.6612 -4.9806;1.2163 -0.68284 -0.79934 -0.99329;1.2976 -0.57811 -0.91636 -0.91978;1.4197 -0.73260 -0.79985 -0.83164;1.0915 -0.75093 -0.92111 -1.1312;1.0140 -0.84386 -1.0104 -1.0169;1.1192 -0.77309 -0.92950 -0.99563;1.2815 -0.65695 -0.87776 -0.90997;1.1542 -0.81919 -0.74055 -1.2878;1.2805 -0.74434 -0.64220 -1.2418;1.2713 -0.67029 -0.80153 -0.94157;1.1236 -0.76230 -0.83288 -1.0281;1.0948 -0.73187 -0.94784 -0.99302;-1.6238 -3.6825 0.38621 -5.2077;1.1515 -0.82429 -0.67365 -1.1932;1.0485 -0.84902 -0.83695 -1.4483;1.0176 -1.0278 -0.74596 -1.7333;1.1408 -0.75383 -0.80735 -1.1990;1.1829 -0.68239 -0.82859 -1.0189;1.2268 -0.64409 -0.75424 -1.0761;1.4918 -0.93587 -0.59513 -0.82296;1.2404 -0.64984 -0.82851 -0.95917;1.4035 -0.59041 -0.62100 -0.97936;1.1601 -0.84260 -0.71079 -1.3982;1.1327 -0.79850 -0.80057 -1.0217;1.2542 -0.84336 -0.88325 -0.86627;1.2897 -0.63084 -0.79299 -0.94653;1.6424 -1.9091 -0.45698 -0.75824;1.3259 -0.56774 -0.77157 -0.96902;1.0357 -0.85837 -0.82249 -1.4466;1.2803 -0.57507 -0.73668 -1.0670;1.4170 -0.88360 -0.74299 -0.81349;1.1512 -0.85123 -0.96758 -0.89704;1.2837 -0.62272 -0.79572 -0.95563;1.2965 -0.59852 -0.84756 -0.93057;0.59672 -7.0929 -3.4473 5.9746;1.0618 -0.87755 -1.0278 -0.94480;1.2640 -0.51449 -0.84705 -0.98832;0.89391 -1.3231 -0.55259 -2.0931;1.1068 -0.79711 -0.80286 -1.1879;0.53051 -1.6330 -0.48272 -2.5023;1.3325 -0.89735 -1.0262 -0.80753;1.2320 -0.81917 -1.1325 -0.87877;1.2783 -1.2176 -0.86626 -0.80631;1.2864 -0.69151 -0.82426 -0.91163];
V=[-0.062510;-0.042501;-2.7267;-0.052732;-0.049972;0.30090;-0.14944;0.052052;-0.045457;-0.046629;-0.025493;-0.11727;1.1193;-0.14276;-0.073538;0.085671;-0.095354;-2.7027;-0.13543;0.33592;0.59580;-0.17341;0.083260;-0.084395;0.28505;-0.19139;-0.095231;-0.11178;0.30798;0.19552;-0.32550;0.0082690;0.096679;-0.015858;-0.18581;-0.22251;1.6286;0.028806;-0.39105;-0.22125;-0.15682;0.27761;-0.016366;-0.019946;-0.046897;0.041497;0.023456;-0.021692;-0.24694;-0.25104;-0.12302;-1.6912;-0.041042;-0.16347;-0.13294;0.084024;0.020838;0.014220;-0.13113;0.18283;0.12432;-0.074277;-0.0015771;-0.032951;3.6212;0.12194;0.24531;0.47493;0.10405;-0.026330;-0.020385;-0.071584;-0.085352;-0.080516;0.23125;0.010423;-0.10448;-0.096528;-0.37727;-0.088315;0.24006;-0.059023;-0.10514;-0.064381;-0.089416;-0.12894;1.4731;-0.016956;-0.14311;0.82629;0.11428;1.1877;-0.22473;-0.23326;-0.21221;-0.094646];

% Dados para serem utilizados na rede
xteste=[log_N_nor K_nor pl_nor]; 

%%%%%%% ANN %%%%%%%
% Camada de entrada
x=[-1 xteste]';
nety=W*x;

% CAMADA OCULTA
y=(1)./(1+exp(-nety'));
netz=[-1 y]*V;

% CAMADA DE SAIDA
ta_nor=netz';
ta=(ta_nor).*tult;

% Saida de EI
EI=EI./1000000;

%%% Correcao de erros %%%
if (ta_nor<0)
    K=Kmax;
    ta=0;
elseif ta>tult
   ta='NAN';
   warndlg('O valor de Sa nao e confiavel','Erro em Sa');
end


