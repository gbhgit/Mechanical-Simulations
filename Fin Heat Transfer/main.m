% INSTITUTO FEDERAL DE EDUCACAO,CIENCIA E TECNOLOGIA DO CEARA
% GRADUACAO EM ENGENHARIA MECANICA: Gabriel Bandeira Holanda

clc; clear all; close all;

function GraficoTemp(matriz, deltaX, deltaY)

copia = matriz;

for i = 1: size(matriz,1)
        position = size(matriz,1)+1-i;
        copia(i,:) = matriz(position,:);
end

for i = 1: size(matriz, 1)
    for j = 1: size(matriz, 2)
        if copia(i,j) == 0
            copia(i,j) = 73.5;
        end
    end
end

pcolor(copia)
colorbar

end

%% Parâmetros da Malha
linesize  = 10;
colsize   = 10;

%% Dados da Questão
L = 0.048; % [m]
w = 0.006; % [m]
k = 50;    % [W/m*k]
h = 500;   % [W/m^2*k]
T_fluid = 303.15; % [k]
T_s     = 373.15; % [k]

%% Calcula o tamanho de cada espaço da malha
dx = L/colsize;
dy = w/linesize;

%% Matriz das temperaturas 
T = zeros(linesize,colsize);

%% Calcula e Gera as Equações da Malha usando abordagem de Gauss Seidel
maxerror  = 1;
error = 0.001;
while ( error < maxerror )
  TANTERIOR = T;
  for j=1:linesize 
   for i=1:colsize
    if (i==1 && j==1)
      A_ip1_j = k*(dy/dx);
      A_im1_j = k*(dy/(dx/2));
      A_i_jp1 = k*(dx/dy);
      A_i_jm1 = k*(dx/(dy/2));
      A_i_j   = A_ip1_j + A_im1_j + A_i_jp1 + A_i_jm1;
      % T(i-1,j) -> T_s AND T(i,j-1) -> T_P_up
      T_P_up = (k*T(j,i) + h*0.5*dy*T_fluid)/(k + h*0.5*dy);
      T(j,i)  = (A_ip1_j*T(j,i+1) + A_im1_j*T_s + A_i_jp1*T(j+1,i) + A_i_jm1*T_P_up)/A_i_j;
    elseif ((j==1) && (i~=1 && i~=colsize))
      A_ip1_j = k*(dy/dx);
      A_im1_j = k*(dy/dx);
      A_i_jp1 = k*(dx/dy);
      A_i_jm1 = k*(dx/(dy/2));
      A_i_j   = A_ip1_j + A_im1_j + A_i_jp1 + A_i_jm1;
      % T(i-1,j) -> T_s AND T(i,j-1) -> T_P_up
      T_P_up = (k*T(j,i) + h*0.5*dy*T_fluid)/(k + h*0.5*dy);
      T(j,i)  = (A_ip1_j*T(j,i+1) + A_im1_j*T(j,i-1) + A_i_jp1*T(j+1,i) + A_i_jm1*T_P_up)/A_i_j;
    elseif (j==1 && i==colsize)
      A_ip1_j = 0;          % É zero, pois não existe fluxo
      A_im1_j = k*(dy/dx);
      A_i_jp1 = k*(dx/dy);
      A_i_jm1 = k*(dx/(dy/2));
      A_i_j   = A_ip1_j + A_im1_j + A_i_jp1 + A_i_jm1;
      % A_ip1_j*T(i+1,j) -> 0
      T_P_up = (k*T(j,i) + h*0.5*dy*T_fluid)/(k + h*0.5*dy);
      T(j,i) = (A_im1_j*T(j,i-1) + A_i_jp1*T(j+1,i) + A_i_jm1*T_P_up)/A_i_j;
    elseif ((j~=1 && j~=linesize) && (i==1))
      A_ip1_j = k*(dy/dx);
      A_im1_j = k*(dy/(dx/2));
      A_i_jp1 = k*(dx/dy);
      A_i_jm1 = k*(dx/dy);
      A_i_j   = A_ip1_j + A_im1_j + A_i_jp1 + A_i_jm1;
      T(j,i)  = (A_ip1_j*T(j,i+1) + A_im1_j*T_s + A_i_jp1*T(j+1,i) + A_i_jm1*T(j-1,i))/A_i_j;
    elseif ((j~=1 && j~=linesize) && (i~=1 && i~=colsize))
      A_ip1_j = k*(dy/dx);
      A_im1_j = k*(dy/dx);
      A_i_jp1 = k*(dx/dy);
      A_i_jm1 = k*(dx/dy);
      A_i_j   = A_ip1_j + A_im1_j + A_i_jp1 + A_i_jm1;
      T(j,i)  = (A_ip1_j*T(j,i+1) + A_im1_j*T(j,i-1) + A_i_jp1*T(j+1,i) + A_i_jm1*T(j-1,i))/A_i_j;
    elseif ((j~=1 && j~=linesize) && (i==colsize))
      A_ip1_j = 0; % É zero, pois não existe fluxo
      A_im1_j = k*(dy/dx);
      A_i_jp1 = k*(dx/dy);
      A_i_jm1 = k*(dx/dy);
      A_i_j   = A_ip1_j + A_im1_j + A_i_jp1 + A_i_jm1;
      T(j,i)  = ( A_im1_j*T(j,i-1) + A_i_jp1*T(j+1,i) + A_i_jm1*T(j-1,i))/A_i_j;
    elseif (j==linesize && i==1)
      A_ip1_j = k*(dy/dx);
      A_im1_j = k*(dy/(dx/2));
      A_i_jp1 = k*(dx/(dy/2)); 
      A_i_jm1 = k*(dx/dy);
      A_i_j   = A_ip1_j + A_im1_j + A_i_jp1 + A_i_jm1;
      T_P_low = (k*T(j,i) + h*0.5*dy*T_fluid)/(k + h*0.5*dy);
      T(j,i)  = (A_ip1_j*T(j,i+1) + A_im1_j*T_s + A_i_jp1*T_P_low + A_i_jm1*T(j-1,i))/A_i_j;
    elseif (j==linesize) && (i~=1 && i~=colsize)
      A_ip1_j = k*(dy/dx);
      A_im1_j = k*(dy/dx);
      A_i_jp1 = k*(dx/(dy/2)); 
      A_i_jm1 = k*(dx/dy);
      A_i_j   = A_ip1_j + A_im1_j + A_i_jp1 + A_i_jm1;
      T_P_low = (k*T(j,i) + h*0.5*dy*T_fluid)/(k + h*0.5*dy);
      T(j,i)  = (A_ip1_j*T(j,i+1) + A_im1_j*T(j,i-1) + A_i_jp1*T_P_low + A_i_jm1*T(j-1,i))/A_i_j;
    elseif (j==linesize && i==colsize)
      A_ip1_j = 0;          % É zero, pois não existe fluxo
      A_im1_j = k*(dy/dx);
      A_i_jp1 = k*(dx/(dy/2)); 
      A_i_jm1 = k*(dx/dy);
      A_i_j   = A_ip1_j + A_im1_j + A_i_jp1 + A_i_jm1;
      T_P_low = (k*T(j,i) + h*0.5*dy*T_fluid)/(k + h*0.5*dy);
      T(j,i)  = (A_im1_j*T(j,i-1) + A_i_jp1*T_P_low + A_i_jm1*T(j-1,i))/A_i_j;
    end
    
   end
  end 
  TDIF = T - TANTERIOR;
  maxerror = max (max (TDIF));
  maxerror
endwhile

%% Converte as temperaturas encontradas de kelvin para Celsius:
T = T - 273.15;

%% Gera distribuição da temperatura na aleta:
GraficoTemp(T,dx,dy);

%% Gera distribuição da temperatura no linha mais central da malha 
%% para verificar a temperatura em função de x quando a linha for a central:
%Y_SHOW(1) = T_s - 273.15;
%X_SHOW(1) = 0;
%for i=1:colsize
% Y_SHOW(i+1) = T(round(linesize/2),i);
% X_SHOW(i+1) = i*dx*1000;
%end 
%plot(X_SHOW,Y_SHOW,'k-');
%ylim([30, 100]); % deixar temperatura de 30-100 visivel nos graficos






