## @Gabriel Bandeira Holanda
## Simula��o de um Pist�o Alternativo 15/08/2019
clc; clear all; close all;

#================================================#
#===== Dados Adquiridos em Artigos e Livros =====#
#================================================#
Dis_eixos = 0;          # DesalinhArea_Mento Entre Eixos [m]
L_Manivela = 0.07;      # Comprimento Manivela  [m]
L_Biela  = 0.24;        # Comprimento Biela     [m]
D_Piston = 0.1;         # DiArea_Metro do Pist�o    [m]
L_Piston = 0.08;        # Comprimento do Pist�o [m]
RPM = 2000;             # Rotacao por Minuto    [RPM]
T0 = 315.15;            # Temperatura Inicial   [k] [42 C]
Coef_Pol = 1.4; # Coeficiente Politropico (Ar)
P_Desc = 820000;        # Press�o de Descarga [Pa]
P_Suc  = 200000;        # Press�o de Suc��o   [Pa]

#================================================#
#==================== C�lculos ==================#
#================================================#
R_Piston = D_Piston/2;      # Raio do Pist�o  [m]
Cur_Piston = L_Manivela*2;  # Curso do Pist�o [m]
Cl = (pi*(R_Piston^2))*Cur_Piston;  # Cilindradas [m^(3)]
Vol_M  = Cl*0.3;              # Volume Morto (30% da Cilindradas)  [m^(3)] 
Area_M = (2*Vol_M)/R_Piston;  # �rea com Ar devido ao Volume Morto [m^2]
Cpms = L_Biela - L_Manivela;  # Dist�ncia Entre o Ponto Morto Superior 
                              #(PMS) e Eixo da Manivela [m]
w = (2*pi*RPM)/60;            # Velocidade Angular

n = 1000;            # Numero de Vezes em que se Divide uma Rota��o
Gap_T = (60/RPM)/n;  # Gap do Tempo Realizar uma Rota��o Dividido 'n' Partes
# Fun��o do Tempo [t(i)] e o Angulo em Fun��o do tempo [teta(t)]
for i=1:n
  t(i)    = i*Gap_T;
  teta(i) = w*t(i);
end

for i=1:n 
  # Y_t (Posi��o Pist�o em Fun��o do Tempo)
  Y_t(i) = abs(Cpms + L_Manivela*cos(w*t(i)) - sqrt( (L_Biela^2) - (L_Manivela*sin(w*t(i))-Dis_eixos)^2 ));
  # V_t (Volume em Fun��o do Tempo)
  V_t(i) = pi*(R_Piston^2)*Y_t(i) + Vol_M;
  # A_t (Area em Fun��o do Tempo)
  A_t(i) = 2*pi*R_Piston*Y_t(i) + Area_M;
end

# Velocidade M�dia do Pist�o em Fun��o do Tempo
for i=1:(n-1)
Vm_Piston(i) = (Y_t(i)-Y_t(i+1))/Gap_T;
end
Vm_Piston(n) = (Y_t(n)-Y_t(1  ))/Gap_T;

P_t(1) = P_Desc;
for i=2:n
  P_t(i) = P_t(i-1)/((V_t(i)/V_t(i-1))^Coef_Pol);
  if  P_t(i) <= P_Suc
    # considerar o Fluxo de Massa
    P_t(i) = P_Suc;
  end
  if  P_t(i) >= P_Desc
    # considerar o Fluxo de Massa
    P_t(i) = P_Desc;
  end
end

T_t(1) = T0;
for i=2:n
  T_t(i) = T_t(i-1)/(P_t(i-1)/P_t(i))^( (Coef_Pol - 1)/Coef_Pol );
end

figure (1)
plot(teta*180/pi, Vm_Piston);
title ('Velocidade M�dia Pist�o vs. Angulo');
xlabel ('Angulo [ � ]' ); ylabel ('Velocidade M�dia Pist�o [ m^3 ]' );

figure (2)
plot(teta*180/pi, T_t);
title ('Volume vs. Angulo');
xlabel ('Angulo [ � ]' ); ylabel ('Volume [ m^3 ]' );



figure (3)
plot(V_t,P_t/1000);
title ('Press�o vs. Volume');
xlabel (' Volume [ m^3 ]' ); ylabel (' Press�o [ KPa ]' );
maxvalue = max(V_t) + (max(V_t)*0.2);
axis([0, maxvalue, 0, 1000]);
 