%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERADOR DE DATOS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Este programa genera datos normalizados para el entrenamiento de la red 
% neuronal:
% DESCRIPCI�N: Variables
%   inpunt_inf_v_sensor: limite inferior de entrada en voltios 
%   inpunt_sup_v_sensor: limite superior de entrada en voltios 
%   step_v: paso entre una y la siguiente entrada (con criterio)
%     
%   ejemplo: 
%   Rango de voltaje emitito por el sensor [0.005 : 0.05]V
%   inpunt_inf_v_sensor: 0.005
%   inpunt_sup_v_sensor: 0.05
%   step_v=0.001 
%
%   Autor: waposat
%   Lima, Abril 04 de 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;
%sensor temperatura: [0.005 : 0.05]V --> [10�C : 80�C]
mv_s1_C=(0.005:0.0001:0.01);
cant_dat=length(mv_s1_C);

for i=1:cant_dat
    s1_C(1,i)=(1555.56*mv_s1_C(1,i)+2.222)+2*sin(0.05*i)+0.1*rand();
    %s1_C(1,i)=sin(mv_s1_C(1,i)*100);
    
end
max_mv=max(mv_s1_C);
min_mv=min(mv_s1_C);
max_C=max(s1_C);
min_C=min(s1_C);

% normalizaci�n de datos
P001=((mv_s1_C'-min_mv)*(1))/(max_mv-min_mv)
T001=((s1_C'-min_C)*(1))/(max_C-min_C)
plot(P001,T001)

% guardar datos para ser usado desde otros archivos sin perderse aunque se
% se utilize "clear all".
save P001;  % ventor de entrada para la red neuronal
save T001;  % vector de objetivos (target), para hallar los errores
save max_mv;    
save min_mv;    
save max_C;
save min_C;
