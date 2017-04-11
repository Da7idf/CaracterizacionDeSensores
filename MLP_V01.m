%%%%%%%%%%%%%%%%%%%%%%%%% RED NEURONAL - TEMPERATURA %%%%%%%%%%%%%%%%%%%%%%
% Este programa n
% DESCRIPCIÓN: Variables
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
clear all, close all, clc

%load('IndiceBursatilRegressionBP.mat')
load('P001'); load('T001'); 
P=P001; 
T=T001;
plot(P001,T001)
Q = size(P,1);
%plot(P001,T001)
%%
n1 = 5;
n2 = 30;
cf=1;
W1 = cf*rand(n1,1);
b1 = cf*rand(n1,1);
W2 = cf*rand(n2,n1);
b2 = cf*rand(n2,1);
W3 = cf*rand(1,n2);
b3 = cf*rand;


alfa = 0.1;
for Epocas = 1:1000
    error_epoca = 0;
    for q = 1:Q
        clc
        q = randi(Q);
        % Propagación de la entrada hacia la salida
        a1 = sigmoid(W1*P(q) + b1);
        a2 = sigmoid(W2*a1 + b2);
        a3 = W3*a2 + b3;
        % Retropropagación de la sensibilidades
        e = T(q)-a3;
        s3 = -2*1*e;
        s2 = diag((1-a2).*a2)*W3'*s3;
        s1 = diag((1-a1).*a1)*W2'*s2;
        % Actualización de pesos sinapticos y polarizaciones
        W3 = W3 - alfa*s3*a2'
        b3 = b3 - alfa*s3;
        W2 = W2 - alfa*s2*a1'
        b2 = b2 - alfa*s2;
        W1 = W1 - alfa*s1*P(q)'
        b1 = b1 - alfa*s1;  
        % Sumando el error cuadratico 
        error_epoca = e^2 + error_epoca;
    end
    % Error cuadratico medio
    emedio(Epocas) = error_epoca/Q;

end
subplot(2,1,1)
plot(emedio)

p = 0:0.01:1;
for q = 1:length(p)
        a1 = sigmoid(W1*p(q) + b1);
        a2 = sigmoid(W2*a1 + b2);
        a3(q) = W3*a2 + b3;
end

p_new=p*(max_mv-min_mv)+min_mv;
a3_new=a3*(max_C-min_C)+min_C;
P_new=P*(max_mv-min_mv)+min_mv;
T_new=T*(max_C-min_C)+min_C;

subplot(2,1,2)
plot(p_new,a3_new,P_new,T_new,'r')
%plot(p,a3,P,T,'r')


