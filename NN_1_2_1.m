%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Este programa ejecuta perfectamente 31 datos de entrada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

%load('IndiceBursatilRegressionBP.mat')
load('P001'); load('T001'); P=P001; T=T001;
format 
% P=-2:0.1:1.5; 
% P=P';
% T=1+sin(pi*P/4);
%T=pi*P/4+0.5*rand;
Q = size(P);

n1=2;
n2=1;
W1 = [-0.27; -0.41];
b1 = [-0.48; -0.13];
W2 = [0.09 -0.17];
b2 = [0.48];

alfa = 0.1;
for Epocas = 1:1000
    sum = 0;
     for q = 1:Q
         clc
%         %q = randi(Q);
        % Propagación de la entrada hacia la salida
        a1 = sigmoid(W1*P(q) + b1);
        a2 = W2*a1 + b2;
        % Retropropagación 
        e = T(q)-a2;
        s2 = -2*1*e;
        s1 = diag((1-a1).*a1)*W2'*s2;
%       % Actualización de pesos 
        
        W2 = W2 - alfa*s2*a1'
        b2 = b2 - alfa*s2;
        W1 = W1 - alfa*s1*P(q)'
        b1 = b1 - alfa*s1;
        % Sumando el error cuadratico 
        sum = e^2 + sum;
        %pause(0.5)
     end
    % Error cuadratico medio
    %emedio(Epocas) = sum/Q;
end
%figure, plot(emedio)
 
%Verificación grafica de la respuesta de la multicapa
%p = 5:1:50;
%p = -2:0.1:4;
%%
p = 0:0.1:1;
for q = 1:length(p)
        a1 = sigmoid(W1*p(q) + b1);
        a2(q) = W2*a1 + b2;
        %a3(q) = W3*a2 + b3;
end
plot(p,a2,P,T,'r.')

  

