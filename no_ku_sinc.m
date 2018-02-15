% MÁQUINA DE INDUCCIÓN TRIFÁSICA ROTATIVA
% Antitransformación de Ku referencia SINCRONISMO
function [ist,irt,iso,isf,isb,par_mec,par_res] = no_ku_sinc(t,y)
%Declaracion de variables globales
global Lr M mdi Bw0 Bw1 Bw2 ws Vmax Lso Lro
global t_phi w_phi wm phi_Vs a120 a
%Referencia de Sincronismo
t_phi = ws*t;
w_phi = ws;
%Transformación inversa de las intensidades 
%Referencia KU SINCRONISMO
iso=y(:,1);
isf=real(y(:,3));
isb=imag(y(:,3));
ist = zeros(length(y),3);
irt = zeros(length(y),3);
%Variables del estator
ist(:,1)=(y(:,1)+exp(i*t_phi).*y(:,3)+exp(-i*t_phi). *conj(y(:,3)))/sqrt(3);
ist(:,2)=(y(:,1)+a^2*exp(i*t_phi).*y(:,3)+a*exp(-i*t_phi). *conj(y(:,3)))/sqrt(3);
ist(:,3)=(y(:,1)+a*exp(i*t_phi).*y(:,3)+a^2*exp(-i*t_phi).*conj(y(:,3)))/sqrt(3);
%Variables del rotor
irt(:,1)=(y(:,2)+exp(i*(t_phi-p*y(:,6))).*y(:,4)+exp(-i*(t_phi-p*y(:,6))). *conj(y(:,4)))/sqrt(3);
irt(:,2)=(y(:,2)+a^2*exp(i*(t_phi-p*y(:,6))).*y(:,4)+a*exp(-i*(t_phi-p*y(:,6))).*conj(y(:,4)))/sqrt(3);
irt(:,3)=(y(:,2)+a*exp(i*(t_phi-p*y(:,6))).*y(:,4)+a^2*exp(-i*(t_phi-p*y(:,6))).*conj(y(:,4)))/sqrt(3);
%Cálculo del par
par_mec=2*p*M*imag(y(:,3).*conj(y(:,4)));
par_res=(Bw0+p_resis)+Bw1*y(:,5)+Bw2*y(:,5).^2;