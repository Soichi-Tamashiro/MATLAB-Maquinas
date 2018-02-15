% KU_SINC
%Funcion que realiza la transformacion KU de los parametros del motor en la referencia de Sincronismo    
function dy=ku_sinc(t,y)
% Variables Globales    
global nmotor In Mn Vrms frec p Rs Ls Rr Lr M mdi Bw0 Bw1 Bw2 ws Vmax Lso Lro    
global t_phi w_phi wm phi_Vs a120 a p_resis
%Referencia Sincronismo %
t_phi = ws*t;
w_phi = ws;
wm = p*y(5);
% Tensiones de estator
Vsa = Vmax*cos(ws*t+phi_Vs);
Vsb = Vmax*cos(ws*t+phi_Vs-a120);
Vsc = Vmax*cos(ws*t+phi_Vs+a120);
% Tensiones transformadas
Vso = (Vsa+Vsb+Vsc)/sqrt(3);
Vsf = (Vsa+a*Vsb+a^2*Vsc)*exp(-i*t_phi)/sqrt(3);
Vro = 0; Vrf = 0;
% Par motor y par resistente
par_mec=2*p*M*imag(y(3)*conj(y(4)));
par_res =(Bw0+p_resis)+Bw1*y(5)+Bw2*y(5)^2;
%Ecuaciones diferenciales
%iso=y(1);iro=y(2);isf=y(3);irf=y(4);w=y(5); ang_mec=y(6);
dy(1,1) = 1/Lso*(Vso-Rs*y(1));
dy(2,1) = 1/Lro*(Vro - Rr*y(2));
dy(3:4,1)= - [Lr -M; -M   Ls]/(Lr*Ls-M^2)*([Rs+i*w_phi*Ls    i*w_phi*M  i*(w_phi-wm)*M    Rr+i*(w_phi-wm)*Lr]*[ y(3) ; y(4) ]-[Vsf ; Vrf]);
dy(5,1)=1/mdi*(par_mec-par_res);
dy(6,1) = y(5);