%PROGRAMA MOTOR%
%Referencia Estator%
t_phi=0;
w_phi=0;
%Referencia Sincronismo%
t_phi=ws*t;
w_phi=ws;
%Referencia Rotor%
t_phi=p*y(6);
w_phi=p*y(5);
%Tensiones del estator
Vsa=Vmax*cso(ws*t+phi_Vs);
Vsb=Vmax*cso(ws*t+phi_Vs-a120);
Vsc=Vmax*cso(ws*t+phi_Vs+a120);
%Constantes
a=exp(i*2*pi/3);
a120=2*pi/3;
%Tensiones transformadas
Vso=(Vsa+Vsb+Vsc)/sqrt(3);
Vsf =(Vsa+a*Vsb+a^2*Vsc)*exp(-i*t_phi)/sqrt(3);
Vro = 0; Vrf = 0;
% Ecuaciones diferenciales
%iso=y(1); iro=y(2); isf=y(3); irf=y(4); w=y(5); ang_mec=y(6);
dy(1,1)=1/Lso*(Vso-Rs*y(1));
dy(2,1)=1/Lro*(Vro-Rr*y(2));
dy(3:4,1)=-[Lr -M;-M Ls]/(Lr*Ls-M^2)*([Rs+i*w_phi*Ls  i*w_phi*M ; i*(w_phi-wm)*M  Rr+i*(w_phi-wm)*Lr]*[ y(3) ; y(4) ] - [Vsf ; Vrf]);
dy(5,1) = 1/mdi * (par_mec - par_res);
dy(6,1) = y(5);
% Par motor y par resistente    
par_mec = 2 * p * M * imag(y(3) * conj(y(4)));    
par_res = Bw0 + Bw1 * y(5) + Bw2 * y(5)^2;
% Vector de tiempo
t=[to:inc:tf];





% Calculo del sistema de ecuaciones ordinarias (funcion ode45)
[t,y] = ode45('ku_sinc',[to:inc:tf],[0 0 is ir y5 ang_meco]);



plot(t,real(ist(:,1)),'r',t,real(ist(:,2)),'b',t,real(ist(:,3)),'k')
title('isa, isb, isc')
plot(t,real(irt(:,1)),'r',t,real(irt(:,2)),'b',t,real(irt(:,3)))
title('ira, irb, irc')

plot(t,y(:,1),'r',t,real(y(:,3)),'b',t,imag(y(:,3)),'k')
title('iso, isf, isb')

plot(t,y(:,5)*60/2/pi,'r')
title('Velocidad')
ylabel('[r.p.m.]')

plot(t,par_mec,'r',t,par_res,'b')
title('Par mecánico (rojo) - Par resistente (azul)')
ylabel('[Nm]')







