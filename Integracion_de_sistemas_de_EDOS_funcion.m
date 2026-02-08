function dx = Integracion_de_sistemas_de_EDOS_funcion(~,x)
Ip  = 0.2;
Mc  = 0.2;
lp  = 0.3;
Mp  = 0.1;
Fc  = 2;
Beq = 0.25;
g   = 9.81;
Bp  = 0.2;

xc     = x(1);
xcdot  = x(2);
alpha  = x(3);
adot   = x(4);

D = (Mc+Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(alpha)^2;

xdd = ( (Ip + Mp*lp^2)*Fc ...
      + Mp^2*lp^2*g*cos(alpha)*sin(alpha) ...
      - (Ip + Mp*lp^2)*Beq*xcdot ...
      - (Ip*Mp*lp + Mp^2*lp^3)*adot^2*sin(alpha) ...
      - Mp*lp*Bp*adot*cos(alpha) ) / D;

add = ( (Mc + Mp)*Mp*g*lp*sin(alpha) ...
      - (Mc + Mp)*Bp*adot ...
      + Fc*Mp*lp*cos(alpha) ...
      - Mp^2*lp^2*adot^2*sin(alpha)*cos(alpha) ...
      - Beq*Mp*lp*xcdot*cos(alpha) ) / D;

dx = zeros(4,1);
dx(1) = xcdot;
dx(2) = xdd;
dx(3) = adot;
dx(4) = add;


end
