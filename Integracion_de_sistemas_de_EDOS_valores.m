clear; clc; close all

alpha0 = 1*pi/180;
x0     = 0;
xdot0  = 0;
adot0  = 0;

X0 = [x0; xdot0; alpha0; adot0];

tspan = [0 10];

[t,X] = ode45(@(t,x) Integracion_de_sistemas_de_EDOS_funcion(t,x),tspan,X0);

figure
subplot(2,1,1)
plot(t,X(:,1),'LineWidth',1.5)
grid on
ylabel('x_c (m)')

subplot(2,1,2)
plot(t,X(:,3)*180/pi,'LineWidth',1.5)
grid on
ylabel('\alpha (deg)')
xlabel('Tiempo (s)')