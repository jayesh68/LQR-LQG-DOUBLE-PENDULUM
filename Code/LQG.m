%clearing all the previous outputs
clc
clear all
close all

% Given information
global M m1 m2 l1 l2 g
M=1000;%Mass of the cart
m1=100;%mass of Pendulum 1
m2=100;%mass of Pendulum 2
l1=20;%length of the string of Pendulum 1
l2=10;%length of the string of Pendulum 2
g=9.81; %declaring the value of the accelertaion due to gravity in m/

global A
A=[0 1 0 0 0 0;
0 0 -(m1*g)/M 0 -(m2*g)/M 0;
0 0 0 1 0 0;
0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
0 0 0 0 0 1;
0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];

global B
B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

global C
C = [1 0 0 0 0 0];

global D
D = 0; % Initialising the D matrix to be Zero

%Initial Condition vector
y0 = [5; 0; 30; 0; 60; 0];

t_int = 0:0.001:1000;%defining the timespan

Bd = 0.1*eye(6);
Vn = 0.01;

global Q
Q=[100 0 0 0 0 0;
0 0 0 0 0 0;
0 0 0 0 0 0;
0 0 0 0 0 0;
0 0 0 0 0 0;
0 0 0 0 0 0];

global R
R=1;

vd = 0.1*eye(6);
vn = 0.01;

global K_val
[K_val, P_mat, Poles] = lqr(A, B, Q, R);
Ac1 = A-(K_val'*C);
e_sys1 = ss(Ac1,[B K_val'],C, 0);

global kalman_gain
kalman_gain = lqe(A, vd, C, vd, vn);
%kalman_gain = kalman(e_sys1, 0.1, 0.1, 0.1)

[t2,y2] = ode45(@pendnonlinear,t_int,y0);

figure();
hold on
plot(t_int,y2(:,1),'r','Linewidth',2)
ylabel('state variables')
xlabel('time (sec)')
title('Feedback controller using smallest output vector x(t)')
legend('x(t)')
grid on
hold off

function dydt = pendnonlinear(t,y)

global K_val g m1 m2 l1 l2 M kalman_gain C
F =-K_val*y;
sum = kalman_gain * (y(1) - C * y(1));
dydt=zeros(6,1);
dydt(1) = y(2) + sum(1);
dydt(2)=(F-(g/2)*(m1*sind(y(3))+m2*sind(2*y(5)))-(m1*l1*(y(4)^2)*sind(y(3)))-(m2*l2*(y(6)^2)*sind(y(5))))/(M+m1*((sind(y(3)))^2)+m2*((sind(y(5)))^2));%xDD
dydt(3)= y(4) + sum(2);
dydt(4)= (dydt(2)*cosd(y(3))-g*(sind(y(3))))/l1';
dydt(5)= y(6) + sum(3);
dydt(6)= (dydt(2)*cosd(y(5))-g*(sind(y(5))))/l2;
end