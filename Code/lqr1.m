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

% Checking for the controllability of the given system
if (rank(ctrb(A,B))==6)
disp("Rank of ctrb matches order of A, system is controllable")
else
disp("Rank of ctrb doesnt matche order of A, system is uncontrollable")
end

global Q
Q=[100 0 0 0 0 0;
0 100 0 0 0 0;
0 0 30000 0 0 0;
0 0 0 30000 0 0;
0 0 0 0 30000 0;
0 0 0 0 0 30000];

global R
R=1;

global C
C = eye(6);% To form a 6 X 6 identity matrix

global D
D = 0; % Initialising the D matrix to be Zero

global K_val
disp("Now, seeing the results using an LQR controller")
[K_val, P_mat, Poles] = lqr(A,B,Q,R);%In-built MATLAB code
Poles

y0 = [0; 0; 30; 0; 60; 0];
t_int = 0:0.001:1000;%defining the timespan

[t1,y1] = ode45(@pendlinear,t_int,y0); %Linearization with initial conditions
[t2,y2] = ode45(@pendnonlinear,t_int,y0); %Non-linear systems
figure
plot(t1,y1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','x-dot(t)', 'theta_1(t)', 'theta-dot_1(t)', 'theta_2(t)', 'theta-dot_2(t)')
title('Response of a linear system')
grid on

figure
plot(t2,y2)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','x-dot(t)', 'theta_1(t)', 'theta-dot_1(t)', 'theta_2(t)', 'theta-dot_2(t)')
title('Response of a nonlinear system')
grid on

function dydt = pendlinear(t,y)
global A B K_val
u = -K_val * y;
dydt = A*y + B*u;
end

function dydt = pendnonlinear(t,y)
global K_val g m1 m2 l1 l2 M
F =-K_val*y;
dydt=zeros(6,1);
dydt(1) = y(2);
dydt(2)=(F-(g/2)*(m1*sind(y(3))+m2*sind(2*y(5)))-(m1*l1*(y(4)^2)*sind(y(3)))-(m2*l2*(y(6)^2)*sind(y(5))))/(M+m1*((sind(y(3)))^2)+m2*((sind(y(5)))^2));%xDD
dydt(3)= y(4);%theta 1D;
dydt(4)= (dydt(2)*cosd(y(3))-g*(sind(y(3))))/l1';%theta 1 Ddot;
dydt(5)= y(6);%theta 2D
dydt(6)= (dydt(2)*cosd(y(5))-g*(sind(y(5))))/l2;%theta 2Ddot;
end


