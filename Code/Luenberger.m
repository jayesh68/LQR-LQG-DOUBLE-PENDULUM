%clearing all the previous outputs
clc
clear all
close all

% Given information
global M m1 m2 L1 L2
M=1000;%Mass of the cart
m1=100;%mass of Pendulum 1
m2=100;%mass of Pendulum 2
L1=20;%length of the string of Pendulum 1
L2=10;%length of the string of Pendulum 2

global g
g=9.81; %declaring the value of the accelertaion due to gravity in m/

global A
A=[0 1 0 0 0 0;
0 0 -(m1*g)/M 0 -(m2*g)/M 0;
0 0 0 1 0 0;
0 0 -((M+m1)*g)/(M*L1) 0 -(m2*g)/(M*L1) 0;
0 0 0 0 0 1;
0 0 -(m1*g)/(M*L2) 0 -(g*(M+m2))/(M*L2) 0];

global B
B=[0; 1/M; 0; 1/(M*L1); 0; 1/(M*L2)];

% Checking for the controllability of the given system
if (rank(ctrb(A,B))==6)
disp("Rank of ctrb matches order of A, system is controllable")
else
disp("Rank of ctrb doesnt matche order of A, system is uncontrollable")
end

global C
C = eye(6);% To form a 6 X 6 identity matrix
global D
D = 0; % Initialising the D matrix to be Zero
y0 = [5; 0; 30; 0; 60; 0];
t_int = 0:0.001:1000;%defining the timespan

Bd = 0.1*eye(6);
Vn = 0.01;

global c1
c1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
global c3
c3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
global c4
c4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];

tspan = 0:0.1:500;
q0 = [2 0 deg2rad(30) 0 deg2rad(60) 0];

Bd = 0.1*eye(6); %Process Noise
Vn = 0.01; %Measurement Noise
[Lue1,~,~] = lqe(A,Bd,c1,Bd,Vn*eye(3));
[Lue3,~,~] = lqe(A,Bd,c3,Bd,Vn*eye(3));
[Lue4,~,~] = lqe(A,Bd,c4,Bd,Vn*eye(3));
Ac1 = A-(Lue1*c1);
Ac3 = A-(Lue3*c3);
Ac4 = A-(Lue4*c4);
e_sys1 = ss(Ac1,[B Lue1],c1,0);
e_sys3 = ss(Ac3,[B Lue3],c3,0);
e_sys4 = ss(Ac4,[B Lue4],c4,0);

unitStep = 0*tspan;
unitStep(200:length(tspan)) = 1;
d = [1;0;0];
sys1 = ss(A,B,c1,d);
sys3 = ss(A,B,c3,d);
sys4 = ss(A,B,c4,d);
[y1,t] = lsim(sys1,unitStep,tspan);
[x1,t] = lsim(e_sys1,[unitStep;y1'],tspan);
[y3,t] = lsim(sys3,unitStep,tspan);
[x3,t] = lsim(e_sys3,[unitStep;y3'],tspan);
[y4,t] = lsim(sys4,unitStep,tspan);
[x4,t] = lsim(e_sys4,[unitStep;y4'],tspan);


figure();
hold on
plot(t,y1(:,1),'r','Linewidth',2)
plot(t,x1(:,1),'k--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','Estimated x(t)')
title('Response for output vector at step input: x(t)')
hold off

figure();
hold on
plot(t,y3(:,1),'r','Linewidth',2)
plot(t,y3(:,3),'b','Linewidth',2)
plot(t,x3(:,1),'k--','Linewidth',1)
plot(t,x3(:,3),'m--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','theta_2(t)','Estimated x(t)','Estimated theta_2(t)')
title('Response for output vector at step input: (x(t),theta_2(t))')
hold off

figure();
hold on
plot(t,y4(:,1),'r','Linewidth',2)
plot(t,y4(:,2),'g','Linewidth',2)
plot(t,y4(:,3),'b','Linewidth',2)
plot(t,x4(:,1),'k--','Linewidth',1)
plot(t,x4(:,2),'r--','Linewidth',1)
plot(t,x4(:,3),'m--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','theta_1(t)','theta_2(t)','Estimated x(t)','Estimated theta_1(t)','Estimated theta_2(t)')
title('Response for output vector at step input: (x(t),theta_1(t),theta_2(t))')
hold off

[t,q1] = ode45(@(t,q)linearObs1(t,q,Lue1),tspan,q0);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer for output vector: x(t)')
legend('x')
hold off

[t,q3] = ode45(@(t,q)linearObs3(t,q,Lue3),tspan,q0);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer for output vector: (x(t),theta_2(t))')
legend('x','theta_2')
hold off


[t,q4] = ode45(@(t,q)linearObs4(t,q,Lue4),tspan,q0);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer for output vector: (x(t),theta_1(t),theta_2(t))')
legend('x','theta_1','theta_2')
hold off

[t,q1] = ode45(@(t,q)nonLinearObs1(t,q,1,Lue1),tspan,q0);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: x(t)')
legend('x')
hold off


[t,q3] = ode45(@(t,q)nonLinearObs3(t,q,1,Lue3),tspan,q0);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: (x(t),theta_2(t))')
legend('x','theta_2')
hold off


[t,q4] = ode45(@(t,q)nonLinearObs4(t,q,1,Lue4),tspan,q0);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: (x(t),theta_1(t),theta_2(t))')
legend('x','theta_1','theta_2')
hold off

function dQe = linearObs4(t,Qe,Lue4)
global A B c4
y4 = [Qe(1); Qe(3); Qe(5)];
K = 1; % feedback = 1;
dQe = (A+B*K)*Qe + Lue4*(y4 - c4*Qe);
end

function dQe = linearObs1(t,Qe,Lue1)
global A B c1
y1 = [Qe(1); 0; 0];
K = 1; % feedback = 1;
dQe = (A+B*K)*Qe + Lue1*(y1 - c1*Qe);
end

function dQe = linearObs3(t,Qe,Lue3)
global A B c3
y3 = [Qe(1); 0; Qe(5)];
K = 1; % feedback = 1;
dQe = (A+B*K)*Qe + Lue3*(y3 - c3*Qe);
end

function dQ = nonLinear(t,y,F)
global M m1 m2 L1 L2 g
x = y(1);
dx = y(2);
t1 = y(3);
dt1 = y(4);
t2 = y(5);
dt2 = y(6);
dQ=zeros(6,1);
dQ(1) = dx;
dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)));
dQ(3) = dt1;
dQ(4) = (cos(t1)*dQ(2)-g*sin(t1))/L1;
dQ(5) = dt2;
dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2;
end

function dQ = nonLinearObs1(t,y,F,Lue1)
global M m1 m2 L1 L2 g
x = y(1);
dx = y(2);
t1 = y(3);
dt1 = y(4);
t2 = y(5);
dt2 = y(6);
dQ=zeros(6,1);
y1 = [x; 0; 0];
c1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
sum = Lue1*(y1-c1*y);
dQ(1) = dx + sum(1);
dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
dQ(3) = dt1+sum(3);
dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
dQ(5) = dt2 + sum(5);
dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end

function dQ = nonLinearObs3(t,y,F,Lue3)
global M m1 m2 L1 L2 g c3
x = y(1);
dx = y(2);
t1 = y(3);
dt1 = y(4);
t2 = y(5);
dt2 = y(6);
dQ=zeros(6,1);
y3 = [x; 0; t2];
sum = Lue3*(y3-c3*y);
dQ(1) = dx + sum(1);
dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
dQ(3) = dt1+sum(3);
dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
dQ(5) = dt2 + sum(5);
dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end

function dQ = nonLinearObs4(t,y,F,Lue4)
global M m1 m2 L1 L2 g c4
x = y(1);
dx = y(2);
t1 = y(3);
dt1 = y(4);
t2 = y(5);
dt2 = y(6);
dQ=zeros(6,1);
y4 = [x; t1; t2];
sum = Lue4*(y4-c4*y);
dQ(1) = dx + sum(1);
dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
dQ(3) = dt1+sum(3);
dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
dQ(5) = dt2 + sum(5);
dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end