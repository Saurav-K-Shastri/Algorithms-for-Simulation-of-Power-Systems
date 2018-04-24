clc;clear all;close all;
data1 =[1 2 0 0.1 0.002];
data2 = [ 1.02 1 -1.5 -0.8];
x = 0.2;

[Ybus,nbus] = Bus_Adm(data1);

[E delta V2] = load_flow_analysis(Ybus,nbus,data2,x);

Pm = E*abs(V2)*sin(delta)/(0.2+0.1);

Pmax1 = E*abs(V2)/(x+data1(4));

Pmax2 = Pmax1/6;

Pmax3 = Pmax1*4/3;

H = 5;
w_s = 314;
delta_t = 0.001;
t_clear = 0.496; % 0.496 is critical clearing time
t_fault = 0;
t_start = 0;
t_final = 20;
d_1 = delta;
w_1 = 0;
[d wr t] = RK4(Pm,Pmax1,Pmax2,Pmax3,H,w_s,delta_t,t_clear,t_fault,t_start,t_final,d_1,w_1);

figure(1)
subplot(1,2,1);
plot(t,d);
grid on
title('delta vs time');
xlabel('time');
ylabel('delta');

subplot(1,2,2);
plot(t,wr);
grid on
title('omega vs time');
xlabel('time');
ylabel('omega');



