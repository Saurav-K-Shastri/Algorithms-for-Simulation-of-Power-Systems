clc;clear all;close all;
data1 =[1 2 0 0.1 0.002];
data2 = [ 1.02 1 -1.5 -0.8];
x = 0.2;

[Ybus,nbus] = Bus_Adm(data1);

[E delta V2] = load_flow_analysis(Ybus,nbus,data2,x);

Pm = E*abs(V2)*sin(delta)/(0.2+0.1);

Pmax1 = E*abs(V2)/(x+data1(4));

Pmax2 = Pmax1/6;

Pmax3 = Pmax1*3/4;

H = 5;
w_s = 314;
delta_t = 0.001;
t_clear = 0.2; % 0.496 is critical clearing time
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

%% Area calculations

DD = 1:0.1:180;

P1 = Pmax1*sind(DD);
P2 = Pmax2*sind(DD);
P3 = Pmax3*sind(DD);
Pm1=Pm*ones(size(DD));

figure(2)

plot(DD,P1)
hold on
plot(DD,P2)
plot(DD,P3)
plot(DD,Pm1)

delta_0 = d(find(round(t,4) == t_fault)+1)*180/pi;
delta_clear = d(find(round(t,4) == t_clear)+1)*180/pi;
delta_final = asind(Pm/Pmax3);

delta_final = 180 - delta_final;

A1 = Pm*(delta_clear - delta_0) - Pmax2*(cosd(delta_clear) - cosd(delta_0))
A2 = Pmax3*(cosd(delta_final) - cosd(delta_clear)) - Pm*(delta_final - delta_clear)

DDD = zeros(size(DD));
DDD(:) = nan;

DDD(find(DD==round(delta_0,1))) = Pm;

DDD(find(DD==round(delta_clear,1))) = Pmax3*sind(delta_clear);

DDD(find(DD==round(delta_final,1))) = Pm;

DDD(find(DD==round(delta_clear,1))) = Pmax3*sind(delta_clear);
stem(DD,DDD)
