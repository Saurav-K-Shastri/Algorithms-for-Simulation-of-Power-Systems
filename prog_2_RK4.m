%% Program2_RK4

clc;clear all;close all;

d_1 = 0.4605;%input('Enter delta_0 in rad: ');
w_1 = 0;%input('Enter w_0 in rad/sec: ');
P =  [1.8 0.65 1.33 0.8];%input('Enter Pmax1,Pmax2,Pmax3,Pm: ');
Pmax1 = P(1);
Pmax2 = P(2);
Pmax3 = P(3);
Pm = P(4) ;

f_s = 60;%input('Enter frequency in Hz: ');
H = 5;%input('Enter Machine Constant in MW/mva: ');
t_start = 0;%input('Enter t_start in sec: ');
t_fault = 0;%input('Enter t_fault in sec: ');
t_clear = 0.9;%input('Enter t_clear in sec: ');
t_final = 5;%input('Enter t_final in sec: ');
delta_t = 0.01;% input('Enter delta_t in sec; ');

w_s = f_s*pi*2;
M = w_s/(2*H);

N = ceil((t_final-t_start)/delta_t);

d(1) = d_1;
wr(1) = w_1;

t(1) = t_start;

for n = 1:1:N
    if t(n)<t_fault
        Pmax = Pmax1;
    elseif t(n)<t_clear
        Pmax = Pmax2;
    else 
        Pmax = Pmax3;
    end
    
    k1 = wr(n)*delta_t;
    l1 = M*(Pm - Pmax*sin(d(n)))*delta_t;
    
    k2 = (wr(n)+(l1/2))*delta_t;
    l2 = M*(Pm - Pmax*sin(d(n)+(k1/2)))*delta_t;
    
    k3 = (wr(n)+(l2/2))*delta_t;
    l3 = M*(Pm - Pmax*sin(d(n)+(k2/2)))*delta_t;
    
    k4 = (wr(n)+(l3))*delta_t;
    l4 = M*(Pm - Pmax*sin(d(n)+(k3)))*delta_t;
    
    d(n+1) = d(n) + (1/6)*(k1+(2*k2)+(2*k3)+k4);
    wr(n+1) = wr(n) + (1/6)*(l1+(2*l2)+(2*l3)+l4);
    
    t(n+1) = t(n) + delta_t;
end

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



    
    


