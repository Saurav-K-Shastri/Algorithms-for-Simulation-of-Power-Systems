%% Program1_Modefied Euler Method
clc;clear all;close all;

d_1 = 0.44854;%input('Enter delta_0 in rad: ');
w_1 = 0;%input('Enter w_0 in rad/sec: ');
P =  [1.714 0.63 1.33 0.8];%input('Enter Pmax1,Pmax2,Pmax3,Pm: ');
Pmax1 = P(1);
Pmax2 = P(2);
Pmax3 = P(3);
Pm = P(4) ;

f_s = 60;%input('Enter frequency in Hz: ');
H = 5;%input('Enter Machine Constant in MW/mva: ');
t_start = 0;%input('Enter t_start in sec: ');
t_fault = 0;%input('Enter t_fault in sec: ');
t_clear = 0.3;%input('Enter t_clear in sec: ');
t_final = 5;%input('Enter t_final in sec: ');
delta_t = 0.01;% input('Enter delta_t in sec; ');

w_s = f_s*pi*2;
M = w_s/(2*H);

N = ceil(t_final-t_start)/delta_t;

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
    
    d_slope1(n) = wr(n);
    wr_slope1(n) = M*(Pm - Pmax*sin(d(n)));
    
    dp = d(n) + d_slope1(n)*delta_t;
    wrp = wr(n) + wr_slope1(n)*delta_t;
    
    d_slope2(n) = wrp;
    wr_slope2(n) = M*(Pm - Pmax*sin(dp));    
    
    d(n+1) = d(n) + 0.5*(d_slope1(n)+d_slope2(n))*delta_t;
    wr(n+1) = wr(n) + 0.5*(wr_slope1(n)+wr_slope2(n))*delta_t; 
    
    t(n+1) = t(n) + delta_t;
end

figure(1)
plot(t,d);
grid on
title('delta vs time');
xlabel('time');
ylabel('delta');


figure(2)
plot(t,wr);
grid on
title('omega vs time');
xlabel('time');
ylabel('omega');



    
    


