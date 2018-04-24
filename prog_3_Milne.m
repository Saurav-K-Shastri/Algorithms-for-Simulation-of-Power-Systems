%% Program3_Milne
% Predictor Corrector

clc;clear all;close all;

d = [0.4856 0.4866 0.4894 0.4942];%input('Enter delta_0 in rad: ');
wr = [0 0.1907 0.3810 0.5704];%input('Enter w_0 in rad/sec: ');
P =  [1.8 0.65 1.33 0.8];%input('Enter Pmax1,Pmax2,Pmax3,Pm: ');
Pmax1 = P(1);
Pmax2 = P(2);
Pmax3 = P(3);
Pm = P(4) ;

f_s = 60;%input('Enter frequency in Hz: ');
H = 5;%input('Enter Machine Constant in MW/mva: ');
t_start = 0;%input('Enter t_start in sec: ');
t_fault = 0;%input('Enter t_fault in sec: ');
t_clear = 0.38;%input('Enter t_clear in sec: ');
t_final = 5;%input('Enter t_final in sec: ');
delta_t = 0.01;% input('Enter delta_t in sec; ');

w_s = f_s*pi*2;
M = w_s/(2*H);

N = ceil((t_final-t_start)/delta_t);

t(1) = t_start;
for o=2:1:4
    t(o) = t(o-1) + delta_t;
end

count = 4;
for n = 1:1:N
    if t(n)<t_fault
        Pmax = Pmax1;
    elseif t(n)<t_clear
        Pmax = Pmax2;
    else 
        Pmax = Pmax3;
    end
    
    Dd = zeros(1,5);
    Dw = zeros(1,5);
    
    for p = 1:1:4
        Dd(p) = wr(count+p-4);
        Dw(p) = M*(Pm - Pmax*sin(d(count+p-4)));
    end
    
    d(count+1) = d(count-3) + (4/3)*(2*Dd(2) - Dd(3) +2*Dd(4))*delta_t;
    wr(count+1) = wr(count-3) + (4/3)*(2*Dw(2) - Dw(3) +2*Dw(4))*delta_t;
    
    Dd(5) = wr(count+1);
    Dw(5) = M*(Pm - Pmax*sin(d(count+1)));
    
    %Corrector
    d(count+1) = d(count-1) + (1/3)*(Dd(3) + 4*Dd(4) + Dd(5))*delta_t;
    wr(count+1) = wr(count-1) + (1/3)*(Dw(3) + 4*Dw(4) + Dw(5))*delta_t;
       
    t(count+1) = t(count) + delta_t;
    count = count+1;
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
