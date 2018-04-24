%% Transient Analysis with system having PV and PQ buses
clc;clear all;close all;

%% inputs
data1 = [1 3 0 1/40 0 ;  1 2 0 1/20 0 ; 2 3 0 1/20 0];%input('fbus tbus R X HLC: ');
[Ybus,nbus] = Bus_Adm(data1);

V = [1 1 1.05];%input('Enter voltages either given or assumed (in sequence): ');
Psp = [-5 4];%input('Psp from 2 onwards: ');
Qsp = -4;%input('Qsp from 2 onwards: ');
pq = length(Qsp);
pv = nbus  - pq -1;

DELTA = zeros(1,nbus-1);
deld = zeros(1,nbus-1);
delv = ones(1,pq);

del = [deld delv];

eps = 0.001;

del_corrected = DELTA;
%% Calculations

for n=1:100
    
    I = zeros(1,nbus - 1); S = zeros(1,nbus - 1);
    
    for k1=2:nbus
        for k2=1:nbus
            I(k1-1) = I(k1-1) + Ybus(k1,k2)*V(k2);
        end
        S(k1-1) = conj(V(k1))*I(k1-1);
    end
    
    Pcal = real(S);
    Qcal = -imag(S);
    Qcal = Qcal(1:pq);
    
    delp = Psp-Pcal;
    delq = Qsp-Qcal;
    
    delpq = [delp delq];
    
   
    for k1 = 2:nbus
        for k2 = 2:nbus
            if k2==k1
                PQdii = (-j)*conj(V(k1))*I(k1-1)+conj(V(k1))*Ybus(k1,k2)*V(k2)*j;
                PQVii = conj(V(k1))*I(k1-1)+conj(V(k1))*Ybus(k1,k2)*V(k1);
                J1(k1-1,k1-1)= real(PQdii);
                J3(k1-1,k1-1) = -imag(PQdii);
                J2(k1-1,k1-1)=real(PQVii);
                J4(k1-1,k1-1)=-imag(PQVii);
            else
                PQdim = conj(V(k1))*Ybus(k1,k2)* V(k2)*j;
                PQVim = conj(V(k1))*Ybus(k1,k2)* V(k2);
                J1(k1-1,k2-1) = real(PQdim);
                J3(k1-1,k2-1) = -imag(PQdim);
                J2(k1-1,k2-1) = real(PQVim);
                J4(k1-1,k2-1) = -imag(PQVim);
            end
        end
    end
    J = [J1 J2(:,1:pq); J3(1:pq,:) J4(1:pq,1:pq)];
    
    del_cap = (inv(J)*(delpq.')).';
    
    del_corrected = del_corrected + del_cap(1:nbus-1);
    
    V_correction = del_cap(nbus:nbus+pq-1)*abs(V(nbus-1:nbus+pq-2));
    
    V_corrected = abs(V(2:nbus-1)) + V_correction;
    
    V(2:nbus-1) = V_corrected;
    
    del_corrected = [0 del_corrected];
    
    V = abs(V).*(cos(del_corrected)+1j*sin(del_corrected));   
    del_corrected = del_corrected(2:end);
       
end

xd1 = 1.2;
xd3 = 0.2;

S1 = zeros(nbus,1); S = zeros(nbus,1);

    for k1=1:nbus
        for k2=1:nbus
            S1(k1) = S1(k1) + Ybus(k1,k2)*V(k2);
        end
        S(k1) = conj(conj(V(k1))*S1(k1));
    end
    
    
E1 = V(1) + S1(1)*1j*xd1;
E3 = V(3) + S1(3)*1j*xd3;

Pm1 = abs(E1)*abs(V(1))*sin(angle(E1) - angle(V(1)))/xd1;
Pm3 = abs(E3)*abs(V(3))*sin(angle(E3) - angle(V(3)))/xd3;

delta_t = 0.0001;
Pm = zeros(2,1);
Pm(1) = Pm1;
Pm(2) = Pm3;

H1 = 7.5;
H2 = 5;
w_s = 314;

t_clear = 0.162;

t_fault = 0;
t_start = 0;
t_final = 5;
d_1 = zeros(2,1);
d_1(1) = angle(E1);
d_1(2) = angle(E3);

w_1 =zeros(2,1);


E = zeros(3,1);
E(1) = E1;
E(3) = E3;

[d wr t] = RK4_transient(V,S,S1,E,Ybus,Pm,H1,H2,w_s,delta_t,t_clear,t_fault,t_start,t_final,d_1,w_1,xd1,xd3);

figure(1)

subplot(221)
plot(t,d(1,:))

subplot(222)
plot(t,wr(1,:))

subplot(223)
plot(t,d(2,:))

subplot(224)
plot(t,wr(2,:))




