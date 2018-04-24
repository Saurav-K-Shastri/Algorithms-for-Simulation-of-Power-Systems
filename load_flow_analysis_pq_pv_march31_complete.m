%% Load Flow Analysis with system having PV and PQ buses
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

    
    
    %%