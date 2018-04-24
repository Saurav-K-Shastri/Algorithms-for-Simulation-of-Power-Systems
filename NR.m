clc;clear all;close all;

data =[1 2 0.12 0.16 0];% [1 2 0.028 0.096 0];% input('Enter fbus tbus R X hlc: ');
%%
[Ybus,nbus] = Bus_Adm(data);
Ybus
%%
V = [1+0j; 1+0j];%input('Enter V:');
Psp = -1 ;%-1.5;%input('Enter Psp');
Qsp = -0.5 ;%-0.5;%input('Enter Qsp');
S1 = zeros(nbus - 1,1); S = zeros(nbus - 1,1);
input1 = [0 ;1];
eps = 0.001;
delta = 0;
for n=1:10%iterations
    
    S1 = zeros(nbus - 1,1); S = zeros(nbus - 1,1);
    for k1=2:nbus
        for k2=1:nbus
            S1(k1-1) = S1(k1-1) + Ybus(k1,k2)*V(k2);
        end
        S(k1-1) = conj(V(k1))*S1(k1-1);
    end
    
    Pcal = real(S);
    Qcal = -imag(S);
    delp = Psp-Pcal;
    delq = Qsp-Qcal;
    delpq = [delp; delq];
    
    if abs(delp)<eps && abs(delq)<eps
        break
    else      
    for k1 = 2:nbus
        for k2 = 2:nbus
            if k2==k1
                PQdii = (-j)*conj(V(k1))*S1(k1-1)+conj(V(k1))*Ybus(k1,k2)*V(k2)*j;
                PQVii = conj(V(k1))*S1(k1-1)+conj(V(k1))*Ybus(k1,k2)*V(k1);
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
    J = [ J1 J2; J3 J4];
    
    corr = inv(J)*(delpq);
    corr(2) = corr(2)*abs(V(2));
    delta = delta + corr(1);
    V(2) = corr(2) + abs(V(2));
    V(2) = V(2)*(cos(delta)+1j*sin(delta));
    end    
    
end

S1 = zeros(nbus,1); S = zeros(nbus,1);

    for k1=1:nbus
        for k2=1:nbus
            S1(k1) = S1(k1) + Ybus(k1,k2)*V(k2);
        end
        S(k1) = conj(conj(V(k1))*S1(k1));
    end
   

S
Sloss = sum((S))
V