function [d wr t] = RK4_transient(V,S,I,E,Ybus,Pm,H1,H2,w_s,delta_t,t_clear,t_fault,t_start,t_final,d_1,w_1,xd1,xd3)
V = V.'; 
S2 = V(2)*conj(I(2));

Y2 = conj(S(2))/abs(V(2))^2;

M1 = w_s/(2*H1);
M2 = w_s/(2*H2);

N = ceil((t_final-t_start)/delta_t);

d(1:2,1) = d_1;
wr(1:2,1) = w_1;

t(1) = t_start;

Ybus_new = zeros(5,5);

Ybus_new(1:3,1:3) = Ybus;
Ybus_new(5,5) = 1/(j*xd3);
Ybus_new(4,4) = 1/(j*xd1);

Ybus_new(1,1) = Ybus(1,1)+(1/(j*xd1));
Ybus_new(3,3) = Ybus(3,3) + (1/(j*xd3));

Ybus_new(1,4) = -1/(j*xd1);
Ybus_new(4,1) = Ybus_new(1,4);

Ybus_new(3,5) = -1/(j*xd3);

Ybus_new(5,3) = Ybus_new(3,5);

Ybus_new(2,2) = Ybus(2,2) + Y2;

Ybus_new1 = Ybus_new;

Emag1 = abs(E(1));
Emag3 = abs(E(3));

for n = 1:1:N
    
    if t(n)<t_clear
        Ybus_new1 = Ybus_new;
        Ybus_new1(:,1) = 0;
        Ybus_new1(1,:) = 0;
        Y = Ybus_new1;
        V(1,n+1) = 0;
        V(2,n+1) = -(((Y(2,1)*V(1,n+1))/Y(2,2)) + ((Y(2,3)*V(3,n))/Y(2,2)));
        V(3,n+1) = -(((Y(3,1)*V(1,n+1))/Y(3,3)) + ((Y(3,2)*V(2,n+1))/Y(3,3)) + ((Y(3,5)*E(3,n))/Y(3,3)));
        
    else
        Ybus_new1 = Ybus_new;        
        Ybus_new1(1,3) = 0;
        Ybus_new1(3,1) = 0;
        Y = Ybus_new1;
        Y(1,1) =  Y(1,1) +j*40;
        Y(3,3) =  Y(3,3) +j*40;
        V(1,n+1) = -(((Y(1,2)*V(1,n))/Y(1,1)) + ((Y(1,3)*V(3,n))/Y(1,1)) + ((Y(1,4)*E(1,n))/Y(3,3)));
        V(2,n+1) = -(((Y(2,1)*V(1,n+1))/Y(2,2)) + ((Y(2,3)*V(3,n))/Y(2,2)));
        V(3,n+1) = -(((Y(3,1)*V(1,n+1))/Y(3,3)) + ((Y(3,2)*V(2,n+1))/Y(3,3)) + ((Y(3,5)*E(3,n))/Y(3,3)));  
    end
        
    Pmax1 = abs(E(1,n))*abs(V(1,n+1))/xd1;
    Pmax2 = abs(E(3,n))*abs(V(3,n+1))/xd3;
    
    k11 = wr(1,n)*delta_t;
    l11 = M1*(Pm(1) - Pmax1*sin(d(1,n)))*delta_t;

    k21 = (wr(1,n)+(l11/2))*delta_t;
    l21 = M1*(Pm(1) - Pmax1*sin(d(1,n)+(k11/2)))*delta_t;
    
    k31 = (wr(1,n)+(l21/2))*delta_t;
    l31 = M1*(Pm(1) - Pmax1*sin(d(1,n)+(k21/2)))*delta_t;
    
    k41 = (wr(1,n)+(l31))*delta_t;
    l41 = M1*(Pm(1) - Pmax1*sin(d(1,n)+(k31)))*delta_t;
     
    d(1,n+1) = d(1,n) + (1/6)*(k11+(2*k21)+(2*k31)+k41);
    wr(1,n+1) = wr(1,n) + (1/6)*(l11+(2*l21)+(2*l31)+l41);
    
    E(1,n+1) = Emag1*(cos(d(1,n+1)) + 1j*sin(d(1,n+1)));

    k12 = wr(2,n)*delta_t;
    l12 = M2*(Pm(2) - Pmax2*sin(d(2,n)))*delta_t;

    k22 = (wr(2,n)+(l12/2))*delta_t;
    l22 = M2*(Pm(2) - Pmax2*sin(d(2,n)+(k12/2)))*delta_t;
    
    k32 = (wr(2,n)+(l22/2))*delta_t;
    l32 = M2*(Pm(2) - Pmax2*sin(d(2,n)+(k22/2)))*delta_t;
    
    k42 = (wr(2,n)+(l32))*delta_t;
    l42 = M2*(Pm(2) - Pmax2*sin(d(2,n)+(k32)))*delta_t;
           
    d(2,n+1) = d(2,n) + (1/6)*(k12+(2*k22)+(2*k32)+k42);
    wr(2,n+1) = wr(2,n) + (1/6)*(l12+(2*l22)+(2*l32)+l42);
    
    E(3,n+1) = Emag3*(cos(d(2,n+1)) + 1j*sin(d(2,n+1)));
        
    t(n+1) = t(n) + delta_t;
end

