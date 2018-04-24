function [d wr t] = RK4(Pm,Pmax1,Pmax2,Pmax3,H,w_s,delta_t,t_clear,t_fault,t_start,t_final,d_1,w_1)

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




    
    


