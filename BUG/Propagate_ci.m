function   [Rho,Ni]=Propagate_ci(dt,t_end,Rho_0,Ni_0)

global W

time_span=0:dt:t_end;

[~,N_time]=size(time_span);
Rho=zeros(N_time,3);
Ni=zeros(N_time,3);

Rho_Ni_0=[Rho_0(1);Rho_0(2);Rho_0(3);Ni_0(1);Ni_0(2);Ni_0(3)];

for i=1:1:N_time
tau=time_span(i);

Ct=cos(W*tau);
St=sin(W*tau);



PHI=[4-3*Ct       0      0       St/W          2/W*(1-Ct)          0;
    6*(St-W*tau)  1      0    -2/W*(1-Ct)    1/W*(4*St-3*W*tau)    0;
        0         0      Ct       0                0           1/W*St;
     3*W*St       0      0        Ct             2*St              0;
   -6*W*(1-Ct)    0      0      -2*St           4*Ct-3             0;
        0         0    -W*St      0                0              Ct ];
    
    
    PHI=[4-3*Ct       0      0       St/W          2/W*(1-Ct)          0;
    6*(St-W*tau)  1      0    -2/W*(1-Ct)    1/W*(4*St-3*W*tau)    0;
        0         0      Ct       0                0                0;
     3*W*St       0      0        Ct             2*St              0;
   -6*W*(1-Ct)    0      0      -2*St           4*Ct-3             0;
        0         0      0     0                0              0 ];
    
    
    
Rho_Ni=PHI*Rho_Ni_0;
    
Rho(i,:)=Rho_Ni(1:3)';
Ni(i,:)=Rho_Ni(4:6)';

end


return