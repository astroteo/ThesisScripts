function   [Ni_0_plus,DV]=Compute_ci(tau,Rho_0,Ni_0_minus,X_goal)

global W


Ct=cos(W*tau);
St=sin(W*tau);



det=(4*St)/(W^3) -(8*Ct*St)/(W^3) +(4*Ct^2*St)/(W^3)+(4*St^3)/(W^3) -(3*St^2*tau)/(W^2);

N_tau_inv=1/det*[(4*St^2)/(W^2)-(3*St*tau)/W,     -((2*St)/(W^2))+(2*Ct*St)/(W^2),                        0;
                 (2*St)/(W^2)-(2*Ct*St)/(W^2),              St^2/(W^2),                                   0;
                           0,                                  0,              4/(W^2)-(8*Ct)/(W^2)+(4*Ct^2)/(W^2)+(4*St^2)/(W^2)-(3*St*tau)/W];
                   
M_tau=[-3*Ct+4,        0,   0;
        6*St-6*W*tau,  1,   0;                   
               0    ,  0,  Ct];


Ni_0_plus=N_tau_inv*(X_goal-M_tau*Rho_0);
DV=norm(Ni_0_plus-Ni_0_minus);

return