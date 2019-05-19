clear all
close all



%orbit
Mu=(3.9856e+14);
r_Earth=6380*10^3;
h_orb=400*10^3;

T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)

t_end_CW=1/100*T_orb
d_TC=200
omega=sqrt(Mu/((r_Earth+h_orb)^3));
w=omega;
setGlobalW(w);


%target
R_0t=[(r_Earth+h_orb);0;0];
V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];
R_T=r_Earth+h_orb;
V_T=norm(V_0t);



%actuator & cubesat.
m_cubesat=5;
f_sat=2/5*1e-3*m_cubesat;% ref1: a_sat=4*1e-6N && m_cubesat=3Kg 
                         % ref2[PLASMA]: f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3 http://pepl.engin.umich.edu/thrusters/CAT.html
                         % ref3[COMMERCIAL]: http://www.busek.com/cubesatprop__main.htm
                         % ref4[NANOSAT]:f_sat=100*1e-6 N &&  Imp_duration= 2*1e-3 s && Isp= 50:100 s http://www.cubesatshop.com/index.php?page=shop.product_details&flypage=flypage.tpl&product_id=74&vmcchk=1&option=com_virtuemart&Itemid=65f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3
F_sat=f_sat
M_cubesat=m_cubesat

setGlobalf_sat(F_sat)
setGlobalm_cubesat(M_cubesat)

%LQR controller.
a_sat=f_sat/m_cubesat
Wu=[1/a_sat , 0,    0;
    0,      1/a_sat, 0;
    0, 0,       1/a_sat];

Wz=eye(6);

wz=0 %<-need a weight to reduce control effort on {x,y} while can stay max on z  

rr_p=1
Cz=1/2*sqrt(1)*[0    0     0    0      0    0;
                0 (1+wz/2) 0    0      0    0;
                0    0     0    0      0    0;
                0    0     0  (1+wz/2) 0    0;
                0    0     0    0      0    0;
                0    0     0    0      0  (1-wz)]-1/2*sqrt(Mu/(norm(R_T)^2))*[(1+wz/2) 0     0      0    0   0;
                                                                              0        0     0      0    0   0;
                                                                              0        0   (1+wz/2) 0    0   0;
                                                                              0        0      0     0    0   0;
                                                                              0        0      0     0 (1-wz) 0;
                                                                              0        0      0     0    0   0];
Dz=zeros(6);




A=[0     1    0     0   0   0;
  3*w^2  0    0    2*w  0   0;
   0     0    0     1   0   0;
   0   -2*w   0     0   0   0;
   0     0    0     0   0   1;
   0     0    0     0 -w^2  0];
   
B=[0 0 0;
   1 0 0; 
   0 0 0;
   0 1 0; 
   0 0 0; 
   0 0 1];


 rr_u=0.01;
 rr_z=2;
   
Q = Cz'*(rr_z*Wz)*Cz;


R = rr_u*Wu;


K = lqr(A,B,Q,R);

K_lqr=-K

setGlobalk_lqr(K_lqr);




kv=-0.05*1000
kp=-7.5*1e-5*1000

Kv=kv;
Kp=kp;

setGlobalKp(Kp)
setGlobalKv(Kv)


%% standard (DOCKING) test


% X_0t=[R_0t;V_0t];
% X_0c=[R_0c;V_0c];
% X_0t_c=[X_0t;X_0c];
% 
% options = odeset('RelTol',1e-13,'AbsTol', 1e-45);
% [t_LQR,X_LQR]= ode113('r2BP_t_c_LQR_sat',t_end,X_0t_c,options);
% 
% R_t_LQR=[X_LQR(:,1),X_LQR(:,2),X_LQR(:,3)];
% R_c_LQR=[X_LQR(:,7),X_LQR(:,8),X_LQR(:,9)];
% V_t_LQR=[X_LQR(:,4),X_LQR(:,5),X_LQR(:,6)];
% V_c_LQR=[X_LQR(:,10),X_LQR(:,11),X_LQR(:,12)];
% 
% Rrel_I_LQR=R_c_LQR-R_t_LQR;
% Vrel_I_LQR=V_c_LQR-V_t_LQR;
% 
% [T,~]=size(t_LQR);
% 
% Rrel_LVLH_LQR=zeros(T,3);
% Vrel_LVLH_LQR=zeros(T,3);
% 
% for t=1:T
%     time=t_LQR(t);
%     theta=omega*time;
%     ROT=[ cos(theta), sin(theta),0;
%            -sin(theta), cos(theta),0;
%                  0       ,    0       ,1];
%              
%         Rrel_LVLH_LQR(t,:)=(ROT*Rrel_I_LQR(t,:)')';
%         Vrel_LVLH_LQR(t,:)=(ROT*(Vrel_I_LQR(t,:)'-cross([0;0;omega],Rrel_I_LQR(t,:)')))';
% end
% 
% 
% options = odeset('RelTol',1e-13,'AbsTol', 1e-7);
% [t_LQR_SAT,X_LQR_SAT]= ode113('r2BP_t_c_LQR',t_end,X_0t_c,options);
% 
% R_t_LQR_SAT=[X_LQR_SAT(:,1),X_LQR_SAT(:,2),X_LQR_SAT(:,3)];
% R_c_LQR_SAT=[X_LQR_SAT(:,7),X_LQR_SAT(:,8),X_LQR_SAT(:,9)];
% V_t_LQR_SAT=[X_LQR_SAT(:,4),X_LQR_SAT(:,5),X_LQR_SAT(:,6)];
% V_c_LQR_SAT=[X_LQR_SAT(:,10),X_LQR_SAT(:,11),X_LQR_SAT(:,12)];
% 
% Rrel_I_LQR_SAT=R_c_LQR_SAT-R_t_LQR_SAT;
% Vrel_I_LQR_SAT=V_c_LQR_SAT-V_t_LQR_SAT;
% 
% [S,~]=size(t_LQR_SAT);
% 
% Rrel_LVLH_LQR_SAT=zeros(S,3);
% Vrel_LVLH_LQR_SAT=zeros(S,3);
% 
% for s=1:S
%     t=t_LQR_SAT(s);
%     theta=omega*t;
%     ROT=[ cos(theta), sin(theta),0;
%            -sin(theta), cos(theta),0;
%                  0       ,    0       ,1];
%              
%         Rrel_LVLH_LQR_SAT(s,:)=(ROT*Rrel_I_LQR_SAT(s,:)')';
%         Vrel_LVLH_LQR_SAT(s,:)=(ROT*(Vrel_I_LQR_SAT(s,:)'-cross([0;0;omega],Rrel_I_LQR_SAT(s,:)')))';
% end
% figure()
% plot3(Rrel_LVLH_LQR(:,1),Rrel_LVLH_LQR(:,2),Rrel_LVLH_LQR(:,3),'--b',Rrel_LVLH_LQR(1,1),Rrel_LVLH_LQR(1,2),Rrel_LVLH_LQR(1,3),'og',Rrel_LVLH_LQR(end,1),Rrel_LVLH_LQR(end,2),Rrel_LVLH_LQR(end,3),'or')
% title('controlloed LQR ,NO-SAT(--b) in LVLH ref ')
% hold on
% grid on 
% 
% figure()
% plot3(Rrel_LVLH_LQR_SAT(:,1),Rrel_LVLH_LQR_SAT(:,2),Rrel_LVLH_LQR_SAT(:,3),'--b',Rrel_LVLH_LQR_SAT(1,1),Rrel_LVLH_LQR_SAT(1,2),Rrel_LVLH_LQR_SAT(1,3),'og',Rrel_LVLH_LQR_SAT(end,1),Rrel_LVLH_LQR_SAT(end,2),Rrel_LVLH_LQR_SAT(end,3),'or')
% title('controlloed LQR, SAT(--b) in LVLH ref ')
% hold on
% grid on 
% 
% 
% 
% disp('Rrel_LVLH_LQR_DOCK, END:')
% Rrel_LVLH_LQR(end,:)
% 
% disp('Vrel_LVLH_LQR_DOCK, END:')
% Vrel_LVLH_LQR(end,:)
% 
% disp('Rrel_LVLH_LQR_DOCK, SAT:')
% Rrel_LVLH_LQR_SAT(end,:)
% 
% disp('Vrel_LVLH_LQR_DOCK, SAT:')
% Vrel_LVLH_LQR_SAT(end,:)
%% x_obs with null velocity w.r.t the target.
               %SP--> {Rho-x_obs=0; Ni-cross(omega_t,(x_obs-Rho))};
               %SK--> {Rho-x_obs=0; Ni=0} 
               %CTR--> {SK: for  -l_box < Rho(j) < l_box [with brake] SP: for translation}

x_obs=20
y_obs=70
z_obs=10

L_box=2;

X_obs=[x_obs;y_obs;z_obs];%-> observation point given in LVLH r.f.

setGlobalx_obs(X_obs)
setGloball_box(L_box)


%chaser

% R_0c=(R_0t+X_obs-[rand(1,1)*0.1*1e2;rand(1,1)*0.1*1e2;0]); %-> pert. only in plane: circular orbit
% % V_0c=(V_0t+[rand(1,1)*0.3*1e-1;rand(1,1)*0.1*1e-1;0]);
% V_0c=V_0t;

%  R_0c=R_0t+X_obs; %--> maintaince (stationary keeping) | valid in this form since LVLH && I coincident @ t=0.|
%  V_0c=V_0t+zeros(3,1);
%  
R_0c=(R_0t+X_obs+rand(3,1)*0.25*1e2); 
V_0c=(V_0t-rand(3,1)*0.25*1e-2);



t_end=1/0.5*T_orb
 X_0t=[R_0t;V_0t];
  X_0c=[R_0c;V_0c];
   X_0t_c=[X_0t;X_0c];

Flag_DV_SK=0;
setGlobalflag_DV_SK(Flag_DV_SK)

options = odeset('RelTol',1e-13,'AbsTol', 1e-6);
[t_LQR_FT,X_LQR_FT]= ode113('r2BP_t_c_LQR_SP',t_end,X_0t_c,options);



R_t_LQR_FT=[X_LQR_FT(:,1),X_LQR_FT(:,2),X_LQR_FT(:,3)];
R_c_LQR_FT=[X_LQR_FT(:,7),X_LQR_FT(:,8),X_LQR_FT(:,9)];
V_t_LQR_FT=[X_LQR_FT(:,4),X_LQR_FT(:,5),X_LQR_FT(:,6)];
V_c_LQR_FT=[X_LQR_FT(:,10),X_LQR_FT(:,11),X_LQR_FT(:,12)];

Rrel_I_LQR_FT=R_c_LQR_FT-R_t_LQR_FT;
Vrel_I_LQR_FT=V_c_LQR_FT-V_t_LQR_FT;

[F,~]=size(t_LQR_FT);

Rrel_LVLH_LQR_FT=zeros(F,3);
Vrel_LVLH_LQR_FT=zeros(F,3);

for f=1:F
    t=t_LQR_FT(f);
    theta=-omega*t;
    ROT=[ cos(theta), sin(theta),0;
           -sin(theta), cos(theta),0;
                 0       ,    0       ,1];
             
        Rrel_LVLH_LQR_FT(f,:)=(ROT*Rrel_I_LQR_FT(f,:)')';
        Vrel_LVLH_LQR_FT(f,:)=(ROT*(Vrel_I_LQR_FT(f,:)'-cross([0;0;omega],Rrel_I_LQR_FT(f,:)')))';
end


%% Free Drift
      % reaching X_obs through C-W equations for different tf

tau_span=500:10:5*3600;
[~,N_tau]=size(tau_span)
DV_tot=zeros(N_tau,1);
for i=1:N_tau
    tau=tau_span(i);
    
Ct=cos(w*tau);
St=sin(w*tau);

det=(4*St)/(w^3) -(8*Ct*St)/(w^3) +(4*Ct^2*St)/(w^3)+(4*St^3)/(w^3) -(3*St^2*tau)/(w^2);
N_tau_inv=1/det*[(4*St^2)/(w^2)-(3*St*tau)/w,     -((2*St)/(w^2))+(2*Ct*St)/(w^2),                        0;
                 (2*St)/(w^2)-(2*Ct*St)/(w^2),              St^2/(w^2),                                   0;
                           0,                                  0,              4/(w^2)-(8*Ct)/(w^2)+(4*Ct^2)/(w^2)+(4*St^2)/(w^2)-(3*St*tau)/w];
                   
M_tau=[-3*Ct+4,        0,   0;
        6*St-6*w*tau,  1,   0;                   
               0    ,  0,  Ct];
 S_tau=[3*w*St,     0 , 0;
        6*w*Ct-6*w, 0 , 0 ;   
            0     ,  0, St/w];
         
 T_tau=[Ct   , 2*St   , 0;
        -2*St, 4*Ct-3 , 0;
          0  ,  0     , Ct];
                 
 Ni_X0_plus=N_tau_inv*(X_obs-M_tau*(R_0c-R_0t)); %<-- ROT_0*(R_0c-R_0c) in general !! assumed LVLH & I coincident @ t=0 !!
 Ni_X0_minus=V_0c-V_0t-cross([0;0;omega],(R_0c-R_0t));
 
 Ni_X_obs_minus=S_tau*(R_0c-R_0t)+T_tau*Ni_X0_plus;
 Ni_X_obs_plus=zeros(3,1);
 
 DV_tot(i)= norm(Ni_X0_plus-Ni_X0_minus)+norm((Ni_X_obs_plus-Ni_X_obs_minus));
 
end

[~,i_min]=min(DV_tot);
t_end_FD=tau_span(i_min)
tau=t_end_FD;

Ct=cos(w*tau);
St=sin(w*tau);

N_tau_inv=1/det*[(4*St^2)/(w^2)-(3*St*tau)/w,     -((2*St)/(w^2))+(2*Ct*St)/(w^2),                        0;
                 (2*St)/(w^2)-(2*Ct*St)/(w^2),              St^2/(w^2),                                   0;
                           0,                                  0,              4/(w^2)-(8*Ct)/(w^2)+(4*Ct^2)/(w^2)+(4*St^2)/(w^2)-(3*St*tau)/w];

M_tau=[-3*Ct+4,        0,   0;
        6*St-6*w*tau,  1,   0;                   
               0    ,  0,  Ct];
           
Ni_X0_plus=N_tau_inv*(X_obs-M_tau*(R_0c-R_0t))

X_0c_FD=[R_0c;V_0t+Ni_X0_plus+cross([0;0;omega],(R_0c-R_0t))];
   X_0t_c_FD=[X_0t;X_0c];

options = odeset('RelTol',1e-13,'AbsTol', 1e-30);
[t_FD,X_FD]= ode113('r2BP_t_c',t_end_FD,X_0t_c_FD,options);
[FD,~]=size(t_FD);



R_t_FD=[X_LQR_FT(:,1),X_LQR_FT(:,2),X_LQR_FT(:,3)];
R_c_FD=[X_LQR_FT(:,7),X_LQR_FT(:,8),X_LQR_FT(:,9)];
V_t_FD=[X_LQR_FT(:,4),X_LQR_FT(:,5),X_LQR_FT(:,6)];
V_c_FD=[X_LQR_FT(:,10),X_LQR_FT(:,11),X_LQR_FT(:,12)];

Rrel_I_FD=R_c_FD-R_t_FD;
Vrel_I_FD=V_c_FD-V_t_FD;

Rrel_LVLH_FD=zeros(FD,3);
Vrel_LVLH_FD=zeros(FD,3);

for fd=1:FD
    t=t_FD(fd);
    theta=-omega*t;
    ROT=[ cos(theta), sin(theta),0;
           -sin(theta), cos(theta),0;
                 0       ,    0       ,1];
             
        Rrel_LVLH_FD(fd,:)=(ROT*Rrel_I_FD(fd,:)')';
        Vrel_LVLH_FD(fd,:)=(ROT*(Vrel_I_FD(fd,:)'-cross([0;0;omega],Rrel_I_FD(fd,:)')))';
end


t_end_CW=t_end_FD
x_0t_c_FD=[R_0c-R_0t;Ni_X0_plus];
options = odeset('RelTol',1e-13,'AbsTol', 1e-30);
[t_FD,x_FD]= ode113('Chol_Wilt_Hill',t_end_CW,x_0t_c_FD,options);


Rho_FD=[x_FD(:,1),x_FD(:,2),x_FD(:,3)];
%% Plots
disp('Rrel_LVLH_LQR_FT, END:')
Rrel_LVLH_LQR_FT(end,:)

disp('Vrel_LVLH_LQR_FT, END:')
Vrel_LVLH_LQR_FT(end,:)

figure()
plot3(Rrel_LVLH_LQR_FT(:,1),Rrel_LVLH_LQR_FT(:,2),Rrel_LVLH_LQR_FT(:,3),'--b',Rrel_LVLH_LQR_FT(1,1),Rrel_LVLH_LQR_FT(1,2),Rrel_LVLH_LQR_FT(1,3),'og',Rrel_LVLH_LQR_FT(end,1),Rrel_LVLH_LQR_FT(end,2),Rrel_LVLH_LQR_FT(end,3),'or',x_obs,y_obs,z_obs,'*g')
title('controlloed LQR, FT(--b) in LVLH ref ')
hold on
grid on

figure()
plot(t_LQR_FT,Rrel_LVLH_LQR_FT(:,1))
title('time <-->x (LVLH)')
hold on
grid on

figure()
plot(t_LQR_FT,Rrel_LVLH_LQR_FT(:,2))
title('time <-->y (LVLH)')
hold on
grid on

figure()
plot(t_LQR_FT,Rrel_LVLH_LQR_FT(:,3))
title('time <-->z (LVLH)')
hold on
grid on


figure()
plot(t_LQR_FT,Rrel_I_LQR_FT(:,1))
title('time <-->x (IN)')
hold on
grid on

figure()
plot(t_LQR_FT,Rrel_I_LQR_FT(:,2))
title('time <-->y (IN)')
hold on
grid on

figure()
plot(t_LQR_FT,Rrel_I_LQR_FT(:,3))
title('time <-->z (IN)')
hold on
grid on

figure()
plot3(Rrel_LVLH_FD(:,1),Rrel_LVLH_FD(:,2),Rrel_LVLH_FD(:,3),'--b',Rrel_LVLH_FD(1,1),Rrel_LVLH_FD(1,2),Rrel_LVLH_FD(1,3),'og',Rrel_LVLH_FD(end,1),Rrel_LVLH_FD(end,2),Rrel_LVLH_FD(end,3),'or',x_obs,y_obs,z_obs,'*g',Rho_FD(:,1),Rho_FD(:,2),Rho_FD(:,3),'--k')
title('Free Drift(--b) in LVLH ref ')
hold on
grid on
