%% cholessy-wilthshare test IN HILL R.F.
close all
clear all
%orbit
Mu=(3.9856e+14);
r_Earth=6380*10^3;
h_orb=800*10^3;

T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)

t_end_CW=1/0.5*T_orb
d_TC=200
omega=sqrt(Mu/((r_Earth+h_orb)^3))

%target
R_0t=[(r_Earth+h_orb);0;0];
V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];
R_T=r_Earth+h_orb;
V_T=norm(V_0t);

%chaser
R_0c=(R_0t-[d_TC;0;0]); %-> pert. only in plane: circular orbit
R_C=R_T-d_TC
omega_c=sqrt(Mu/(norm(R_0c)^3));
V_0c=[0;sqrt(Mu/norm(R_0c));0];
V_C=norm(V_0c);



%Cholessy-Wiltshare integration
Rho_0=R_0c-R_0t;



w=omega;
setGlobalW(w);



Ni_0=(V_0c-V_0t)-cross([0;0;w],(R_0c-R_0t));

x_0=[Rho_0;Ni_0];
options = odeset('RelTol',1e-13,'AbsTol', 1e-40);
[t_cw,x]= ode113('Chol_Wilt_Hill',t_end_CW,x_0,options);

Rho=[x(:,1),x(:,2),x(:,3)];
 Ni=[x(:,4),x(:,5),x(:,6)];
 
%analytical Cholessy-Wiltshare
[N,~]=size(t_cw);
x_AN=zeros(N,3);
for i=1:N
    t=t_cw(i);
   
    x_AN(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*t)+Ni_0(1)/w*sin(w*t)+4*Rho_0(1)+2/w*Ni_0(2);
    x_AN(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*t)-2/w*Ni_0(1)*cos(w*t)+(6*w*Rho_0(1)+3*Ni_0(2))*t-Rho_0(2)+2/w*Ni_0(1));
    x_AN(i,3)=Rho_0(3)*cos(w*t)-Ni_0(3)*w*sin(w*t);
end

% target and chaser positions
R_t=zeros(N,3);
R_c=zeros(N,3);
V_t=zeros(N,3);
V_c=zeros(N,3);

for k=1:N
    tiempo=t_cw(k);
    theta=omega*tiempo;
    theta_c=omega_c*tiempo;
    
    R_t(k,:)=[R_T*cos(theta),R_T*sin(theta),0];
    R_c(k,:)=[R_C*cos(theta_c),R_C*sin(theta_c),0];
    V_t(k,:)=[V_T*(-sin(theta)),V_T*cos(theta),0];
    V_c(k,:)=[V_C*(-sin(theta_c)),V_C*cos(theta_c),0];
    
end



%relative positions & velocities

Rrel_I=R_c-R_t;
Vrel_I=V_c-V_t;

Rrel_LVLH=zeros(N,3);
Vrel_LVLH=zeros(N,3);
Rrel_BASTA=zeros(N,3);
Rrel_I_rec=zeros(N,3);
Vrel_I_rec=zeros(N,3);

%LVLH-->IN
for l=1:N
    t=t_cw(l);
    theta=omega*t;
    ROT=[cos(theta), sin(theta),0;
         -sin(theta), cos(theta),0;
             0       ,    0     ,1];
        
        Rrel_I_rec(l,:)=(ROT'*Rho(l,:)')';
        Vrel_I_rec(l,:)=(ROT'*(Ni(l,:)'+cross([0;0;omega],Rho(l,:)')))';
        Rrel_BASTA(l,:)=(ROT*Rrel_I_rec(l,:)')';
        
end

%IN-->LVLH
for j=1:N
    t=t_cw(j);
    theta=omega*t;
    ROT=[ cos(theta), sin(theta),0;
           -sin(theta), cos(theta),0;
                 0       ,    0       ,1];
             
        Rrel_LVLH(j,:)=(ROT*Rrel_I(j,:)')';
        Vrel_LVLH(j,:)=(ROT*(Vrel_I(j,:)'-cross([0;0;omega],Rrel_I(j,:)')))';
        
end
Rrel_puffo=zeros(N,3);
Vrel_puffo=zeros(N,3);

%PUFFO
for p=1:N
    t=t_cw(p);
    theta=omega*t;
        C=[ cos(theta), sin(theta),0;
              -sin(theta), cos(theta),0;
                 0       ,    0       ,1];
    Rrel_puffo(p,:)=(C*Rrel_I_rec(p,:)')';
         C_dot= omega*[-sin(theta), cos(theta),0
                          -cos(theta) -sin(theta),0
                          0, 0 , 0];
                      
    Vrel_puffo(p,:)=(C_dot*Rrel_I_rec(p,:)'+C*Vrel_I_rec(p,:)')';
                      
                      
end
 
% figure()
% plot3(Rho(:,1),Rho(:,2),Rho(:,3),'b',x_AN(:,1),x_AN(:,2),x_AN(:,3),'--k',Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or',x_AN(1,1),x_AN(1,2),x_AN(1,3),'og',x_AN(end,1),x_AN(end,2),x_AN(end,3),'or')
% title('!CW-INT! vs !CW-AN!')
% hold on
% grid on 


% figure()
% plot3(R_t(:,1),R_t(:,2),R_t(:,3),'--b',R_c(:,1),R_c(:,2),R_c(:,3),'--r',R_t(1,1),R_t(1,2),R_t(1,3),'og',R_t(end,1),R_t(end,2),R_t(end,3),'or',R_c(1,1),R_c(1,2),R_c(1,3),'og',R_c(end,1),R_c(end,2),R_c(end,3),'or')
% title('!CW-INT! vs !CIRC-diff! in LVLH ref')
% hold on
% grid on 

figure()
plot3(Rho(:,1),Rho(:,2),Rho(:,3),'--r',Rrel_LVLH(:,1),Rrel_LVLH(:,2),Rrel_LVLH(:,3),'--b',Rrel_LVLH(1,1),Rrel_LVLH(1,2),Rrel_LVLH(1,3),'og',Rrel_LVLH(end,1),Rrel_LVLH(end,2),Rrel_LVLH(end,3),'or',Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or')
title('!CW-INT! vs !CIRC-diff! in LVLH ref')
hold on
grid on 

figure()
plot3(Ni(:,1),Ni(:,2),Ni(:,3),'--r',Vrel_LVLH(:,1),Vrel_LVLH(:,2),Vrel_LVLH(:,3),'--b')
title('VELOCITY COMPARISON:!CW-INT! vs !CIRC-diff! in LVLH ref')
hold on
grid on 


figure()
plot3(Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'--b',Rrel_I_rec(:,1),Rrel_I_rec(:,2),Rrel_I_rec(:,3),'--r',Rrel_I(1,1),Rrel_I(1,2),Rrel_I(1,3),'og',Rrel_I(end,1),Rrel_I(end,2),Rrel_I(end,3),'or',Rrel_I_rec(1,1),Rrel_I_rec(1,2),Rrel_I_rec(1,3),'og',Rrel_I_rec(end,1),Rrel_I_rec(end,2),Rrel_I_rec(end,3),'or')
title('!CW-INT! vs !CIRC-diff! in IN ref')
hold on
grid on 
figure()
plot3(Vrel_I(:,1),Vrel_I(:,2),Vrel_I(:,3),'--b',Vrel_I_rec(:,1),Vrel_I_rec(:,2),Vrel_I_rec(:,3),'--r')
title('VELOCITY COMPARISON:!CW-INT! vs !CIRC-diff! in IN ref')
hold on
grid on 

figure()
plot3(Rrel_puffo(:,1),Rrel_puffo(:,2),Rrel_puffo(:,3),'--b',Rho(:,1),Rho(:,2),Rho(:,3),'--r')
title('R-puffo')
hold on
grid on 

figure()
plot3(Vrel_puffo(:,1),Vrel_puffo(:,2),Vrel_puffo(:,3),'--b',Ni(:,1),Ni(:,2),Ni(:,3),'--r')
title('V-puffo')
hold on
grid on 

%% PD-DOCKING control 

%paper "the Applicability of Pulsed Plasma Thrusters to Rendezvous and
%Docking of Cubesats" [continous trust is assumed]



%RANDOM chaser
R_0c=(R_0t-rand(3,1)*0.5*1e2);
V_0c=(V_0t+rand(3,1)*0.5*1e-1);




X_0t=[R_0t;V_0t];
X_0c=[R_0c;V_0c];
X_0t_c=[X_0t;X_0c];

%free motion
options = odeset('RelTol',1e-13,'AbsTol', 1e-30);
[t_int,Xt_c]= ode113('r2BP_t_c',t_end_CW,X_0t_c,options);

R_t=[Xt_c(:,1),Xt_c(:,2),Xt_c(:,3)];
R_c=[Xt_c(:,7),Xt_c(:,8),Xt_c(:,9)];
V_t=[Xt_c(:,4),Xt_c(:,5),Xt_c(:,6)];
V_c=[Xt_c(:,10),Xt_c(:,11),Xt_c(:,12)];

Rrel_I=R_c-R_t; %--> free motion

%PD control
kv=-0.05
kp=-7.5*1e-5

m_cubesat=5;
f_sat=2/5*1e-3*m_cubesat;% ref1: a_sat=4*1e-6N && m_cubesat=3Kg 
                         % ref2[PLASMA]: f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3 http://pepl.engin.umich.edu/thrusters/CAT.html
                         % ref3[COMMERCIAL]: http://www.busek.com/cubesatprop__main.htm
                         % ref4[NANOSAT]:f_sat=100*1e-6 N &&  Imp_duration= 2*1e-3 s && Isp= 50:100 s http://www.cubesatshop.com/index.php?page=shop.product_details&flypage=flypage.tpl&product_id=74&vmcchk=1&option=com_virtuemart&Itemid=65f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3

Kv=kv;
Kp=kp;
F_sat=f_sat
M_cubesat=m_cubesat

setGlobalKp(Kp)
setGlobalKv(Kv)
setGlobalf_sat(F_sat)
setGlobalm_cubesat(M_cubesat)

t_end_PD=t_end_CW;

% options = odeset('RelTol',1e-13,'AbsTol', 1e-30);
% [t_PD,X_PD]= ode113('r2BP_t_c_PD',t_end_PD,X_0t_c,options);

options = odeset('RelTol',1e-13,'AbsTol', 1e-45);
[t_PD,X_PD]= ode113('r2BP_t_c_PD_sat',t_end_PD,X_0t_c,options);


R_t_PD=[X_PD(:,1),X_PD(:,2),X_PD(:,3)];
R_c_PD=[X_PD(:,7),X_PD(:,8),X_PD(:,9)];
V_t_PD=[X_PD(:,4),X_PD(:,5),X_PD(:,6)];
V_c_PD=[X_PD(:,10),X_PD(:,11),X_PD(:,12)];



Rrel_I_PD=R_c_PD-R_t_PD;
Vrel_I_PD=V_c_PD-V_t_PD;

[Q,~]=size(t_PD);

Rrel_LVLH_PD=zeros(Q,3);
Vrel_LVLH_PD=zeros(Q,3);

%IN-->LVLH
for q=1:Q
    t=t_PD(q);
    theta=omega*t;
    ROT=[ cos(theta), sin(theta),0;
           -sin(theta), cos(theta),0;
                 0       ,    0       ,1];
             
        Rrel_LVLH_PD(q,:)=(ROT*Rrel_I_PD(q,:)')';
        Vrel_LVLH_PD(q,:)=(ROT*(Vrel_I_PD(q,:)'-cross([0;0;omega],Rrel_I_PD(q,:)')))';
        
end
% figure()
% plot3(R_t_PD(:,1),R_t_PD(:,2),R_t_PD(:,3),'--b',R_c_PD(:,1),R_c_PD(:,2),R_c_PD(:,3),'--r',R_t_PD(1,1),R_t_PD(1,2),R_t_PD(1,3),'og',R_t_(end,1),R_t(end,2),R_t(end,3),'or',R_c(1,1),R_c(1,2),R_c(1,3),'og',R_c(end,1),R_c(end,2),R_c(end,3),'or')
% title('!CW-INT! vs !CIRC-diff! in LVLH ref')
% hold on
% grid on 

disp('Rrel_I_PD, END:')
Rrel_I_PD(end,:)
disp('Rrel_LVLH_PD, END:')
Rrel_LVLH_PD(end,:)

disp('Vrel_I_PD, END:')
Vrel_I_PD(end,:)
disp('Vrel_I_PD, END:')
Vrel_LVLH_PD(end,:)

figure()
plot3(Rrel_I_PD(:,1),Rrel_I_PD(:,2),Rrel_I_PD(:,3),'--b',Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'--k',Rrel_I(1,1),Rrel_I(1,2),Rrel_I(1,3),'og',Rrel_I(end,1),Rrel_I(end,2),Rrel_I(end,3),'or',Rrel_I_PD(1,1),Rrel_I_PD(1,2),Rrel_I_PD(1,3),'og',Rrel_I_PD(end,1),Rrel_I_PD(end,2),Rrel_I_PD(end,3),'or')
title('Uncontrolloed (--k) VS PD (--b) in IN ref ')
hold on
grid on 

figure()
plot3(Rrel_I_PD(:,1),Rrel_I_PD(:,2),Rrel_I_PD(:,3),'--b',Rrel_I_PD(1,1),Rrel_I_PD(1,2),Rrel_I_PD(1,3),'og',Rrel_I_PD(end,1),Rrel_I_PD(end,2),Rrel_I_PD(end,3),'or')
title('controlled PD (--b) in IN ref ')
hold on
grid on 

figure()
plot3(Rrel_LVLH_PD(:,1),Rrel_LVLH_PD(:,2),Rrel_LVLH_PD(:,3),'--b',Rrel_LVLH_PD(1,1),Rrel_LVLH_PD(1,2),Rrel_LVLH_PD(1,3),'og',Rrel_LVLH_PD(end,1),Rrel_LVLH_PD(end,2),Rrel_LVLH_PD(end,3),'or')
title('controlloed PD (--b) in LVLH ref ')
hold on
grid on 

figure()
plot3(Rrel_LVLH_PD(:,1),Rrel_LVLH_PD(:,2),Rrel_LVLH_PD(:,3),'--r',Rrel_I_PD(:,1),Rrel_I_PD(:,2),Rrel_I_PD(:,3),'--b')
title('controlled PD (--b) in IN ref && (--r) in LVLH  ')
hold on
grid on 

%% LQR docking control

a_sat=f_sat/m_cubesat
Wu=[1/a_sat , 0,    0;
    0,      1/a_sat, 0;
    0, 0,       1/a_sat];

% Wu=[1 , 0,    0;
%     0,  1,    0;
%     0,  0,    1];
% 
Wz=eye(6);
Wz_red=eye(5);


Cz=1/2*m_cubesat*eye(6);
Cz=1/2*sqrt(1)*[0 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 1 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 1]+1/2*sqrt(Mu/(norm(R_T)^2))*[1 0 0 0 0 0;
                                                         0 0 0 0 0 0;
                                                         0 0 1 0 0 0;
                                                         0 0 0 0 0 0;
                                                         0 0 0 0 1 0;
                                                         0 0 0 0 0 0];% potential energy -[Mu/Rt-Mu/R_c]=-Mu/(R_t~R_c)^2*(R_c-R_t)=-Mu/R_t*[Rrel_I]
                                                     
 %need a weight to reduce control effort on {x,y} while can stay max on z
Wz=0.15   

Cz=1/2*sqrt(1)*[0    0     0    0      0    0;
                0 (1-Wz)/2 0    0      0    0;
                0    0     0    0      0    0;
                0    0     0  (1-Wz)/2 0    0;
                0    0     0    0      0    0;
                0    0     0    0      0   Wz]-1/2*sqrt(Mu/(norm(R_T)^2))*[(1-Wz)/2 0     0      0    0   0;
                                                                           0        0     0      0    0   0;
                                                                           0        0   (1-Wz)/2 0    0   0;
                                                                           0        0      0     0    0   0;
                                                                           0        0      0     0    Wz  0;
                                                                           0        0      0     0    0   0];
                                                                       
                                                                       
%need a weight to encrease control effort on {x,y,z} while can stay lower on  {x_dot,y_dot,z_dot}     

% Wv=0.1; 
% 
% Cz=1/2*sqrt(1)*(1+Wv)*[0    0     0    0      0    0;
%                        0 (1-Wz)/2 0    0      0    0;
%                        0    0     0    0      0    0;
%                        0    0     0  (1-Wz)/2 0    0;
%                        0    0     0    0      0    0;
%                        0    0     0    0      0   Wz]-1/2*sqrt(Mu/(norm(R_T)^2))* (1-Wv)*[(1-Wz)/2 0     0      0    0   0;
%                                                                                           0        0     0      0    0   0;
%                                                                                           0        0   (1-Wz)/2 0    0   0;
%                                                                                           0        0      0     0    0   0;
%                                                                                           0        0      0     0    Wz  0;
%                                                                                           0        0      0     0    0   0];


Cz_red=1/2*m_cubesat*eye(5);


Dz=zeros(6);

Dz_red=zeros(5);



A=[0     1    0     0   0   0;
  3*w^2  0    0    2*w  0   0;
   0     0    0     1   0   0;
   0   -2*w   0     0   0   0;
   0     0    0     0   0   1;
   0     0    0     0 -w^2  0];


A_red=[0     1    0   0   0;
    3*w^2   0   2*w  0   0;
       0   -2*w   0   0   0;
       0     0    0   0   1;
       0     0    0 -w^2  0];
   
B=[0 0 0;
   1 0 0; 
   0 0 0;
   0 1 0; 
   0 0 0; 
   0 0 1];

B_red=[0 0 0;
       1 0 0;
       0 1 0;
       0 0 0;
       0 0 1];
   
   
 rr_u=0.99;
 rr_z=0.01;
   
Q = Cz'*(rr_z*Wz)*Cz;
Q_red=Cz_red'*Wz_red*Cz_red;

R = rr_u*Wu;
R_red=rr_u*Wu;

K = lqr(A,B,Q,R);
K_red= lqr(A_red,B_red,Q_red,R_red);

K_lqr=-K

K_lqr_red=-K_red

setGlobalk_lqr(K_lqr);



options = odeset('RelTol',1e-13,'AbsTol', 1e-45);
[t_LQR,X_LQR]= ode113('r2BP_t_c_LQR_sat',t_end_PD,X_0t_c,options);

R_t_LQR=[X_LQR(:,1),X_LQR(:,2),X_LQR(:,3)];
R_c_LQR=[X_LQR(:,7),X_LQR(:,8),X_LQR(:,9)];
V_t_LQR=[X_LQR(:,4),X_LQR(:,5),X_LQR(:,6)];
V_c_LQR=[X_LQR(:,10),X_LQR(:,11),X_LQR(:,12)];

Rrel_I_LQR=R_c_LQR-R_t_LQR;
Vrel_I_LQR=V_c_LQR-V_t_LQR;

[S,~]=size(t_LQR);

Rrel_LVLH_LQR=zeros(S,3);
Vrel_LVLH_LQR=zeros(S,3);

for s=1:S
    t=t_LQR(s);
    theta=omega*t;
    ROT=[ cos(theta), sin(theta),0;
           -sin(theta), cos(theta),0;
                 0       ,    0       ,1];
             
        Rrel_LVLH_LQR(s,:)=(ROT*Rrel_I_LQR(s,:)')';
        Vrel_LVLH_LQR(s,:)=(ROT*(Vrel_I_LQR(s,:)'-cross([0;0;omega],Rrel_I_LQR(s,:)')))';
        
end

disp('Rrel_I_LQR, END:')
Rrel_I_LQR(end,:)
disp('Rrel_LVLH_LQR, END:')
Rrel_LVLH_LQR(end,:)

disp('Vrel_I_LQR, END:')
Vrel_I_LQR(end,:)
disp('Vrel_LVLH_LQR, END:')
Vrel_LVLH_LQR(end,:)


figure()
plot3(Rrel_I_LQR(:,1),Rrel_I_LQR(:,2),Rrel_I_LQR(:,3),'--b',Rrel_I_LQR(1,1),Rrel_I_LQR(1,2),Rrel_I_LQR(1,3),'og',Rrel_I_LQR(end,1),Rrel_I_LQR(end,2),Rrel_I_LQR(end,3),'or')
title('controlled LQR (--b) in IN ref ')
hold on
grid on 

figure()
plot3(Rrel_LVLH_LQR(:,1),Rrel_LVLH_LQR(:,2),Rrel_LVLH_LQR(:,3),'--b',Rrel_LVLH_LQR(1,1),Rrel_LVLH_LQR(1,2),Rrel_LVLH_LQR(1,3),'og',Rrel_LVLH_LQR(end,1),Rrel_LVLH_LQR(end,2),Rrel_LVLH_LQR(end,3),'or')
title('controlloed LQR (--b) in LVLH ref ')
hold on
grid on 


figure()
plot3(Rrel_I_LQR(:,1),Rrel_I_LQR(:,2),Rrel_I_LQR(:,3),'--b',Rrel_I_LQR(1,1),Rrel_I_LQR(1,2),Rrel_I_LQR(1,3),'og',Rrel_I_LQR(end,1),Rrel_I_LQR(end,2),Rrel_I_LQR(end,3),'or',Rrel_I_PD(:,1),Rrel_I_PD(:,2),Rrel_I_PD(:,3),'--r',Rrel_I_PD(1,1),Rrel_I_PD(1,2),Rrel_I_PD(1,3),'og',Rrel_I_PD(end,1),Rrel_I_PD(end,2),Rrel_I_PD(end,3),'or')
title('COMPARISON LQR (--b) VS PD (--r) in IN ref ')
hold on
grid on 
   
%% simulink test launch (perturbations included )

SIM_TEST=0; %<- 1 to run simulink

t_now=datetime('now');
t_2000=datetime(2000,1,1,00,00,00);
DAYs=juliandate(t_now)-juliandate(t_2000);

[E_Moon,~]=uplanet_mick( DAYs, 11);
[E,~]=uplanet_mick( DAYs, 3);

if SIM_TEST==1 
out=sim('test_hill_sim');

figure()
plot3(Rrel_I_sim.signals.values(:,1),Rrel_I_sim.signals.values(:,2),Rrel_I_sim.signals.values(:,3))
title('sdfdzf')
grid on
hold on
end

% !!! Note: Moon perturbation O(1e-6) comparable to the actuator's limit
% O(1e-4 : 1e-6) !!!
