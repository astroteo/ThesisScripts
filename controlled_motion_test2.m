clear all
close all

%% desired observation point in LVLH reference frame.
x_obs=20
y_obs=70
z_obs=0

L_box=2;

X_obs=[x_obs;y_obs;z_obs];%-> observation point given in LVLH r.f.


%% target & chaser initial conditions
%orbit
Mu=(3.9856e+14);
r_Earth=6380*10^3;

h_orb=400*10^3;
%i_orb=51.65*pi/180;
i_orb=0;
om_orb=0;
OM_orb=0;
%theta0_orb=90*pi/180;
theta0_orb=0;
a_orb=r_Earth+h_orb;
e_orb=0;


T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)



omega=sqrt(Mu/((r_Earth+h_orb)^3));
w=omega;
setGlobalW(w);

d_TC=200;


ROT_plane0=[cos(om_orb)*cos(OM_orb)-sin(om_orb)*cos(i_orb)*sin(OM_orb) -sin(om_orb)*cos(OM_orb)-cos(om_orb)*cos(i_orb)*sin(OM_orb) sin(i_orb)*sin(OM_orb)
            cos(om_orb)*sin(OM_orb)+sin(om_orb)*cos(i_orb)*cos(OM_orb) -sin(om_orb)*sin(OM_orb)+cos(om_orb)*cos(i_orb)*cos(OM_orb) -sin(i_orb)*cos(OM_orb)
              sin(om_orb)*sin(i_orb)                                                 cos(om_orb)*sin(i_orb)                                 cos(i_orb)]  ;

ROT_theta0=[ cos(theta0_orb), sin(theta0_orb),0;
           -sin(theta0_orb), cos(theta0_orb),0;
                 0       ,    0       ,1];
             
ROT_0=ROT_plane0*ROT_theta0;
setGlobalROTplane_0(ROT_plane0);



%target:
[R_0t,V_0t]=OrbPar2RV_MY(a_orb,e_orb,i_orb,OM_orb,om_orb,theta0_orb,Mu);
w_v=cross(R_0t,V_0t)/(norm(R_0t)^2);
W_v=w_v;
setGlobalw_v(W_v);

%chaser:
% R_0c=(R_0t+rand(3,1)*0.1*1e2); 
% V_0c=(V_0t-rand(3,1)*0.1*1e-2);

% R_0c=(R_0t+X_obs-[rand(1,1)*0.1*1e1;rand(1,1)*0.1*1e1;0]); %-> pert. only in plane: circular orbit
% V_0c=(V_0t+[rand(1,1)*0.1*1e-2;rand(1,1)*0.1*1e-3;0]);

%  R_0c=R_0t+X_obs; %--> maintaince (stationary keeping) | valid in this form since LVLH && I coincident @ t=0.|
%  V_0c=V_0t+zeros(3,1);

% R_0c=(R_0t+[70;120;0]); 
%  V_0c=(V_0t-rand(3,1)*0.1*1e-2);

%  R_0c=(R_0t+[-140;-240;0]); 
%  V_0c=(V_0t-rand(3,1)*0.1*1e-2);

%  R_0c=(R_0t+[-70;120;0]); 
%  V_0c=(V_0t-rand(3,1)*0.1*1e-2);

%   R_0c=(R_0t+[-150;120;10]); 
%   V_0c=(V_0t-rand(3,1)*0.1*1e-2);

  R_0c=(R_0t+[-70;120;30]); 
  V_0c=(V_0t-rand(3,1)*0.1*1e-2);

  R_0c=(R_0t+[-70;120;-30]); 
  V_0c=(V_0t-rand(3,1)*0.1*1e-2);




%Relative initial position:
Rrel_I0=R_0c-R_0t;
Vrel_I0=V_0c-V_0t;

Rho_0=ROT_0'*Rrel_I0;
Ni_0=ROT_0'*(Vrel_I0-cross(w_v,Rrel_I0));



%% actuators & controllers design

%actuator & cubesat.
m_cubesat=5;
f_sat=1*1e-3*m_cubesat;% ref1: a_sat=4*1e-6N && m_cubesat=3Kg 
                         % ref2[PLASMA]: f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3 http://pepl.engin.umich.edu/thrusters/CAT.html
                         % ref3[COMMERCIAL]: http://www.busek.com/cubesatprop__main.htm
                         % ref4[NANOSAT]:f_sat=100*1e-6 N &&  Imp_duration= 2*1e-3 s && Isp= 50:100 s http://www.cubesatshop.com/index.php?page=shop.product_details&flypage=flypage.tpl&product_id=74&vmcchk=1&option=com_virtuemart&Itemid=65f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3
F_sat=f_sat
M_cubesat=m_cubesat

setGlobalf_sat(F_sat)
setGlobalm_cubesat(M_cubesat)




%LQR design
a_sat=f_sat/m_cubesat
Wu=[1/a_sat , 0,    0;
    0,      1/a_sat, 0;
    0, 0,       1/a_sat];

Wz=eye(6);

wz=0.999 %<-need a weight to reduce control effort on {x,y} while can stay max on z  

rr_p=1
Cz=1/2*sqrt(1)*[0    0     0    0      0    0;
                0 (1+wz/2) 0    0      0    0;
                0    0     0    0      0    0;
                0    0     0  (1+wz/2) 0    0;
                0    0     0    0      0    0;
                0    0     0    0      0  (1-wz)]-1/2*sqrt(Mu/(norm(R_0t)^2))*[(1+wz/2) 0     0      0    0   0;
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


 rr_u=0.999;
 rr_z=0.001;
   
Q = Cz'*(rr_z*Wz)*Cz;


R = rr_u*Wu;


K = lqr(A,B,Q,R);

K_lqr=-K

setGlobalk_lqr(K_lqr);

%PD design
kv=-0.05*1000
kp=-7.5*1e-5*1000

Kv=kv;
Kp=kp;

setGlobalKp(Kp)
setGlobalKv(Kv)

%% motion to observation point:

setGlobalx_obs(X_obs)
setGloball_box(L_box)


t_end=1/0.5*T_orb

 X_0t=[R_0t;V_0t];
  X_0c=[R_0c;V_0c];
   X_0t_c=[X_0t;X_0c];
   
options = odeset('RelTol',1e-13,'AbsTol', 1e-6);
[t_LQR_FT,X_LQR_FT]= ode113('r2BP_t_c_LQR_SP_new',t_end,X_0t_c,options);

R_t_LQR_FT=[X_LQR_FT(:,1),X_LQR_FT(:,2),X_LQR_FT(:,3)];
R_c_LQR_FT=[X_LQR_FT(:,7),X_LQR_FT(:,8),X_LQR_FT(:,9)];
V_t_LQR_FT=[X_LQR_FT(:,4),X_LQR_FT(:,5),X_LQR_FT(:,6)];
V_c_LQR_FT=[X_LQR_FT(:,10),X_LQR_FT(:,11),X_LQR_FT(:,12)];

Rrel_I_LQR_FT=R_c_LQR_FT-R_t_LQR_FT;
Vrel_I_LQR_FT=V_c_LQR_FT-V_t_LQR_FT;

[F,~]=size(t_LQR_FT);

Rrel_LVLH_LQR_FT=zeros(F,3);
Vrel_LVLH_LQR_FT=zeros(F,3);

i_v=[1;0;0];
j_v=[0;1;0];

for f=1:F
    
    t=t_LQR_FT(f);
    
    RT=norm(R_t_LQR_FT(f,:));
    theta=-acos(dot(i_v,(R_t_LQR_FT(f,:)/RT)));
    
      if dot(j_v,(R_t_LQR_FT(f,:)/RT)) <0 
        theta=2*pi-theta;
      end

    
    ROT_theta=[ cos(theta), sin(theta),    0;
               -sin(theta), cos(theta),    0;
                 0       ,    0       ,    1];


    ROT=ROT_plane0*ROT_theta;
    
    Rrel_LVLH_LQR_FT(f,:)=(ROT'*(R_c_LQR_FT(f,:)-R_t_LQR_FT(f,:))')';
    Vrel_LVLH_LQR_FT(f,:)=(ROT'*((V_c_LQR_FT(f,:)-V_t_LQR_FT(f,:))'-cross(w_v,(R_c_LQR_FT(f,:)-R_t_LQR_FT(f,:))')))';%!! TAKE CARE FOR THIS CROSS-PRODUCT !!
   
    
end

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




