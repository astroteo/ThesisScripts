clear all
close all

%% ISS Orbit random chaser.

%orbit
Mu=(3.9856e+14);
r_Earth=6380*10^3;

h_orb=400*10^3;
i_orb=51.65*pi/180;
%i_orb=0;
om_orb=0;
OM_orb=0;
theta0_orb=0;
a_orb=r_Earth+h_orb;
e_orb=0;


T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)

omega=sqrt(Mu/((r_Earth+h_orb)^3));
w=omega;
setGlobalW(w);


%target:
[R_0target,V_0target]=OrbPar2RV_MY(a_orb,e_orb,i_orb,OM_orb,om_orb,theta0_orb,Mu);
w_v=cross(R_0target,V_0target)/(norm(R_0target)^2);
setGlobalw_v(w_v);

 ROT_plane0=[cos(om_orb)*cos(OM_orb)-sin(om_orb)*cos(i_orb)*sin(OM_orb) -sin(om_orb)*cos(OM_orb)-cos(om_orb)*cos(i_orb)*sin(OM_orb) sin(i_orb)*sin(OM_orb)
            cos(om_orb)*sin(OM_orb)+sin(om_orb)*cos(i_orb)*cos(OM_orb) -sin(om_orb)*sin(OM_orb)+cos(om_orb)*cos(i_orb)*cos(OM_orb) -sin(i_orb)*cos(OM_orb)
              sin(om_orb)*sin(i_orb)                                                 cos(om_orb)*sin(i_orb)                                 cos(i_orb)]  ;
          
ROT_plane=eye(3);
setGlobalROT_plane(ROT_plane);


 %% Desired Observation Point & Obstacle (ISS ellipsoid)
x_obs=80+3
y_obs=60+3
z_obs=20


x_coll=0;
y_coll=0;
z_coll=0;

lenght_ISS=100; %--> along y
width_ISS=50; %--> alng z
height_ISS=30;


[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,width_ISS,height_ISS);

%%  Initial Conditions
 
Rho_0_3d=[0;-100;0];
tau=1/3*T_orb;

Ct=cos(w*tau);
St=sin(w*tau);

X_obs_3d=[x_obs;y_obs;z_obs];

det=(4*St)/(w^3) -(8*Ct*St)/(w^3) +(4*Ct^2*St)/(w^3)+(4*St^3)/(w^3) -(3*St^2*tau)/(w^2);

N_tau_inv=1/det*[(4*St^2)/(w^2)-(3*St*tau)/w,     -((2*St)/(w^2))+(2*Ct*St)/(w^2),                        0;
                 (2*St)/(w^2)-(2*Ct*St)/(w^2),              St^2/(w^2),                                   0;
                           0,                                  0,              4/(w^2)-(8*Ct)/(w^2)+(4*Ct^2)/(w^2)+(4*St^2)/(w^2)-(3*St*tau)/w];
                   
M_tau=[-3*Ct+4,        0,   0;
        6*St-6*w*tau,  1,   0;                   
               0    ,  0,  Ct];


Ni_0_3d=N_tau_inv*(X_obs_3d-M_tau*Rho_0_3d);


%% Free Drift

t_end=tau;
X_0_3d=[Rho_0_3d;Ni_0_3d];

options = odeset('RelTol',1e-13,'AbsTol', 1e-12);
[t_CW_FD,x_CW_FD]= ode113('Chol_Wilt_Hill',t_end,X_0_3d,options);%<--[]

[CW_FD,~]=size(t_CW_FD);
zebra_FD=zeros(CW_FD,1);

figure()
surf(x_ISS, y_ISS,z_ISS,'facecolor','g')
hold on
plot3(x_CW_FD(:,1),x_CW_FD(:,2),x_CW_FD(:,3),'--b',x_CW_FD(1,1),x_CW_FD(1,2),x_CW_FD(1,3),'og',x_CW_FD(end,1),x_CW_FD(end,2),x_CW_FD(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('Free Drift')
grid on 


          
%% Chaser Initial Conditions: LVLH-->IN


R_0chaser=R_0target+ROT_plane*Rho_0_3d; % assuming: theta0_orb=0
V_0chaser=V_0target+ROT_plane*(Ni_0_3d+cross([0;0;w],Rho_0_3d));% assuming: theta0_orb=0


%% Position Potential PARAMETERS
m1=1;
 m1=lenght_ISS/lenght_ISS;
m2=1;
 m2=width_ISS/lenght_ISS;
m3=1;
 m3=height_ISS/lenght_ISS;

M_3d=[m1,0, 0;
      0, m2,0;
      0, 0, m3];


%% Spherical Harmonic 2D Potential PARAMETES

X_coll_3d=zeros(3,1);

L_max=norm(Rho_0_3d)

%lambda_2: Potential Skewness
safety_radius=0; % narrower potential ==> steeper gradient ==> safer maneuver

lambda_2_MY=((height_ISS+safety_radius))^2/(1*L_max);
lambda_2_paper=1e2;
lambda_2_2d= 253.1250;
lambda_2_florida=(1e-2)^-1;% paper defines sigma_k=1/lambda_2_k

lambda_2=lambda_2_MY;

%lambda_1: Potential Width
lambda_1_paper_Rich=1/2*(Rho_0_3d./L_max-X_coll_3d./L_max)'*(M_3d*(Rho_0_3d./L_max-X_coll_3d./L_max));
lambda_1_paper=1e0;
lambda_1_logic=1/2*([L_max;L_max;L_max]./L_max-X_coll_3d./L_max)'*(M_3d*([L_max;L_max;L_max]/L_max-X_coll_3d./L_max));
lambda_1_MY=lambda_1_paper;
lambda_1_illogic=1e6;
lambda_1_2d=1.2656;
lambda_1_florida=exp((lenght_ISS/L_max)^2*lambda_2)*((norm(Rho_0_3d)/L_max)^2-1/2*(lenght_ISS/L_max)^2);
lambda_1=lambda_1_florida;




lambda_v=[lambda_1;lambda_2];

n1=(lenght_ISS+safety_radius)/lenght_ISS;
n2=(width_ISS+safety_radius)/lenght_ISS;
n3=(height_ISS+safety_radius)/lenght_ISS;

N_3d=[n1,0,0;
      0,n2,0;
      0,0,n3]; 
  
%N_3d=eye(3)


%% GNC by Laplace Aritificial Potential [ANALYCAL GRADIENT]

%parameters to compute analytical gradient
setGloball_max(L_max)
setGlobalM_3d(M_3d)
setGlobalN_3d(N_3d)
setGloballambda_v(lambda_v)
setGlobalx_obs_3d(X_obs_3d )


% alpha IsSue
alpha_v=norm(Rho_0_3d-X_obs_3d)/t_end
% alpha_v=norm(Ni_0_2d)*1e-4;
setGlobalalpha_v(1*alpha_v );


%% GNC by Laplace Artificial potential CONTROL EFFORT AT DISCRETE TIMES

%parameters to compute analytical gradient
setGloball_max(L_max)
setGlobalM_3d(M_3d)
setGlobalN_3d(N_3d)
setGloballambda_v(lambda_v)
setGlobalx_obs_3d(X_obs_3d )

%step_time to reduce control action
step_time=100;
setGlobalstep_time(step_time)
flag_time=0;
setGlobalflag_time(flag_time )
time_before=0;
setGlobaltime_before(time_before)
time_last_fire=0;
setGlobaltime_last_fire(time_last_fire)
time_before_Sp=0;
setGlobaltime_before_Sp(time_before_Sp)
flag_Sp=1;
setGlobalflag_Sp(flag_Sp)

% alpha IsSue
alpha_v_paper_Rich=norm(Rho_0_3d-X_obs_3d)/(1*t_end)
alpha_v_paper=norm(Ni_0_3d)*1e-4;
f_sat=1*1e-3
m_cubesat=5
alpha_v_sat=(f_sat/m_cubesat)*step_time; % to be used with dt=step_time
alpha_v_MY=alpha_v_paper_Rich
alpha_v=alpha_v_MY
setGlobalalpha_v(1*alpha_v );

% %integration time
t_end_Lap_sp=t_end+5*T_orb




%% Trajectory of  Design 

dt=1;
t_span=[0:dt:t_end_Lap_sp]';

[Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Full_Analytic_Laplace_3d_step_time_New(dt,t_span,X_0_3d);

figure()
surf(x_ISS, y_ISS,z_ISS,'facecolor','g')
hold on
plot3(Rho_an(:,1),Rho_an(:,2),Rho_an(:,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(end,1),Rho_an(end,2),Rho_an(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('collision Avoidance with analytical gradient & analytical trajectory propagation')
grid on 


%% actuators & controllers design

%actuator & cubesat.
m_cubesat=5;
f_satur=250*1e-3;% ref1: a_sat=4*1e-6N && m_cubesat=3Kg 
                         % ref2[PLASMA]: f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3 http://pepl.engin.umich.edu/thrusters/CAT.html
                         % ref3[COMMERCIAL]: http://www.busek.com/cubesatprop__main.htm
                         % ref4[NANOSAT]:f_sat=100*1e-6 N &&  Imp_duration= 2*1e-3 s && Isp= 50:100 s http://www.cubesatshop.com/index.php?page=shop.product_details&flypage=flypage.tpl&product_id=74&vmcchk=1&option=com_virtuemart&Itemid=65f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3
F_sat=f_satur
M_cubesat=m_cubesat

setGlobalf_sat(F_sat)
setGlobalm_cubesat(M_cubesat)




%LQR design
a_sat=f_satur/m_cubesat
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
                0    0     0    0      0  (1-wz)]-1/2*sqrt(Mu/(norm(R_0target)^2))*[(1+wz/2) 0     0      0    0   0;
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


 rr_u=0.9999; % real_act: 0.9999
 rr_z=0.0001;          %:0.00001
   
Q = Cz'*(rr_z*Wz)*Cz;


R = rr_u*Wu;


K = lqr(A,B,Q,R);

K_lqr=-K

setGlobalk_lqr(K_lqr);

%% RealTime 

X_0target=[R_0target;V_0target];
  X_0chaser=[R_0chaser;V_0chaser];
   X_0tar_cha=[X_0target;X_0chaser];
   
t_end_FT=t_end_Lap_sp;

options = odeset('RelTol',1e-13,'AbsTol', 1e-6);
[t_LQR_FT,X_LQR_FT]= ode113('r2BP_t_c_LQR_sat_SP_RealTime',t_end_FT,X_0tar_cha,options);

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


    ROT=ROT_plane*ROT_theta;
    
    Rrel_LVLH_LQR_FT(f,:)=(ROT'*(R_c_LQR_FT(f,:)-R_t_LQR_FT(f,:))')';
    Vrel_LVLH_LQR_FT(f,:)=(ROT'*((V_c_LQR_FT(f,:)-V_t_LQR_FT(f,:))'-cross(w_v,(R_c_LQR_FT(f,:)-R_t_LQR_FT(f,:))')))';%!! TAKE CARE FOR THIS CROSS-PRODUCT !!
   
    
end



figure()
plot3(Rrel_LVLH_LQR_FT(:,1),Rrel_LVLH_LQR_FT(:,2),Rrel_LVLH_LQR_FT(:,3),'--b',Rrel_LVLH_LQR_FT(1,1),Rrel_LVLH_LQR_FT(1,2),Rrel_LVLH_LQR_FT(1,3),'og',Rrel_LVLH_LQR_FT(end,1),Rrel_LVLH_LQR_FT(end,2),Rrel_LVLH_LQR_FT(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('controlloed LQR along reference trajectory in LVLH ref ')
hold on
surf(x_ISS, y_ISS,z_ISS,'facecolor','g')
axis equal
hold on
grid on
%% 2d Potentials Plot

