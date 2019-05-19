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
[R_0t,V_0t]=OrbPar2RV_MY(a_orb,e_orb,i_orb,OM_orb,om_orb,theta0_orb,Mu);
w_v=cross(R_0t,V_0t)/(norm(R_0t)^2);

%chaser:
R_0c=(R_0t+rand(3,1)*0.25*1e2); 
V_0c=(V_0t+rand(3,1)*0.5*1e-1);

 ROT_plane0=[cos(om_orb)*cos(OM_orb)-sin(om_orb)*cos(i_orb)*sin(OM_orb) -sin(om_orb)*cos(OM_orb)-cos(om_orb)*cos(i_orb)*sin(OM_orb) sin(i_orb)*sin(OM_orb)
            cos(om_orb)*sin(OM_orb)+sin(om_orb)*cos(i_orb)*cos(OM_orb) -sin(om_orb)*sin(OM_orb)+cos(om_orb)*cos(i_orb)*cos(OM_orb) -sin(i_orb)*cos(OM_orb)
              sin(om_orb)*sin(i_orb)                                                 cos(om_orb)*sin(i_orb)                                 cos(i_orb)]  ;


%% Desired Observation Point & Obstacle (ISS ellipsoid)
x_obs=-50+3
y_obs=50+3
z_obs=20



lenght_ISS=100; %--> along y
width_ISS=50; %--> alng z
height_ISS=30;


[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,width_ISS,height_ISS);

%% Initial Conditions & Integration time 

% !!!! CRUCIAL TEST: by  givinig initial conditions such that x_obs_2d is
% reached no trajectory corrections are made. ====> OK !!! **

Rho_0_3d=[-40;-50;-20];

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
Ni_0_3d_CT=[Ni_0_3d(1);Ni_0_3d(2);Ni_0_3d(3)]; %** CRUCIAL TEST
Ni_0_3d_RAND=[-rand(1)*1e-1; -rand(1)*1e-1;-rand(1)*1e-1]; %% GENERIC STUFF
Ni_0_3d=Ni_0_3d_CT;





%% Free Drift !! UNforced Problem!!

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
grid on 




%% Position Potential PARAMETERS


m1=1;

%m1=lenght_ISS/lenght_ISS;


m2=1;
%m2=width_ISS/lenght_ISS;

m3=1;
%m2=height_ISS/lenght_ISS;

M_3d=[m1,0, 0;
      0, m2,0;
      0, 0, m3];


%% Spherical Harmonic 2D Potential PARAMETES

X_coll_3d=zeros(3,1);

L_max=norm(Rho_0_3d)

lambda_1_paper_Rich=1/2*(Rho_0_3d./L_max-X_coll_3d./L_max)'*(M_3d*(Rho_0_3d./L_max-X_coll_3d./L_max));
lambda_1_paper=1e6;
lambda_1_logic=1/2*([L_max;L_max;L_max]./L_max-X_coll_3d./L_max)'*(M_3d*([L_max;L_max;L_max]/L_max-X_coll_3d./L_max));
lambda_1_MY=lambda_1_paper;
lambda_1_illogic=1e6;

lambda_1=1*lambda_1_paper;

safety_radius=0; % narrower potential ==> steeper gradient ==> safer maneuver

lambda_2_MY=((height_ISS+safety_radius)/2)^2/(0.1*L_max);
lambda_2_paper=1e-4;
lambda_2=lambda_2_MY;

lambda_v=[lambda_1;lambda_2];

n1=(lenght_ISS+safety_radius)/lenght_ISS;
n2=(width_ISS+safety_radius)/lenght_ISS;
n3=(height_ISS+safety_radius)/lenght_ISS;




N_3d=[n1,0,0;
      0,n2,0;
      0,0,n3]; 


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
setGlobalalpha_v(alpha_v );

Perform_Integration=0;
if Perform_Integration==1
%integration time
t_end_Lap=t_end+3*T_orb

%integration
options = odeset('RelTol',1e-7,'AbsTol', 1e-1);
[t_CW,x_CW]= ode113('Chol_Wilt_Hill_Laplace_2d_analytic',t_end_Lap,X_0_2d,options);%<--[ISSUE: n° of correction strictly dependant on ODE's precision]


%plot
[CW,~]=size(t_CW);
zebra=zeros(CW,1);

figure()
surf(x_ISS, y_ISS,z_ISS,'facecolor','g')
hold on
plot3(x_CW(:,1),x_CW(:,2),x_CW(:,3),'--b',x_CW(1,1),x_CW(1,2),x_CW(1,3),'og',x_CW(end,1),x_CW(end,2),x_CW(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('collision Avoidance with analytical gradient')
grid on 

end


%% GNC by Laplace Artificial potential CONTROL EFFORT AT DISCRETE TIMES

%parameters to compute analytical gradient
setGloball_max(L_max)
setGlobalM_3d(M_3d)
setGlobalN_3d(N_3d)
setGloballambda_v(lambda_v)
setGlobalx_obs_3d(X_obs_3d )

%step_time to reduce control action
step_time=4;
setGlobalstep_time(step_time)
flag_time=0;
setGlobalflag_time(flag_time )
time_before=0;
setGlobaltime_before(time_before)
time_last_fire=0;
setGlobaltime_last_fire(time_last_fire)



% alpha IsSue
alpha_v_paper_rich=norm(Rho_0_3d-X_obs_3d)/(1*t_end)
alpha_v_paper=norm(Ni_0_3d)*1e-4;
f_sat=1*1e-3
m_cubesat=5
alpha_v_sat=(f_sat/m_cubesat)*step_time; % to be used with dt=step_time
alpha_v=alpha_v_paper_rich
setGlobalalpha_v(2*alpha_v );

% %integration time
t_end_Lap_st=t_end+10*T_orb

t_span_Lap_st=linspace(0,t_end_Lap_st,10000);


%% Full Analytical Computation [discrete time control action]

dt=4;
t_span=[0:dt:30*tau]';

[Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Full_Analytic_Laplace_3d_step_time(dt,t_span,X_0_3d);

count_impulse_time
count_impulse_done
[size_ODE,~]=size(t_CW_FD)
%[size_ODE_ctrl,~]=size(t_CW)
[size_SPAN,~]=size(t_span)

figure()
surf(x_ISS, y_ISS,z_ISS,'facecolor','g')
hold on
plot3(Rho_an(:,1),Rho_an(:,2),Rho_an(:,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(end,1),Rho_an(end,2),Rho_an(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('collision Avoidance with analytical gradient & analytical trajectory propagation')
hold on
grid on 

Set_Point_traj_6_3d=[t_span,Rho_an,Ni_an]

save('Set_Point_traj_6_3d')


disp('X_obs_3d:')
X_obs_3d'