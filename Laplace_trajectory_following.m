clear all
close all

%% ISS Orbit random chaser.

%orbit
Mu=(3.9856e+14);
r_Earth=6380*10^3;

h_orb=400*10^3;
i_orb=0;
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
[R_0t,V_0t]=OrbPar2rv(a_orb,e_orb,i_orb,OM_orb,om_orb,theta0_orb,Mu);
w_v=cross(R_0t,V_0t)/(norm(R_0t)^2);


 ROT_plane0=[cos(om_orb)*cos(OM_orb)-sin(om_orb)*cos(i_orb)*sin(OM_orb) -sin(om_orb)*cos(OM_orb)-cos(om_orb)*cos(i_orb)*sin(OM_orb) sin(i_orb)*sin(OM_orb)
            cos(om_orb)*sin(OM_orb)+sin(om_orb)*cos(i_orb)*cos(OM_orb) -sin(om_orb)*sin(OM_orb)+cos(om_orb)*cos(i_orb)*cos(OM_orb) -sin(i_orb)*cos(OM_orb)
              sin(om_orb)*sin(i_orb)                                                 cos(om_orb)*sin(i_orb)                                 cos(i_orb)]  ;


%% Trajectory Load:

load('t_traj1')
load('xxx_traj1')
load('y_traj1')
load('vx_traj1')
load('vy_traj1')

[SP,~]=size(t_traj1);

Set_Point=[t_traj1,xxx_traj1,y_traj1,zeros(SP,1),vx_traj1,vy_traj1,zeros(SP,1)];

%% Chaser Initial Conditions
Rho_0=Set_Point(1,2:4)';
Ni_0=Set_Point(1,5:7)';

R_0c=R_0t+Rho_0; % assuming: i_orb=0 && theta0_orb=0
V_0c=V_0t+Ni_0+cross(w_v,Rho_0);% assuming: i_orb=0 && theta0_orb=0

%% ODE [IN integration]
setGlobalSet_Point(Set_Point )

