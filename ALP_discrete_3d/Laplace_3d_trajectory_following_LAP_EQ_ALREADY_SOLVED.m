close all
clear all

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
x_obs=-20
y_obs=-20
z_obs=0



lenght_ISS=25; %--> along x
width_ISS=15; %--> along y
height_ISS=10;%--> along z


%[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,width_ISS,height_ISS);

%% Initial Conditions & Integration time 

% !!!! CRUCIAL TEST: by  givinig initial conditions such that x_obs_2d is
% reached no trajectory corrections are made. ====> OK !!! **

rho_0_3d=[5;40;0];

tau=1/5*T_orb;
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


Ni_0_3d=N_tau_inv*(X_obs_3d-M_tau*rho_0_3d);
X_0_3d=[rho_0_3d;Ni_0_3d];

ci_0_3d=[rho_0_3d;Ni_0_3d];



%% Free Drift !! UNforced Problem!!

t_end=tau;

options = odeset('RelTol',1e-13,'AbsTol', 1e-12);
[t_CW_FD,r_v_CW_FD]= ode113('Chol_Wilt_Hill',t_end,ci_0_3d,options);%(plot at the end of the script)



%% Obstacle (plot only) !! always refer to Discrete_Laplacae_3d_CI_obstacle.m for data: step_grid, L_box ,....!!


step_grid=1;

L_box=40;

setGlobalL_box_x(L_box)
setGlobalL_box_y(L_box)
setGlobalL_box_z(L_box)



X_max=L_box;
Y_max=L_box;
Z_max=L_box;


R1_boundedx=-X_max:step_grid:X_max;
R1_boundedy=-Y_max:step_grid:Y_max;
R1_boundedz=-Z_max:step_grid:Z_max;

[~,N]=size(R1_boundedx);
[~,M]=size(R1_boundedy);
[~,P]=size(R1_boundedz);

if x_obs >=0
i_obs=floor(x_obs/step_grid)+floor((X_max/step_grid))
else
i_obs=floor((x_obs+floor(X_max/step_grid)*step_grid)/step_grid)
end

if y_obs >=0
j_obs=floor(y_obs/step_grid)+floor((Y_max/step_grid))
else
j_obs=floor((y_obs+floor(Y_max/step_grid)*step_grid)/step_grid)
end

if z_obs >=0
k_obs=floor(z_obs/step_grid)+floor((Z_max/step_grid))
else
k_obs=floor((z_obs+floor(Z_max/step_grid)*step_grid)/step_grid)
end

Phi0_3d_plot_plus=zeros(M,N);
Phi0_3d_plot_minus=zeros(M,N);

for i=1:N
    for j=1:M
        
        x=R1_boundedx(i);
        y=R1_boundedy(j);
       
        
       
          % plot 
         if y>=0
                
                  if abs(x)<=lenght_ISS && abs(y) <= height_ISS
                
                Phi0_3d_plot_plus(i,j)=height_ISS/2;
                Phi0_3d_plot_minus(i,j)=-height_ISS/2;
                
                  elseif abs(x) <=lenght_ISS   && abs(x) >=3/2*width_ISS && abs(y)<=3/2*width_ISS
                     
                  Phi0_3d_plot_plus(i,j)=height_ISS/2;
                  Phi0_3d_plot_minus(i,j)=-height_ISS/2;
                  end
                
            elseif y<0
                
                
                if abs(x)<=height_ISS/2 && abs(y) <= 2*width_ISS
                
                Phi0_3d_plot_plus(i,j)=height_ISS/2;
                Phi0_3d_plot_minus(i,j)=-height_ISS/2;
                
                end
                
               
        end
            
   end
end

%% ALP method

load('Phi_3d') % X_goal cannot be changed
               % Rho_0 can be changed
               
               
               
              
              


setGlobalstep_grid(step_grid)
setGlobalPhi_3d(Phi_3d)

Phi_2d_xy=Phi_3d(:,:,(P-1)/2);


figure()
mesh(R1_boundedx,R1_boundedy,Phi_2d_xy')
colormap summer
hold on
plot3(x_obs,y_obs,0,'*b')
title('claplace potential x-y')
grid on 


figure()
contour(R1_boundedx,R1_boundedy,Phi_2d_xy')
colormap summer
hold on
plot3(x_obs,y_obs,0,'*b')
title('claplace potential x-y')
grid on 



% alpha IsSue
alpha_v=1*norm(rho_0_3d-X_obs_3d)/(1*t_end)
   % alpha_v=norm(Ni_0_2d)*1e-4;
 setGlobalalpha_v(alpha_v );

%step_time to reduce control action
step_time=100;
setGlobalstep_time(step_time)
flag_time=0;
setGlobalflag_time(flag_time )
time_before=0;
setGlobaltime_before(time_before)
time_last_fire=0;
setGlobaltime_last_fire(time_last_fire)
setGloaltau_or(tau)
setGloalx_obs(X_obs_3d)

dt=1;
t_span=[0:dt:2*tau]';
X_0_3d=[rho_0_3d;Ni_0_3d];

% [Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done,i_obs]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_2d_step_time(dt,t_span,X_0_3d,X_obs_3d);
  [rho_an,ni_an,DV_v,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_3d_step_time(dt,t_span,ci_0_3d);


figure()
hold on
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus')
colormap summer
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_minus')
colormap summer
plot3(rho_an(:,1),rho_an(:,2),rho_an(:,3),'--b',rho_an(1,1),rho_an(1,2),rho_an(1,3),'og',rho_an(end,1),rho_an(end,2),rho_an(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('collision Avoidance with discrete gradient & analytical trajectory propagation')
axis equal
grid on 

figure()
hold on
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus')
colormap summer
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_minus')
colormap summer
hold on
plot3(r_v_CW_FD(:,1),r_v_CW_FD(:,2),r_v_CW_FD(:,3),'--b',r_v_CW_FD(1,1),r_v_CW_FD(1,2),r_v_CW_FD(1,3),'og',r_v_CW_FD(end,1),r_v_CW_FD(end,2),r_v_CW_FD(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('free drift')
axis equal
grid on 



%results

DV_tot=sum(norm(DV_v));
Goal_error=norm(rho_an(end,:)'-X_obs_3d)
Total_Impulse=count_impulse_done

DV_tot_obs=sum(norm(DV_v(1:i_obs,:)))
TOF_obs=i_obs*dt

Set_Point_traj_2=[t_span,Rho_an,Ni_an]

save('Set_Point_traj_2','Set_Point_traj_2')
