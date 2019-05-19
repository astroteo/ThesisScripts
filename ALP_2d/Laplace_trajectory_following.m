
clearvars -global
close all
clear all


%% ISS Orbit random chaser.

%orbit
Mu=(3.9856e+14);
r_Earth=6380*10^3;

h_orb=400*10^3;
i_orb=51.65*pi/180;
i_orb=0;
i_ORBITA=0;
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
%[R_0t,V_0t]=OrbPar2rv(a_orb,e_orb,i_orb,OM_orb,om_orb,theta0_orb,Mu)
%[R_0t,V_0t]=OrbPar2RV_MY(a_orb,e_orb,i_ORBITA,OM_orb,om_orb,theta0_orb,Mu)

R_0target=[a_orb;0;0]
V_0target=[0;omega*a_orb;0]

w_v=cross(R_0target,V_0target)/(norm(R_0target)^2);
setGlobalw_v(w_v);

ROT_plane0=[cos(om_orb)*cos(OM_orb)-sin(om_orb)*cos(i_orb)*sin(OM_orb) -sin(om_orb)*cos(OM_orb)-cos(om_orb)*cos(i_orb)*sin(OM_orb) sin(i_orb)*sin(OM_orb)
            cos(om_orb)*sin(OM_orb)+sin(om_orb)*cos(i_orb)*cos(OM_orb) -sin(om_orb)*sin(OM_orb)+cos(om_orb)*cos(i_orb)*cos(OM_orb) -sin(i_orb)*cos(OM_orb)
              sin(om_orb)*sin(i_orb)                                                 cos(om_orb)*sin(i_orb)                                 cos(i_orb)]  ;
ROT_plane=eye(3);

          
setGlobalROT_plane(ROT_plane);



%% 2D Trajectory Load:

% all trajectories coming from file: Obstacle_Avoidance_test1.m


%traj1
load('t_traj1')
load('xxx_traj1')
load('y_traj1')
load('vx_traj1')
load('vy_traj1')

[SP,~]=size(t_traj1);

Set_Point=[t_traj1,xxx_traj1,y_traj1,zeros(SP,1),vx_traj1,vy_traj1,zeros(SP,1)];

%traj2
load('Set_Point_traj_2'); %PASS RIGHT x_obs~[100;100;0]
Set_Point=Set_Point_traj_2;

%traj3
load('Set_Point_traj_3'); % PASS LEFT:  x_obs~[ 0;100;0] 
Set_Point=Set_Point_traj_3;

%% 3D Trajectory Load:


% all trajectories coming from : trajectory_design_3d_test_0.m

%traj4_3d !! achieved imposing profile velocity on actuators limit
                                      %===> to be tested with horizon=0
load('Set_Point_traj_4_3d'); %     %x_obs=[-47    67     0] m_cubesat=5 f_sat=2*10^-3
Set_Point=Set_Point_traj_4_3d;
                           


%traj_2_3d
load('Set_Point_traj_2_3d'); %     X_obs=[-47;67 ;0]
Set_Point=Set_Point_traj_2_3d;





%traj6_3d
 load('Set_Point_traj_6_3d'); %    x_obs=[-47    53    20]
Set_Point=Set_Point_traj_6_3d;

%traj5_3d
 load('Set_Point_traj_5_3d'); %     x_obs=[-47 67  20]
Set_Point=Set_Point_traj_5_3d;


%traj_3_3d
load('Set_Point_traj_3_3d'); %     X_obs=[-47;67 ;0]
Set_Point=Set_Point_traj_3_3d;

%traj1_3d
 load('Set_Point_traj_1_3d'); %     x_obs=[53 67 10]
Set_Point=Set_Point_traj_1_3d;

[size_Set_Point,~]=size(Set_Point)

%% 3D trajectory load from discrete laplace equation

% CI_Iz obstacle
load('Set_Point_traj_2')
Set_Point=Set_Point_traj_2; %x_obs=[0;-10;15]

load('Set_Point_traj_2_Dt_impulse_200s')
Set_Point=Set_Point_traj_2;
L_box=60;

X_max=L_box;
Y_max=L_box;
Z_max=L_box;

setGlobalL_box_x(X_max)
setGlobalL_box_y(Y_max)
setGlobalL_box_z(Z_max)

step_grid=1;

R1_boundedx=-X_max:step_grid:X_max;
R1_boundedy=-Y_max:step_grid:Y_max;
R1_boundedz=-Z_max:step_grid:Z_max;

[~,N]=size(R1_boundedx);
[~,M]=size(R1_boundedy);
[~,P]=size(R1_boundedz);

Phi0_3d_plot_plus=zeros(M,N);
Phi0_3d_plot_plus_Iz=zeros(M,N);
Phi0_3d_plot_minus=zeros(M,N);

height_ISS_CI_Iz=10;
lenght_ISS_CI_Iz=25;
width_ISS_CI_Iz=15;

for i=2:N-1
     for j=2:M-1
         for k=2:P-1;
             
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        z=R1_boundedz(k);
        
        
        
           if z <= 0 && z >= - height_ISS_CI_Iz/2 || z > 0 && z<=height_ISS_CI_Iz/2
            
            if y>=0 
                
                  if abs(x)<=lenght_ISS_CI_Iz && abs(y) <= height_ISS_CI_Iz/2
                      
                  Phi0_3d_plot_plus(i,j)=height_ISS_CI_Iz/2;
                  Phi0_3d_plot_minus(i,j)=-height_ISS_CI_Iz/2;
                
                
                
                  elseif abs(x) <=lenght_ISS_CI_Iz   && abs(x) >=2*height_ISS_CI_Iz && abs(y)<=2*width_ISS_CI_Iz
                      
                 Phi0_3d_plot_plus(i,j)=height_ISS_CI_Iz/2;
                 Phi0_3d_plot_minus(i,j)=-height_ISS_CI_Iz/2;
                     
                  
                      
                  end
                
            elseif y<0 
                
                
                if abs(x)<=height_ISS_CI_Iz/2 && abs(y) <= 1*width_ISS_CI_Iz
                    
                Phi0_3d_plot_plus(i,j)=height_ISS_CI_Iz/2;
                Phi0_3d_plot_minus(i,j)=-height_ISS_CI_Iz/2;
                
                
                
                
                end
                
               
            end
            
            
            
            
           elseif z > height_ISS_CI_Iz/2 && z <=2*width_ISS_CI_Iz
               
               
               if abs(x)< height_ISS_CI_Iz/2 && abs(y) < height_ISS_CI_Iz/2
               
                Phi0_3d_plot_plus_Iz(i,j)=2*width_ISS_CI_Iz;
               
               end
               
           
             
        
           end
           
           
           
  
  
     
         end    
     end
 end



%% observation Point [COHERENT WITH TRAJECTORY]& Obstacle (ISS ellipsoid)
x_obs=0;
y_obs=-10;
z_obs=15;




X_obs=[x_obs;y_obs;z_obs];
setGlobalx_obs(X_obs);


lenght_ISS=100; %--> along y
width_ISS=50; %--> alng z
height_ISS=30;


[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,width_ISS,height_ISS);






%% Free Drift trajectory


% !!!! CRUCIAL TEST: by  givinig initial conditions such that x_obs_2d is
% reached no trajectory corrections are made. ====> OK !!! **

Rho_0_3d=[Set_Point(1,2);Set_Point(1,3);Set_Point(1,4)];

tau=1/4*T_orb;
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

%% Chaser Initial Conditions: LVLH-->IN
Rho_0=Set_Point(1,2:4)';
Ni_0=Set_Point(1,5:7)';

R_0chaser=R_0target+ROT_plane*Rho_0; % assuming: theta0_orb=0
V_0chaser=V_0target+ROT_plane*(Ni_0+cross([0;0;w],Rho_0));% assuming: theta0_orb=0


%% actuators & controllers design

%actuator & cubesat.
m_cubesat=5;
f_satur=2*1e-3*m_cubesat;% ref1: a_sat=4*1e-6N && m_cubesat=3Kg 
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


 rr_u=0.9999;
 rr_z=0.0001;
   
Q = Cz'*(rr_z*Wz)*Cz;


R = rr_u*Wu;


K = lqr(A,B,Q,R);

K_lqr=-K

setGlobalk_lqr(K_lqr);
%% ODE [ Integration in IN r.f.]
setGlobalSet_Point(Set_Point);
i_0Sp=1;
setGlobali_Sp(i_0Sp);

Horizon=1; % 0--> static orizon  1--> trajectory tracking 
setGlobalHorizon(Horizon)

if Horizon==1;
t_end=Set_Point(end,1)+1000;%-> trajectory tracking
  
else
t_end=tau+10*T_orb;%-->static horizon 
time_before=0;
setGlobaltime_before(time_before)
end

X_0target=[R_0target;V_0target];
  X_0chaser=[R_0chaser;V_0chaser];
   X_0tar_cha=[X_0target;X_0chaser];
   
options = odeset('RelTol',1e-13,'AbsTol', 1e-6);
[t_LQR_FT,X_LQR_FT]= ode113('r2BP_t_c_LQR_sat_SP_new_WorK_Folder',t_end,X_0tar_cha,options);

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

test_3d=3;%0 for a 2D-test on  x-y ISS ellipsoid's contour
          %1 for a 3D-test on ISS ellipsoid
          %3 for a 3d-test for discete  laplace equation solved problem
          

if test_3d==0
    
figure()
plot3(Rrel_LVLH_LQR_FT(:,1),Rrel_LVLH_LQR_FT(:,2),Rrel_LVLH_LQR_FT(:,3),'--b',Rrel_LVLH_LQR_FT(1,1),Rrel_LVLH_LQR_FT(1,2),Rrel_LVLH_LQR_FT(1,3),'og',Rrel_LVLH_LQR_FT(end,1),Rrel_LVLH_LQR_FT(end,2),Rrel_LVLH_LQR_FT(end,3),'or',x_obs,y_obs,z_obs,'*g',Set_Point(:,2),Set_Point(:,3),Set_Point(:,4),'k')
title('controlloed LQR along reference trajectory in LVLH ref ')
hold on
surf(x_ISS, y_ISS,z_ISS,'facecolor','g')
hold on
grid on

figure()
plot3(Set_Point(:,2),Set_Point(:,3),Set_Point(:,4),'k')
title('refernce trajectory in LVLH ref ')
hold on
grid on

    
    
elseif test_3d==1

figure()
plot3(Rrel_LVLH_LQR_FT(:,1),Rrel_LVLH_LQR_FT(:,2),Rrel_LVLH_LQR_FT(:,3),'--b',Rrel_LVLH_LQR_FT(1,1),Rrel_LVLH_LQR_FT(1,2),Rrel_LVLH_LQR_FT(1,3),'og',Rrel_LVLH_LQR_FT(end,1),Rrel_LVLH_LQR_FT(end,2),Rrel_LVLH_LQR_FT(end,3),'or',x_obs,y_obs,z_obs,'*g',Set_Point(:,2),Set_Point(:,3),Set_Point(:,4),'k')
title('controlloed LQR along reference trajectory in LVLH ref ')
hold on
surf(x_ISS, y_ISS,z_ISS,'facecolor','g')
axis equal
hold on
grid on

figure()
plot3(Set_Point(:,2),Set_Point(:,3),Set_Point(:,4),'k')
title('refernce trajectory in LVLH ref ')
hold on
grid on

elseif test_3d==3

figure()
plot3(Rrel_LVLH_LQR_FT(:,1),Rrel_LVLH_LQR_FT(:,2),Rrel_LVLH_LQR_FT(:,3),'--b',Rrel_LVLH_LQR_FT(1,1),Rrel_LVLH_LQR_FT(1,2),Rrel_LVLH_LQR_FT(1,3),'og',Rrel_LVLH_LQR_FT(end,1),Rrel_LVLH_LQR_FT(end,2),Rrel_LVLH_LQR_FT(end,3),'or',x_obs,y_obs,z_obs,'*g',Set_Point(:,2),Set_Point(:,3),Set_Point(:,4),'k')
hold on
plot3(Set_Point(:,2),Set_Point(:,3),Set_Point(:,4),'k')
hold on
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus')
colormap summer
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_minus')
colormap summer
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus_Iz',30)
colormap autumn
title('controlloed LQR along reference trajectory in LVLH ref ')
hold on
axis equal
hold on
grid on





end


