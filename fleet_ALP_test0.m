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
x_obs=100+3
y_obs=100-3
z_obs=0



lenght_ISS=100; %--> along y
width_ISS=50; %--> alng z
height_ISS=30;


[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,height_ISS,width_ISS);

%% Initial Conditions & Integration time 

% !!!! CRUCIAL TEST: by  givinig initial conditions such that x_obs_2d is
% reached no trajectory corrections are made. ====> OK !!! **
Rho_0_2d=[300;-300];
Rho_0_3d=[Rho_0_2d(1);Rho_0_2d(2);0];

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
Ni_0_2d_CT=[Ni_0_3d(1);Ni_0_3d(2)]; %** CRUCIAL TEST
Ni_0_2d_RAND=[-rand(1)*1e-1; -rand(1)*1e-1]; %% GENERIC STUFF
Ni_0_2d=Ni_0_2d_CT-Ni_0_2d_RAND;

X_0_2d=[Rho_0_2d;Ni_0_2d];


%% Free Drift !! UNforced Problem!!

t_end=tau;
X_0_3d=[Rho_0_3d;Ni_0_3d];

options = odeset('RelTol',1e-13,'AbsTol', 1e-12);
[t_CW_FD,x_CW_FD]= ode113('Chol_Wilt_Hill',t_end,X_0_3d,options);%<--[]

[CW_FD,~]=size(t_CW_FD);
zebra_FD=zeros(CW_FD,1);

figure()
surf(x_ISS, y_ISS,zeros(size(x_ISS)),'facecolor','g')
hold on
plot3(x_CW_FD(:,1),x_CW_FD(:,2),zebra_FD,'--b',x_CW_FD(1,1),x_CW_FD(1,2),zebra_FD,'og',x_CW_FD(end,1),x_CW_FD(end,2),zebra_FD,'or',x_obs,y_obs,0,'*b')
grid on 

%% Grid Properties 

L_max=1600;  

step_grid=20;

R1_bounded=-L_max:step_grid:(L_max-step_grid);
    
R1_boundedx=R1_bounded;
R1_boundedy=R1_bounded;

[~,N]=size(R1_boundedx);
[~,M]=size(R1_boundedy);



%% Position Potential PARAMETERS

% ??? m1=lenght_ISS; ???
m1=1;

% ??? m2=width_ISS; ???
m2=1;


M_2d=[m1, 0;
      0, m2];

%% Position Potential
X_obs_2d=[x_obs;y_obs];%-> observation point given in LVLH r.f.

phi_pos_2d=zeros(N,M);

for i=1:N
    
    for j=1:M
        
       
            
            x=R1_boundedx(1,i)/L_max;
            y=R1_boundedy(1,j)/L_max;
           
            X_2d= [x;y];
            
           phi_pos_2d(i,j)=1/2*(X_2d-X_obs_2d/L_max)'*(M_2d*(X_2d-X_obs_2d/L_max));
            
            
   
    end
end
%% Spherical Harmonic 2D Potential PARAMETES

X_coll_2d=[-500;-500];

lambda_1_paper_Rihch=1/2*(Rho_0_2d./L_max-X_coll_2d./L_max)'*(M_2d*(Rho_0_2d./L_max-X_coll_2d./L_max));
lambda_1_paper=1e6;
lambda_1_logic=1/2*([L_max;L_max]./L_max-X_coll_2d./L_max)'*(M_2d*([L_max;L_max]/L_max-X_coll_2d./L_max));
lambda_1_MY=max(max(phi_pos_2d));
lambda_1_illogic=1e6;

lambda_1=1*lambda_1_MY;

safety_radius=-10;

lambda_2_MY=((lenght_ISS+safety_radius)/2)^2/(0.001*L_max);
lambda_2_paper=1e-4;

lambda_2=lambda_2_MY;

lambda_v=[lambda_1;lambda_2];

n1=(lenght_ISS+safety_radius)/lenght_ISS;
n2=(width_ISS+safety_radius)/lenght_ISS;

N_2d=[n1,0;
      0, n2]; % spherical-shape of potential around obstacle.
  
%N_2d=eye(2);

%% Spherical Armonic Potential [with 1 swarm]

phi_harm_2d=zeros(N,M);

lambda_1_f=lambda_1;

lambda_2_f_MY=1e+3;
lambda_2_f=lambda_2_f_MY;

N_2d_f=eye(2)


for i=1:N
    
    for j=1:M
        
       
            
            x=R1_boundedx(1,i)/L_max;%--> necessity to adimensionalize
            y=R1_boundedy(1,j)/L_max;
           
            X_2d= [x;y];
            
           
           phi_harm_2d(i,j)=lambda_1*exp(-lambda_2*((X_2d-X_coll_2d/L_max))'*N_2d_f*((X_2d-X_coll_2d/L_max)))+lambda_1_f*exp(-lambda_2_f*X_2d'*N_2d*X_2d);
           
           
    end
end


%% Total Potential & Plots

phi_2d=phi_harm_2d+phi_pos_2d;

R1_boundedx=R1_bounded;
R1_boundedy=R1_bounded;

Disp_Potentials=1;
if Disp_Potentials==1

figure()
mesh(R1_boundedx,R1_boundedy,phi_pos_2d)
hold on
title('"position" potential')
grid on  

figure()
contour(R1_boundedx,R1_boundedy,phi_pos_2d)
hold on
plot3(x_obs,y_obs,0,'*b',Rho_0_2d(1),Rho_0_2d(2),0,'ob')
title('Phi_{pos} contour') 
hold on

figure()
mesh(R1_boundedx,R1_boundedy,phi_harm_2d)
title('harmonic potential')
grid on  

figure()
contour(R1_boundedx,R1_boundedy,phi_harm_2d)
hold on
plot3(x_obs,y_obs,0,'*b',Rho_0_2d(1),Rho_0_2d(2),0,'ob')
title('Phi_{harm} contour') 
hold on



figure()
contour(R1_boundedx,R1_boundedy,phi_2d)
hold on
title('total potential contour')
plot3(x_obs,y_obs,0,'*b')
hold on

figure()
hold on
mesh(R1_boundedx,R1_boundedy,phi_2d)
title('total potential')
grid on  

% quiver(R1_bounded,R1_bounded,px,py)
% hold off
end

%% CHECK: analytical gradient [INSIDE ODE]

x_t=100;% !! MUST BE  A GRID POINT TO DO THE CHECK!!
y_t=100;

% x_t=X_obs_2d(1);% !! MUST BE MINIMUM GADIENT!!
% y_t=X_obs_2d(2);

total_stepx=N;
total_stepy=M;

step_gridx=step_grid;
step_gridy=step_grid;


if x_t >=0
i_bottom=total_stepx/2+floor(x_t/step_gridx);
i_up=i_bottom+1;
else     
i_up=ceil((total_stepx/2*step_grid-abs(x_t))/step_gridx);
i_bottom=i_up-1;
end


if y_t >=0
j_bottom=total_stepy/2+floor(y_t/step_gridy);
j_up=j_bottom+1;
else     
j_up=ceil((total_stepy/2*step_gridy-abs(y_t))/step_gridy);
j_bottom=j_up-1;
end

[px_harm,py_harm] = gradient(phi_harm_2d,step_grid,step_grid);
[px_pos,py_pos] = gradient(phi_pos_2d,step_grid,step_grid);

px=px_harm(i_bottom,j_bottom)+px_pos(i_bottom,j_bottom);% up/bottom is arbitrary since i'm choosing a knot
py=py_harm(i_bottom,j_bottom)+py_pos(i_bottom,j_bottom);


p_dir_MAT=[px/norm([px;py]);py/norm([px;py])]


%analylitically computed gradient (inside ODE)
f_harm_x=lambda_1*(2*n1*(lambda_2)*(x_t/L_max)*exp((-lambda_2)*(n1*(x_t/L_max)^2+n2*(y_t/L_max)^2)));
f_harm_y=lambda_1*(2*n2*(lambda_2)*(y_t/L_max)*exp((-lambda_2)*(n1*(x_t/L_max)^2+n2*(y_t/L_max)^2)));

%phi_2d(i,j)=1/2*(X_2d/L_max-X_obs_2d/L_max)'*M_2d*(X_2d/L_max-X_obs_2d/L_max);
f_pos_x=m1/L_max*(x_t/L_max-x_obs/L_max);
f_pos_y=m2/L_max*(y_t/L_max-y_obs/L_max);

%d/dx(f(x,y)+g(x,y))=df(x,y)/dx+dg(x,y)/dx 
%d/dy(f(x,y)+g(x,y))=df(x,y)/dy+dg(x,y)/dy
fx=f_harm_x+f_pos_x;
fy=f_harm_y+f_pos_y;

f_dir=[fx/norm([fx;fy]);fy/norm([fx;fy])]




%% GNC by Laplace Aritificial Potential [ANALYCAL GRADIENT]

%parameters to compute analytical gradient
setGloball_max(L_max)
setGlobalM_2d(M_2d)
setGlobalN_2d(N_2d)
setGloballambda_v(lambda_v)
setGlobalx_obs_2d(X_obs_2d )


% alpha IsSue
alpha_v=norm(Rho_0_2d-X_obs_2d)/t_end
% alpha_v=norm(Ni_0_2d)*1e-4;
setGlobalalpha_v(alpha_v );

Perform_Integration=0;
if Perform_Integration==1
%integration time
t_end_Lap=t_end+1.5*T_orb

%integration
options = odeset('RelTol',1e-7,'AbsTol', 1e-1);
[t_CW,x_CW]= ode113('Chol_Wilt_Hill_Laplace_2d_analytic',t_end_Lap,X_0_2d,options);%<--[ISSUE: n° of correction strictly dependant on ODE's precision]


%plot
[CW,~]=size(t_CW);
zebra=zeros(CW,1);

figure()
surf(x_ISS, y_ISS,zeros(size(x_ISS)),'facecolor','g')
hold on
plot3(x_CW(:,1),x_CW(:,2),zebra,'--b',x_CW(1,1),x_CW(1,2),zebra,'og',x_CW(end,1),x_CW(end,2),zebra,'or',x_obs,y_obs,0,'*b')
title('collision Avoidance with analytical gradient')
grid on 

end


%% GNC by Laplace Artificial potential CONTROL EFFORT AT DISCRETE TIMES

%parameters to compute analytical gradient
setGloball_max(L_max)
setGlobalM_2d(M_2d)
N_2d=eye(2);
setGlobalN_2d(N_2d)
setGloballambda_v(lambda_v)
setGlobalx_obs_2d(X_obs_2d )

%step_time to reduce control action
step_time=100;
setGlobalstep_time(step_time)
flag_time=0;
setGlobalflag_time(flag_time )
time_before=0;
setGlobaltime_before(time_before)
time_last_fire=0;
setGlobaltime_last_fire(time_last_fire)



% alpha IsSue
alpha_v=norm(Rho_0_2d-X_obs_2d)/(1*t_end)
% alpha_v=norm(Ni_0_2d)*1e-4;
setGlobalalpha_v(alpha_v );

% %integration time
t_end_Lap_st=t_end+40*T_orb

t_span_Lap_st=linspace(0,t_end_Lap_st,10000);

%  options = odeset('RelTol',1e-10,'AbsTol', 1e-2);
%  [t_CW_st,x_CW_st]= ode113('Chol_Wilt_Hill_Laplace_2d_analytic_step_time',t_end_Lap_st,X_0_2d,options);%<--[FIXED STEP SOLVER]
% 
% %plot
% [CW_st,~]=size(t_CW_st);
% zebra_st=zeros(CW_st,1);
% 
% 
% figure()
% surf(x_ISS, y_ISS,zeros(size(x_ISS)),'facecolor','g')
% hold on
% plot3(x_CW_st(:,1),x_CW_st(:,2),zebra_st,'--b',x_CW_st(1,1),x_CW_st(1,2),zebra_st,'og',x_CW_st(end,1),x_CW_st(end,2),0,'or',x_obs,y_obs,0,'*b')
% title('collision Avoidance with analytical gradient: TIME STEP')
% grid on 
%% Full Analytical Computation [discrete time control action]

dt=10;
t_span=[0:dt:2*tau]';

[Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Full_Analytic_Laplace_2d_step_time(dt,t_span,X_0_3d);

count_impulse_time
count_impulse_done
[size_ODE,~]=size(t_CW_FD)
%[size_ODE_ctrl,~]=size(t_CW)
[size_SPAN,~]=size(t_span)

figure()
surf(x_ISS, y_ISS,zeros(size(x_ISS)),'facecolor','g')
hold on
plot3(Rho_an(:,1),Rho_an(:,2),Rho_an(:,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(end,1),Rho_an(end,2),Rho_an(end,3),'or',x_obs,y_obs,0,'*b')
title('collision Avoidance with analytical gradient & analytical trajectory propagation')
grid on 

Set_Point_traj_2=[t_span,Rho_an,Ni_an]

save('Set_Point_traj_2')