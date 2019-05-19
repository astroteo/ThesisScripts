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
x_obs=20
y_obs=-40+3
z_obs=0



lenght_ISS=50; %--> along x
width_ISS=25; %--> along y
height_ISS=15;%--> along z


%[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,width_ISS,height_ISS);

%% Initial Conditions & Integration time 

% !!!! CRUCIAL TEST: by  givinig initial conditions such that x_obs_2d is
% reached no trajectory corrections are made. ====> OK !!! **
Rho_0_2d=[20;30];
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

%X_0_2d=[Rho_0_2d;Ni_0_2d];
X_0_3d=[Rho_0_3d;Ni_0_3d];

%% Free Drift !! UNforced Problem!!

t_end=tau;

options = odeset('RelTol',1e-13,'AbsTol', 1e-12);
[t_CW_FD,x_CW_FD]= ode113('Chol_Wilt_Hill',t_end,X_0_3d,options);%(plot at the end of the script)

[CW_FD,~]=size(t_CW_FD);
zebra_FD=zeros(CW_FD,1);

%% 2d Discrete Potential C+I shaped object 
     % usefull to test BUG also
     

L_max=norm(Rho_0_3d); 

L_box=200;

setGlobalL_box_x(L_box)
setGlobalL_box_y(L_box)

X_max=L_box;
Y_max=L_box;

step_grid=1;

R1_boundedx=-X_max:step_grid:X_max;
R1_boundedy=-Y_max:step_grid:Y_max;

[~,N]=size(R1_boundedx);
[~,M]=size(R1_boundedy);

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



Phi0_2d=zeros(N,M);
Phi_BOX=zeros(N,M);

for i=1:N
    for j=1:M
        
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        
       
            
            if y>=0
                
                  if abs(x)<=lenght_ISS && abs(y) <= height_ISS
                
                Phi0_2d(i,j)=1;
                
                  elseif abs(x) <=lenght_ISS   && abs(x) >=3/2*width_ISS && abs(y)<=3/2*width_ISS
                     
                  Phi0_2d(i,j)=1;
                      
                  end
                
            elseif y<0
                
                
                if abs(x)<=height_ISS/2 && abs(y) <= 2*width_ISS
                
                Phi0_2d(i,j)=1;
                end
                
               
            end
            
            if i==N || i==1 || j==M || j==1
           
           Phi0_2d(i,j)=1;
           Phi_BOX(i,j)=0.5;
            end
            
            
            if i==i_obs && j==j_obs
                
            Phi0_2d(i,j)=-1;% only to visualize the goal point to be set to zero in the Laplace solution
            end
            
            
            
     end
end
        
        
        
figure()
hold on
mesh(R1_boundedx,R1_boundedy,Phi0_2d')
title('bounduary conditions')
xlabel('x')
ylabel('y')
grid on  

% Phi_2d=1*ones(M,N)+Phi_BOX;
% Phi_2d_new=1*ones(M,N)+Phi_BOX;


Phi_2d=1*ones(M,N);
Phi_2d_new=1*ones(M,N);
epsilon=1;
count=0;

while epsilon >= 0.1 && count <=10000

 for i=2:N-1
     for j=2:M-1
         
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        
       
            
            if y>=0
                
                  if abs(x)<=lenght_ISS && abs(y) <= height_ISS
                
                  Phi_2d(i,j)=1;
                
                  elseif abs(x) <=lenght_ISS   && abs(x) >=3/2*width_ISS && abs(y)<=3/2*width_ISS
                     
                  Phi_2d(i,j)=1;
                      
                  end
                
            elseif y<0
                
                
                if abs(x)<=height_ISS/2 && abs(y) <= 2*width_ISS
                
                Phi_2d(i,j)=1;
                end
                
               
            end
            
            
            
            
            if i==i_obs && j==j_obs
                
            Phi_2d(i,j)=0;
            end
        
        
           
  Phi_2d_new(i,j)=1/4*(Phi_2d(i-1,j)+Phi_2d(i+1,j)+Phi_2d(i,j-1)+Phi_2d(i,j+1));
     
         
     end
 end
 
 count=count+1;

epsilon=max(max(abs(Phi_2d_new-Phi_2d)));


Phi_2d=Phi_2d_new;

end

figure()
hold on
mesh(R1_boundedx,R1_boundedy,Phi_2d')
title('potential')
xlabel('x')
ylabel('y')
grid on  

figure()
contour(R1_boundedx,R1_boundedy,Phi_2d')
hold on
title('total potential contour')
plot3(x_obs,y_obs,0,'*b')
xlabel('x')
ylabel('y')
grid on
hold on


%% ALP method 

%parameters and global variables

setGlobalstep_grid(step_grid)
setGlobalPhi_2d(Phi_2d)


% alpha IsSue
alpha_v=norm(Rho_0_3d-X_obs_3d)/(1*t_end)
   % alpha_v=norm(Ni_0_2d)*1e-4;
 setGlobalalpha_v(alpha_v );

%step_time to reduce control action
step_time=5;
setGlobalstep_time(step_time)
flag_time=0;
setGlobalflag_time(flag_time )
time_before=0;
setGlobaltime_before(time_before)
time_last_fire=0;
setGlobaltime_last_fire(time_last_fire)

dt=0.1;
t_span=[0:dt:8*tau]';
X_0_3d=[Rho_0_3d;Ni_0_3d];

[Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done,i_obs]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_2d_step_time(dt,t_span,X_0_3d,X_obs_3d);

figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_an(:,1),Rho_an(:,2),Rho_an(:,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(end,1),Rho_an(end,2),Rho_an(end,3),'or',x_obs,y_obs,0,'*b')
title('collision Avoidance with discrete gradient & analytical trajectory propagation')
grid on 

figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_an(1:i_obs,1),Rho_an(1:i_obs,2),Rho_an(1:i_obs,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(i_obs,1),Rho_an(i_obs,2),Rho_an(i_obs,3),'or',x_obs,y_obs,0,'*b')
title('collision Avoidance with discrete gradient & analytical trajectory propagation')
grid on 

figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(x_CW_FD(:,1),x_CW_FD(:,2),zebra_FD,'--b',x_CW_FD(1,1),x_CW_FD(1,2),zebra_FD,'og',x_CW_FD(end,1),x_CW_FD(end,2),zebra_FD,'or',x_obs,y_obs,0,'*b')
title('free drift')
grid on 


%results

DV_tot=sum(norm(DV_v));
Goal_error=norm(Rho_an(end,:)'-X_obs_3d)
Total_Impulse=count_impulse_done

DV_tot_obs=sum(norm(DV_v(1:i_obs,:)))
TOF_obs=i_obs*dt
