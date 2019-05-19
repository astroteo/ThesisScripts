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
x_obs=10
y_obs=-20
z_obs=0



lenght_ISS=25; %--> along x
width_ISS=15; %--> along y
height_ISS=10;%--> along z


%[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,width_ISS,height_ISS);

%% Initial Conditions & Integration time 

% !!!! CRUCIAL TEST: by  givinig initial conditions such that x_obs_2d is
% reached no trajectory corrections are made. ====> OK !!! **

Rho_0_3d=[-20;40;0];

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


%% Laplace equation solution

Phi_3d=1*ones(M,N,P);
Phi_3d_new=1*ones(M,N,P);
epsilon=1;
count=0;

Phi0_3d_plot_plus=zeros(M,N);
Phi0_3d_plot_minus=zeros(M,N);

displaY_count=0;

while epsilon >= 0.1 && count <=2000

 for i=2:N-1
     for j=2:M-1
         for k=2:P-1;
             
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        z=R1_boundedz(k);
        
           if z <= height_ISS/2
            
            if y>=0 
                
                  if abs(x)<=lenght_ISS && abs(y) <= height_ISS/2
                      
                  Phi0_3d_plot_plus(i,j)=height_ISS/2;
                  Phi0_3d_plot_minus(i,j)=-height_ISS/2;
                
                  Phi_3d(i,j,k)=1;
                
                  elseif abs(x) <=lenght_ISS   && abs(x) >=2*height_ISS && abs(y)<=2*width_ISS
                      
                 Phi0_3d_plot_plus(i,j)=height_ISS/2;
                 Phi0_3d_plot_minus(i,j)=-height_ISS/2;
                     
                  Phi_3d(i,j,k)=1;
                      
                  end
                
            elseif y<0 
                
                
                if abs(x)<=height_ISS/2 && abs(y) <= 1*width_ISS
                    
                Phi0_3d_plot_plus(i,j)=height_ISS/2;
                Phi0_3d_plot_minus(i,j)=-height_ISS/2;
                
                Phi_3d(i,j,k)=1;
                
                
                end
                
               
            end
            
            
            if i==i_obs && j==j_obs && k==k_obs
            
            Phi_3d(i,j,k)=0;
            end
            
           
            
               
        
           end
           
           
           
  Phi_3d_new(i,j,k)=1/6*(Phi_3d(i-1,j,k)+Phi_3d(i+1,j,k)+Phi_3d(i,j-1,k)+Phi_3d(i,j+1,k)+Phi_3d(i,j,k-1)+Phi_3d(i,j,k+1));
  
  
     
         end    
     end
 end
 
 
count=count+1;

if count-displaY_count >= 100
    
    displaY_count=displaY_count+100
    epsilon
end

epsilon=max(max(max(abs(Phi_3d_new-Phi_3d))));


Phi_3d=Phi_3d_new;

end

% figure()
% hold on
% contour3(R1_boundedx,R1_boundedy,R1_boundedz,Phi_3d)
% plot3(x_obs,y_obs,0,'*b')
% title('potential')
% xlabel('x')
% ylabel('y')
% grid on  

% figure()
% contour(R1_boundedx,R1_boundedy,Phi_2d')
% hold on
% title('total potential contour')
% plot3(x_obs,y_obs,0,'*b')
% xlabel('x')
% ylabel('y')
% grid on
% hold on

save('Phi_3d')

figure()
hold on
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus')
colormap summer
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_minus')
colormap summer
title('bounduary conditions')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid on  




Phi_2d_xy=Phi_3d(:,:,(P-1)/2);
figure()
mesh(R1_boundedx,R1_boundedy,Phi_2d_xy')
colormap summer
hold on
plot(x_obs,y_obs,'*b')
title('Laplace potential x-y')
grid on 


figure()
contour(R1_boundedx,R1_boundedy,Phi_2d_xy')
colormap summer
hold on
plot(x_obs,y_obs,'*b')
title('Laplace potential contour x-y')
grid on 


Phi_2d_xz=zeros(N,P);
 for I=1:N
for K=1:P
c=Phi_3d(I,j_obs,K);
Phi_2d_xz(I,K)=c;
end
end

figure()
mesh(R1_boundedx,R1_boundedz,Phi_2d_xz')
colormap summer
hold on
plot(x_obs,z_obs,'*b')
title('Laplace potential x-z')
grid on 


figure()
contour(R1_boundedx,R1_boundedz,Phi_2d_xz')
colormap summer
hold on
plot(y_obs,z_obs,'*b')
title('Laplace potential contour x-z')
grid on 



%% ALP method 

%parameters and global variables

setGlobalstep_grid(step_grid)
setGlobalPhi_3d(Phi_3d)

% alpha IsSue
alpha_v=norm(Rho_0_3d-X_obs_3d)/(1*t_end)
   % alpha_v=norm(Ni_0_2d)*1e-4;
 setGlobalalpha_v(alpha_v );

%step_time to reduce control action
step_time=10;
setGlobalstep_time(step_time)
flag_time=0;
setGlobalflag_time(flag_time )
time_before=0;
setGlobaltime_before(time_before)
time_last_fire=0;
setGlobaltime_last_fire(time_last_fire)
setGloaltau_or(tau)

dt=0.1;
t_span=[0:dt:2*tau]';
X_0_3d=[Rho_0_3d;Ni_0_3d];

% [Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done,i_obs]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_2d_step_time(dt,t_span,X_0_3d,X_obs_3d);
  [Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_3d_step_time(dt,t_span,X_0_3d);



figure()
hold on
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus')
colormap summer
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_minus')
colormap summer
plot3(Rho_an(:,1),Rho_an(:,2),Rho_an(:,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(end,1),Rho_an(end,2),Rho_an(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('collision Avoidance with discrete gradient & analytical trajectory propagation')
axis equal
grid on 

figure()
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus')
colormap summer
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_minus')
colormap summer
hold on
plot3(x_CW_FD(:,1),x_CW_FD(:,2),x_CW_FD(:,3),'--b',x_CW_FD(1,1),x_CW_FD(1,2),x_CW_FD(1,3),'og',x_CW_FD(end,1),x_CW_FD(end,2),x_CW_FD(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('free drift')
axis equal
grid on 


%results

DV_tot=sum(norm(DV_v));
Goal_error=norm(Rho_an(end,:)'-X_obs_3d)
Total_Impulse=count_impulse_done

DV_tot_obs=sum(norm(DV_v(1:i_obs,:)))
TOF_obs=i_obs*dt


%% crucial test for gradient

Phi_3d_CT=zeros(M,N,P);

 for i=1:N
     for j=1:M
         for k=2:P;
             
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        z=R1_boundedz(k);
        
          
           
           r=[x;y;z];
           
           
  Phi_3d_CT(i,j,k)=1/2*((r-X_obs_3d)'*(r-X_obs_3d));
  
  
     
         end    
     end
 end
 
setGlobalstep_grid(step_grid)
setGlobalPhi_3d(Phi_3d_CT)

[Rho_an_CT,Ni_an_CT,DV_v_CT,count_impulse_time_CT,count_impulse_done_CT]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_3d_step_time(dt,t_span,X_0_3d);
 
figure()
plot3(Rho_an_CT(:,1),Rho_an_CT(:,2),Rho_an_CT(:,3),'--b',Rho_an_CT(1,1),Rho_an_CT(1,2),Rho_an_CT(1,3),'og',Rho_an_CT(end,1),Rho_an_CT(end,2),Rho_an_CT(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('crucial test for gradient')
axis equal
hold on
grid on 


Phi_2d_xy_CT=Phi_3d_CT(:,:,(P-1)/2);
figure()
mesh(R1_boundedx,R1_boundedy,Phi_2d_xy_CT')
colormap summer
hold on
plot(x_obs,y_obs,'*b')
title('Harmonic (goal) potential x-y')
grid on 


figure()
contour(R1_boundedx,R1_boundedy,Phi_2d_xy_CT')
colormap summer
hold on
plot(x_obs,y_obs,'*b')
title('Harminc (goal) potential contour x-y')
grid on 

