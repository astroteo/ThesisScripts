clc
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
x_obs=30
y_obs=50
z_obs=0



lenght_ISS=20; %--> along x
width_ISS=-5; %--> along y
height_ISS=15;%--> along z


%[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,width_ISS,height_ISS);

%% Initial Conditions & Integration time 

% !!!! CRUCIAL TEST: by  givinig initial conditions such that x_obs_2d is
% reached no trajectory corrections are made. ====> OK !!! **
Rho_0_2d=[-50;-50;0];
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


Static=1; %1--> start with null velocity w.r.t. to target
          %0 --> start with C.I. for reaching the goal in TOF=tau= T_orb/4

if Static==0;
X_0_3d=[Rho_0_3d;Ni_0_3d];
else
X_0_3d=[Rho_0_3d;zeros(3,1)];
end

%% Free Drift !! UNforced Problem!!

t_end=tau;

options = odeset('RelTol',1e-13,'AbsTol', 1e-12);
[t_CW_FD,x_CW_FD]= ode113('Chol_Wilt_Hill',t_end,X_0_3d,options);%(plot at the end of the script)

[CW_FD,~]=size(t_CW_FD);
zebra_FD=zeros(CW_FD,1);

%% 2d Discrete Potential C+I shaped object 
     
lenght_ISS_CI=60; %--> along x !! reduced dymension for have faster computations !!
width_ISS_CI=20; %--> along y
height_ISS_CI=20;


L_box=200;

setGlobalL_box_x(L_box)
setGlobalL_box_y(L_box)

X_max=L_box;
Y_max=L_box;

step_grid=5;

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
                
                  if abs(x)<=lenght_ISS_CI && abs(y) <= height_ISS_CI
                
                Phi0_2d(i,j)=1;
                
                  elseif abs(x) <=lenght_ISS_CI   && abs(x) >=2/3*lenght_ISS_CI && abs(y)<=3*height_ISS_CI
                     
                  Phi0_2d(i,j)=1;
                      
                  end
                
            elseif y<0
                
                
                if abs(x)<=height_ISS_CI && abs(y) <= 3*height_ISS_CI
                
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


%% Laplace equation Solution                                          

Phi_2d=ones(M,N);
Phi_2d_new=ones(M,N);

count=0;
epsilon=1;


displaY_count=0;

tic

while epsilon >= 0.1 && count <=100000

 for i=2:N-1
     for j=2:M-1
         
        
           if Phi0_2d(i,j)==1 %<-- smarter way to iterate bounduary conditions
               Phi_2d(i,j)=1;
           end
            
            
            
            
            
            if i==i_obs && j==j_obs
                
            Phi_2d(i,j)=0;
            end
        
        
           
  Phi_2d_new(i,j)=1/4*(Phi_2d(i-1,j)+Phi_2d(i+1,j)+Phi_2d(i,j-1)+Phi_2d(i,j+1));
     
         
     end
 end
 
 count=count+1;

epsilon=max(max(abs(Phi_2d_new-Phi_2d)));

if count-displaY_count >= 1000
    
    displaY_count=displaY_count+1000
    epsilon
end


Phi_2d=Phi_2d_new;

end
               
            
 toc           
%% Countour (Bounduary-generation)
    
 


Phi_2d_contour=zeros(N,M);
n_el_goal=1

for i=1:N
    for j=1:M
        
       
        
        
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        
        
        
%      OBSTACLE  
       
     if y>=0
                
                  if abs(x)==lenght_ISS_CI  && abs(y) <= height_ISS_CI 
                
                      Phi_2d_contour(i,j)=1;
                
                  elseif abs(x)<=lenght_ISS_CI && abs(x) <=2/3*(lenght_ISS_CI) &&  abs(y) == height_ISS_CI
                      
                      Phi_2d_contour(i,j)=1;
                
                  elseif abs(x) ==lenght_ISS_CI && abs(y)<=3*(height_ISS_CI) && abs(y) >= (height_ISS_CI) || abs(x) ==2/3*(lenght_ISS_CI)  && abs(y)<=3*height_ISS_CI && abs(y) >=(height_ISS_CI)
                     
                  Phi_2d_contour(i,j)=1;
                  
                  elseif abs(x) <=(lenght_ISS_CI)  &&  abs(x) >=2/3*(lenght_ISS_CI)  && abs(y)==3*(height_ISS_CI)
                     
                  Phi_2d_contour(i,j)=1;
                  
                  
                  elseif abs(x)==lenght_ISS_CI  && abs(y) == 0
                
                Phi_2d_contour(i,j)=1;
                  end
                
                if abs(x) <=(lenght_ISS_CI) && abs(x)>=(height_ISS_CI) && y ==0
                     
                  Phi_2d_contour(i,j)=1;
                  
                end
                  
                      
             
                
     elseif y<0
                
                
                if abs(x)==(height_ISS_CI) &&  abs(y) <= 3*(height_ISS_CI) 
                
                Phi_2d_contour(i,j)=1;
                end
                
                
                
                
                if abs(x)<=(height_ISS_CI)  &&  abs(y) == 3*(height_ISS_CI)
                
                Phi_2d_contour(i,j)=1;
                end
                
               
    
     end

    %B0X
    
    
                    if i==N-1 ||i==2|| j==M-1 ||j==2
           
                    Phi_2d_contour(i,j)=1;
                   
                    end
                    
                    
                    
     %GOAL
                   
                    if i==i_obs+n_el_goal  && j>=j_obs-n_el_goal   && j<=j_obs+n_el_goal ||  i==i_obs-n_el_goal  && j>=j_obs-n_el_goal   && j<=j_obs+n_el_goal
                  Phi_2d_contour(i,j)=-1;
                    end
                    
                    if j==j_obs+n_el_goal  && i>=i_obs-n_el_goal   && i<=i_obs+n_el_goal|| j==j_obs-n_el_goal  && i>=i_obs-n_el_goal    && i<=i_obs+n_el_goal
                    Phi_2d_contour(i,j)=-1;
                    end
   end
end
    



figure()
 hold on
mesh(R1_boundedx,R1_boundedy,Phi_2d_contour')
title('contour with goal')
xlabel('x')
ylabel('y')
grid on  

%% element generation from contour
I=0;
J=0;
L=0;

I_goal=0;
J_goal=0;
L_goal=0;   % X_el=[x_el,y_el,1/0 obstacle/goal, 1/-1horizontal/vertical]


for i=2:1:N-1
    for j=2:1:M-1
        
        
            
            
        
        
            if (Phi_2d_contour(i+1,j)==1 && Phi_2d_contour(i-1,j)==1 )&& ( Phi_2d_contour(i,j+1)==0 && Phi_2d_contour(i,j-1)==0)
                I=I+1;
                X_el_Hor(I,:)=[R1_boundedx(i),R1_boundedy(j)];
                
                
                L=L+1;
                X_el(L,1:2)=[R1_boundedx(i),R1_boundedy(j)];
                X_el(L,4)=1;
                X_el(L,3)=1;
            end
            
            if Phi_2d_contour(i,j+1)==1 && Phi_2d_contour(i,j-1)==1 &&  (Phi_2d_contour(i+1,j)==0 && Phi_2d_contour(i-1,j)==0)
                J=J+1;
                X_el_Vert(J,:)=[R1_boundedx(i),R1_boundedy(j)];
                
                L=L+1;
                X_el(L,1:2)=[R1_boundedx(i),R1_boundedy(j)];
                X_el(L,3)=1;
                X_el(L,4)=-1;
                
            end
            
        
        
         if Phi_2d_contour(i,j)==-1
             
            if Phi_2d_contour(i+1,j)==-1 && Phi_2d_contour(i-1,j)==-1 && ( Phi_2d_contour(i,j+1)==0 && Phi_2d_contour(i,j-1)==0)
                I_goal=I_goal+1;
                X_el_Hor_goal(I_goal,:)=[R1_boundedx(i),R1_boundedy(j)];
                
                L=L+1;
                X_el(L,1:2)=[R1_boundedx(i),R1_boundedy(j)];
                X_el(L,4)=1;
                X_el(L,3)=0;
                
          
            end
            
            if Phi_2d_contour(i,j+1)==-1 && Phi_2d_contour(i,j-1)==-1 && (Phi_2d_contour(i+1,j)==0 && Phi_2d_contour(i-1,j)==0)
                J_goal=J_goal+1;
                X_el_Vert_goal(J_goal,:)=[R1_boundedx(i),R1_boundedy(j)];
                
                 L=L+1;
                X_el(L,1:2)=[R1_boundedx(i),R1_boundedy(j)];
                X_el(L,4)=-1;
                X_el(L,3)=0;
                
          
            end
            
            
         end
        
         
        
    end
end

figure()
 hold on
contour(R1_boundedx,R1_boundedy,Phi_2d_contour')
hold on
plot(X_el_Hor(:,1),X_el_Hor(:,2),'*b')
hold on
plot(X_el_Hor_goal(:,1),X_el_Hor_goal(:,2),'*b')
hold on
title('Horizontal elements')
xlabel('x')
ylabel('y')
grid on

figure()
 hold on
contour(R1_boundedx,R1_boundedy,Phi_2d_contour')
hold on
plot(X_el_Vert(:,1),X_el_Vert(:,2),'*r')
hold on
plot(X_el_Vert_goal(:,1),X_el_Vert_goal(:,2),'*r')
hold on
title('Vertical elements')

figure()
 hold on
contour(R1_boundedx,R1_boundedy,Phi_2d_contour')
hold on
plot(X_el(:,1),X_el(:,2),'*k')
hold on
title('all elements')

%% G matrix and Sigma vector

%G=zeros(L,L);

l_el=1/2*step_grid

%phi-BC generation


for p=1:1:L
    
    if X_el(p,3)==1
        
        phi_BC(p,1)=1;
        
    elseif X_el(p,3)==0
        
        phi_BC(p,1)=0;
    
    end
end

%G matrix computation

tic

for i=1:1:L
    for j=1:1:L
        
        x_i=X_el(i,1:2);
        x_j=X_el(j,1:2);
        
        
        
            
        if X_el(j,4)==1 %--> Horizontal element to integrate on 
            
            G(i,j)=1/(4*pi)*((x_i(1)-x_j(1)+l_el)*log((x_i(2)-x_j(2))^2+(x_i(1)-x_j(1)+l_el)^2) - (x_i(1)-x_j(1)-l_el)*log((x_i(2)-x_j(2))^2+(x_i(1)-x_j(1)-l_el)^2) + 2*(x_i(2)-x_j(2))*atan((x_i(1)-x_j(1)+l_el)/(x_i(2)-x_j(2))) - 2*(x_i(2)-x_j(2))*atan((x_i(1)-x_j(1)-l_el)/(x_i(2)-x_j(2)))+4*l_el);            
            
        elseif X_el(j,4)==-1 %--> Vertical element to integrate on
            
            G(i,j)=1/(4*pi)*((x_i(2)-x_j(2)+l_el)*log((x_i(1)-x_j(1))^2+(x_i(2)-x_j(2)+l_el)^2) - (x_i(2)-x_j(2)-l_el)*log((x_i(1)-x_j(1))^2+(x_i(2)-x_j(2)-l_el)^2) + 2*(x_i(1)-x_j(1))*atan((x_i(2)-x_j(2)+l_el)/(x_i(1)-x_j(1))) - 2*(x_i(1)-x_j(1))*atan((x_i(2)-x_j(2)-l_el)/(x_i(1)-x_j(1)))+4*l_el);
        else
            disp('sdadaSDsdAAS$YY$EWEEHY$E$"YQYUQ"U$W"UW"UQ/')
        end
        
       

        
        
      
        
   end
end


G_inv=inv(G);



Sigma=G\phi_BC;

[Lf,Uf,PP] = lu(G);

Sigma_int=Lf\phi_BC;

Sigma_sol=Uf\Sigma_int

nanna=0;
nanna_sol=0;

for i=1:L
    
    TTf_sol=isnan(Sigma_sol(i));
    TTf=isnan(Sigma(i));
    if TTf_sol==1
        nanna_sol=nanna_sol+1;
        
    elseif TTf==1
        nanna=nanna+1
        
    end
end



save('G','G')

%% potential computation:
Phi_2d_BEM=zeros(M-1,N-1);



for i=3:N-2
    for j=3:M-2
        
        
        
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        
        if Phi0_2d(i,j)==0
        
        for l=1:1:L
            
            x_el=X_el(l,1:2);
            
            if X_el(l,4)==1
                
  Int_bound=1/(4*pi)*Sigma(l)*((x-x_el(1)+l_el)*log((y-x_el(2))^2+(x-x_el(1)+l_el)^2) - (x-x_el(1)-l_el)*log((y-x_el(2))^2+(x-x_el(1)-l_el)^2) + 2*(y-x_el(2))*atan((x-x_el(1)+l_el)/(y-x_el(2))) - 2*(y-x_el(2))*atan((x-x_el(1)-l_el)/(y-x_el(2)))+4*l_el);
            
            elseif X_el(l,4)==-1
  Int_bound=1/(4*pi)*Sigma(l)*((y-x_el(2)+l_el)*log((x-x_el(1))^2+(y-x_el(2)+l_el)^2) - (y-x_el(2)-l_el)*log((x-x_el(1))^2+(y-x_el(2)-l_el)^2) + 2*(x-x_el(1))*atan((y-x_el(2)+l_el)/(x-x_el(1))) - 2*(x-x_el(1))*atan((y-x_el(2)-l_el)/(x-x_el(1)))+4*l_el); 
            else
                disp('enculet')
                
            end
            
            
        Phi_2d_BEM(i,j)=Phi_2d_BEM(i,j)+Int_bound;
        
        
        end
        
        end
            
    
       
       
    end
end
toc
Phi_2d_BEM=Phi_2d_BEM./(max(max(Phi_2d_BEM)));
R1_boundedx_plot=R1_boundedx(1:N-1);
R1_boundedy_plot=R1_boundedx(1:M-1);


%BEM
figure()
 hold on
mesh(R1_boundedx_plot,R1_boundedy_plot,Phi_2d_BEM')
hold on
grid on
title('BEM potential')

figure()
 hold on
contour(R1_boundedx_plot,R1_boundedy_plot,Phi_2d_BEM',50)
hold on
grid on
title('BEM contour')

%Laplace
figure()
 hold on
mesh(R1_boundedx,R1_boundedy,Phi_2d')
hold on
grid on
title('Laplace potential')

figure()
 hold on
contour(R1_boundedx,R1_boundedy,Phi_2d',50)
grid on
hold on
title('Laplace contour')


%% ALP 
%parameters and global variables

setGlobalstep_grid(step_grid)
setGlobalPhi_2d(Phi_2d_BEM)


%step_time to reduce control action
step_time=5;
setGlobalstep_time(step_time)
flag_time=0;
setGlobalflag_time(flag_time )
time_before=0;
setGlobaltime_before(time_before)
time_last_fire=0;
setGlobaltime_last_fire(time_last_fire)

dt=1;
t_span=[0:dt:8*tau]';

Shape_Profile=0;

if Shape_Profile==0

alpha_v=norm(Rho_0_3d-X_obs_3d)/(1*t_end)
setGlobalalpha_v(alpha_v );
    
[Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done,I_obs]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_2d_step_time(dt,t_span,X_0_3d,X_obs_3d);
else
    
    F_sat=1*1e-3;
    setGlobalf_sat(F_sat)
    M_cubesat=5;
    setGlobalm_cubesat(M_cubesat)
    
    a_sat=F_sat/M_cubesat
    
[Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_2d_step_time_Shape_Pro(dt,t_span,X_0_3d);

end



figure()
contour(R1_boundedx, R1_boundedy,Phi_2d_contour')
hold on
plot3(Rho_an(:,1),Rho_an(:,2),Rho_an(:,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(end,1),Rho_an(end,2),Rho_an(end,3),'or',x_obs,y_obs,0,'*b')
title('collision Avoidance with discrete gradient & analytical trajectory propagation: BEM method')
grid on 

figure()
contour(R1_boundedx, R1_boundedy,Phi_2d_contour')
hold on
plot3(x_CW_FD(:,1),x_CW_FD(:,2),zebra_FD,'--b',x_CW_FD(1,1),x_CW_FD(1,2),zebra_FD,'og',x_CW_FD(end,1),x_CW_FD(end,2),zebra_FD,'or',x_obs,y_obs,0,'*b')
title('free drift')
grid on 

%results
% DV_tot=0;
%  for i=1:I_obs
% 
% DV_tot=norm(DV_v(i,:))+DV_tot;
% end

Goal_error=norm(Rho_an(end,:)'-X_obs_3d)
Total_Impulse=count_impulse_done

% DV_tot_obs=sum(norm(DV_v(1:I_obs,:)))
% TOF_obs=I_obs*dt

% Guidance Through laplace potential
setGlobalPhi_2d(Phi_2d)


Shape_Profile=0;

if Shape_Profile==0

alpha_v=norm(Rho_0_3d-X_obs_3d)/(1*t_end)
setGlobalalpha_v(alpha_v );
    
[Rho_an_Lap,Ni_an_Lap,DV_v,count_impulse_time,count_impulse_done,I_obs]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_2d_step_time(dt,t_span,X_0_3d,X_obs_3d);
else
    
    F_sat=1*1e-3;
    setGlobalf_sat(F_sat)
    M_cubesat=5;
    setGlobalm_cubesat(M_cubesat)
    
    a_sat=F_sat/M_cubesat
    
[Rho_an_Lap,Ni_an_Lap,DV_v,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_2d_step_time_Shape_Pro(dt,t_span,X_0_3d);

end

figure()
contour(R1_boundedx, R1_boundedy,Phi_2d_contour')
hold on
plot3(Rho_an_Lap(:,1),Rho_an_Lap(:,2),Rho_an_Lap(:,3),'--b',Rho_an_Lap(1,1),Rho_an_Lap(1,2),Rho_an_Lap(1,3),'og',Rho_an_Lap(end,1),Rho_an_Lap(end,2),Rho_an_Lap(end,3),'or',x_obs,y_obs,0,'*b')
title('collision Avoidance with discrete gradient & analytical trajectory propagation:Laplace Guidance')
grid on 

