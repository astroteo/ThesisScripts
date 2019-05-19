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
x_obs=50
y_obs=50
z_obs=0



lenght_ISS=20; %--> along x
width_ISS=-5; %--> along y
height_ISS=15;%--> along z


%[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,width_ISS,height_ISS);

%% Initial Conditions & Integration time 

% !!!! CRUCIAL TEST: by  givinig initial conditions such that x_obs_2d is
% reached no trajectory corrections are made. ====> OK !!! **
Rho_0_2d=[-20;-30;0];
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
     
lenght_ISS_CI=30; %--> along x !! reduced dymension for have faster computations !!
width_ISS_CI=10; %--> along y
height_ISS_CI=10;


L_box=100;

setGlobalL_box_x(L_box)
setGlobalL_box_y(L_box)
setGlobalL_box_z(L_box)

X_max=L_box;
Y_max=L_box;
Z_max=L_box;

step_grid=5;

R1_boundedx=-X_max:step_grid:X_max;
R1_boundedy=-Y_max:step_grid:Y_max;
R1_boundedz=-Z_max:step_grid:Z_max;

[~,N]=size(R1_boundedx);
[~,M]=size(R1_boundedy);
[~,P]=size(R1_boundedz);

if x_obs >=0
i_obs=floor(x_obs/step_grid)+floor((X_max/step_grid))+1
else
i_obs=floor((x_obs+floor(X_max/step_grid)*step_grid)/step_grid)+1
end

if y_obs >=0
j_obs=floor(y_obs/step_grid)+floor((Y_max/step_grid))+1
else
j_obs=floor((y_obs+floor(Y_max/step_grid)*step_grid)/step_grid)+1
end

if z_obs >=0
k_obs=floor(z_obs/step_grid)+floor((Z_max/step_grid))+1
else
k_obs=floor((z_obs+floor(Z_max/step_grid)*step_grid)/step_grid)+1
end



Phi0_3d=zeros(N,M,P);
v=0;
v_obj=0;

for i=1:N
    for j=1:M
        for k=1:P
        
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        z=R1_boundedz(k);
        
       if abs(z) < height_ISS_CI
            
            if y>=0
                
                  if abs(x)<=lenght_ISS_CI && abs(y) <= height_ISS_CI
                
                Phi0_3d(i,j,k)=1;
                
                v=v+1;
                Phi_vol(v,:)=[x,y,z];
                
                v_obj=v_obj+1;
                Phi_vol_obj(v_obj,:)=[x,y,z];
                
                  elseif abs(x) <=lenght_ISS_CI   && abs(x) >=2/3*lenght_ISS_CI && abs(y)<=3*height_ISS_CI
                     
                  Phi0_3d(i,j,k)=1;
                  
                v=v+1;
                Phi_vol(v,:)=[x,y,z];
                
                v_obj=v_obj+1;
                Phi_vol_obj(v_obj,:)=[x,y,z];
                      
                  end
                
            elseif y<0
                
                
                if abs(x)<=height_ISS_CI && abs(y) <= 3*height_ISS_CI
                
                Phi0_3d(i,j,k)=1;
                
                v=v+1;
                Phi_vol(v,:)=[x,y,z];
                
                v_obj=v_obj+1;
                Phi_vol_obj(v_obj,:)=[x,y,z];
                end
                
               
            end
            
            
          end 
            
            if i==N || i==1 || j==M || j==1 || k==1 || k==P
           
           Phi0_3d(i,j,k)=1;
           
           v=v+1;
           Phi_vol(v,:)=[x,y,z];
           
           
            end
            
            
            if i==i_obs && j==j_obs && k==k_obs
                
            Phi0_3d(i,j,k)=-1;% only to visualize the goal point to be set to zero in the Laplace solution
            
            v=v+1;
            Phi_vol(v,:)=[x,y,z];
            
           
            
            end
            
            
        
            
       end      
     end
end
        
        
        

figure()
hold on
plot3 (Phi_vol(:,1),Phi_vol(:,2),Phi_vol(:,3),'*b')
title('volume points')
axis equal
xlabel('x')
ylabel('y')
grid on  
                                        
figure()
hold on
plot3 (Phi_vol_obj(:,1),Phi_vol_obj(:,2),Phi_vol_obj(:,3),'*r')
hold on
plot3 (x_obs,y_obs,z_obs,'*k')
title('object volume points')
axis equal
xlabel('x')
ylabel('y')
grid on  

%% Countour (Bounduary-generation)
    
Phi_3d_contour=zeros(N,M,P);
n_el_goal=1


s=0;
s_obj=0;
s_goal=0;

BOX_check=0;

for i=1:N
    for j=1:M
      for k=1:P  
       
        
        
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        z=R1_boundedz(k);
        
        
        
if abs(z)<height_ISS_CI 
       
     if y>=0
                
                  if abs(x)==lenght_ISS_CI  && abs(y) <= height_ISS_CI 
                
                      Phi_3d_contour(i,j,k)=1;
                      
                      s=s+1;
                      Phi_surf(s,:)=[x,y,z];
                      
                      s_obj=s_obj+1;
                      Phi_surf_obj(s_obj,:)=[x,y,z];
                
                  elseif abs(x)<=lenght_ISS_CI && abs(x) <=2/3*(lenght_ISS_CI) &&  abs(y) == height_ISS_CI
                      
                      Phi_3d_contour(i,j,k)=1;
                      
                      s=s+1;
                      Phi_surf(s,:)=[x,y,z];
                      
                      s_obj=s_obj+1;
                      Phi_surf_obj(s_obj,:)=[x,y,z];
                
                  elseif abs(x) ==lenght_ISS_CI && abs(y)<=3*(height_ISS_CI) && abs(y) >= (height_ISS_CI) || abs(x) ==2/3*(lenght_ISS_CI)  && abs(y)<=3*height_ISS_CI && abs(y) >=(height_ISS_CI)
                     
                   Phi_3d_contour(i,j,k)=1;
                  
                      s=s+1;
                      Phi_surf(s,:)=[x,y,z];
                      
                      s_obj=s_obj+1;
                      Phi_surf_obj(s_obj,:)=[x,y,z];
                  
                  elseif abs(x) <=(lenght_ISS_CI)  &&  abs(x) >=2/3*(lenght_ISS_CI)  && abs(y)==3*(height_ISS_CI)
                     
                    Phi_3d_contour(i,j,k)=1;
                  
                      s=s+1;
                      Phi_surf(s,:)=[x,y,z];
                      
                      s_obj=s_obj+1;
                      Phi_surf_obj(s_obj,:)=[x,y,z];
                  
                  
                  elseif abs(x)==lenght_ISS_CI  && abs(y) == 0
                
                     Phi_3d_contour(i,j,k)=1;
                
                      s=s+1;
                      Phi_surf(s,:)=[x,y,z];
                      
                      s_obj=s_obj+1;
                      Phi_surf_obj(s_obj,:)=[x,y,z];
                  end
                
                if abs(x) <=(lenght_ISS_CI) && abs(x)>=(height_ISS_CI) && y ==0
                     
                 Phi_3d_contour(i,j,k)=1;
                  
                      s=s+1;
                      Phi_surf(s,:)=[x,y,z];
                      
                      s_obj=s_obj+1;
                      Phi_surf_obj(s_obj,:)=[x,y,z];
                  
                end
                  
                      
             
                
     elseif y<0
                
                
                if abs(x)==(height_ISS_CI) &&  abs(y) <= 3*(height_ISS_CI) 
                
                 Phi_3d_contour(i,j,k)=1;
                
                       s=s+1;
                      Phi_surf(s,:)=[x,y,z];
                      
                      s_obj=s_obj+1;
                      Phi_surf_obj(s_obj,:)=[x,y,z];
                end
                
                
                
                
                if abs(x)<(height_ISS_CI)  &&  abs(y) == 3*(height_ISS_CI)
                
                Phi_3d_contour(i,j,k)=1;
                
                 s=s+1;
                 Phi_surf(s,:)=[x,y,z];
                 
                 s_obj=s_obj+1;
                 Phi_surf_obj(s_obj,:)=[x,y,z];
                end
                
               
    
     end
     
elseif abs(z)==height_ISS_CI
    
    if y>=0
                
                  if abs(x)  <=2/3*lenght_ISS_CI   &&  abs(y) <= height_ISS_CI
                
                Phi_3d_contour(i,j,k)=1;
                
                s_obj=s_obj+1;
                Phi_surf_obj(s_obj,:)=[x,y,z];
                
                  elseif abs(x) <=lenght_ISS_CI   && abs(x) >= 2/3*lenght_ISS_CI && abs(y)<=3*height_ISS_CI 
                     
                Phi_3d_contour(i,j,k)=1;
                  
                 s_obj=s_obj+1;
                 Phi_surf_obj(s_obj,:)=[x,y,z];
                      
                  end
                
      elseif y<0
                
                
                if abs(x)<=height_ISS_CI && abs(y) <= 3*height_ISS_CI
                
                Phi_3d_contour(i,j,k)=1;
                
                s_obj=s_obj+1;
                 Phi_surf_obj(s_obj,:)=[x,y,z];
                end
                
               
     end
         
end

    %B0X
    
    
                    if i==N-1 ||i==2|| j==M-1 ||j==2 || k==2 || k==P-1
           
                     Phi_3d_contour(i,j,k)=1;
                     
                     s=s+1;
                    Phi_surf(s,:)=[x,y,z];
                   
                    end
                    
                    
                    if i==N-2 ||i==2|| j==M-2 ||j==2 
                        
                        BOX_check=BOX_check+1;
                        Phi_BOX_check(BOX_check,:)=[x,y,z];
                    end
                    
                    
                    
     %GOAL
                   
                    if ((i==i_obs+n_el_goal || i==i_obs-n_el_goal) && (j>=j_obs-n_el_goal   && j<=j_obs+n_el_goal  && k>=k_obs-n_el_goal   && k<=k_obs+n_el_goal)) 
               
                  Phi_3d_contour(i,j,k)=-1;
                           
                  s_goal=s_goal+1;
                  Phi_surf_goal(s,:)=[x,y,z];
                  
                  
                    elseif ((j==j_obs+n_el_goal || j==j_obs-n_el_goal) && (i>=i_obs-n_el_goal   && i<=i_obs+n_el_goal  && k>=k_obs-n_el_goal   && k<=k_obs+n_el_goal)) 
                        
                        
                  Phi_3d_contour(i,j,k)=-1;
                           
                  s_goal=s_goal+1;
                  Phi_surf_goal(s,:)=[x,y,z];
                  
                    elseif ( (k==k_obs+n_el_goal || k==k_obs-n_el_goal) && (j>=j_obs-n_el_goal   && j<=j_obs+n_el_goal  && i>=i_obs-n_el_goal   && i<=i_obs+n_el_goal))
                  
                  Phi_3d_contour(i,j,k)=-1;
                           
                  s_goal=s_goal+1;
                  Phi_surf_goal(s,:)=[x,y,z];
                  
                    end
                    
                  
                    

    

                    
      end
   end
end
    
figure()
 hold on
plot3 (Phi_surf(:,1),Phi_surf(:,2),Phi_surf(:,3),'*b')
title('contour with goal')
xlabel('x')
ylabel('y')
grid on  
axis equal

figure()
 hold on
plot3 (Phi_surf_obj(:,1),Phi_surf_obj(:,2),Phi_surf_obj(:,3),'*b')
 hold on
plot3 (Phi_surf_goal(:,1),Phi_surf_goal(:,2),Phi_surf_goal(:,3),'*b')
hold on
plot3 (x_obs,y_obs,z_obs,'*k')
title('object contour with goal ')
xlabel('x')
ylabel('y')
grid on  
axis equal

figure()
 hold on
plot3 (Phi_BOX_check(:,1),Phi_BOX_check(:,2),Phi_BOX_check(:,3),'*b')
hold on
plot3 (x_obs,y_obs,z_obs,'*k')
title('box check ')
xlabel('x')
ylabel('y')
grid on  
axis equal

%% element generation from contour
I=0;
J=0;
K=0;

L=0;

I_goal=0;
J_goal=0;
K_goal=0;

L_goal=0; % X_el=[x_el,y_el,z_el,u,v,w,1/0 obstacle/goal] 

L_obj=0;


for i=3:1:N-3
    for j=3:1:M-3
        for k=3:1:P-3 %!! remember to move/double the box for {i,j,k}=2:{i,j,k}=[N-1,M-1;P-1] !!
        
        %  exclude vertices: IF Phi_3d_contour(i,j,k)==1 && (Phi_3d_contour(i+1,j,k)~=0 && Phi_3d_contour(i-1,j,k)~=0) || (Phi_3d_contour(i,j+1,k)~=0 && Phi_3d_contour(i,j-1,k)~=0 )|| (Phi_3d_contour(i,j,k+1)~=0 && Phi_3d_contour(i,j,k-1)~=0)
        %  exclude edges:    IF Phi_3d_contour(i,j,k)==1 && (Phi_3d_contour(i+1,j+1,k)~=0 && Phi_3d_contour(i-1,j-1,k)~=0) || (Phi_3d_contour(i+1,j,k+1)~=0 && Phi_3d_contour(i-1,j,k-1)~=0 )|| (Phi_3d_contour(i,j+1,k+1)~=0 && Phi_3d_contour(i,j-1,k-1)~=0)
        
        %           ALT:     IF Phi_3d_contour(i,j,k)==-1 &&( (Phi_3d_contour(i+1,j+1,k)~=0 && Phi_3d_contour(i-1,j-1,k)~=0 && Phi_3d_contour(i+1,j-1,k)~=0 && Phi_3d_contour(i-1,j+1,k)~=0) || (Phi_3d_contour(i+1,j,k+1)~=0 && Phi_3d_contour(i-1,j,k-1)~=0 && Phi_3d_contour(i+1,j,k-1)~=0 && Phi_3d_contour(i-1,j,k+1)~=0 )|| (Phi_3d_contour(i,j+1,k+1)~=0 && Phi_3d_contour(i,j-1,k-1)~=0 && Phi_3d_contour(i,j-1,k+1)~=0 && Phi_3d_contour(i,j+1,k-1)~=0) ) 
        %   SEEMS 2 WORK:    IF Phi_3d_contour(i,j,k)==1 &&( (Phi_3d_contour(i+1,j+1,k)~=0 && Phi_3d_contour(i-1,j-1,k)~=0 && Phi_3d_contour(i+1,j-1,k)~=0 && Phi_3d_contour(i-1,j+1,k)~=0) || (Phi_3d_contour(i+1,j,k+1)~=0 && Phi_3d_contour(i-1,j,k-1)~=0 && Phi_3d_contour(i+1,j,k-1)~=0 && Phi_3d_contour(i-1,j,k+1)~=0 )|| (Phi_3d_contour(i,j+1,k+1)~=0 && Phi_3d_contour(i,j-1,k-1)~=0 && Phi_3d_contour(i,j-1,k+1)~=0 && Phi_3d_contour(i,j+1,k-1)~=0) )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
                                                                     
        if Phi_3d_contour(i,j,k)==1 &&( (Phi_3d_contour(i+1,j+1,k)~=0 && Phi_3d_contour(i-1,j-1,k)~=0 && Phi_3d_contour(i+1,j-1,k)~=0 && Phi_3d_contour(i-1,j+1,k)~=0) || (Phi_3d_contour(i+1,j,k+1)~=0 && Phi_3d_contour(i-1,j,k-1)~=0 && Phi_3d_contour(i+1,j,k-1)~=0 && Phi_3d_contour(i-1,j,k+1)~=0 )|| (Phi_3d_contour(i,j+1,k+1)~=0 && Phi_3d_contour(i,j-1,k-1)~=0 && Phi_3d_contour(i,j-1,k+1)~=0 && Phi_3d_contour(i,j+1,k-1)~=0) )
                        
                     L=L+1;
        
                        X_el(L,1:3)=[R1_boundedx(i),R1_boundedy(j),R1_boundedz(k)];
                        X_el(L,7)=1;
                        
                     L_obj=L_obj+1;
                        
                        X_el_obj(L_obj,1:3)=[R1_boundedx(i),R1_boundedy(j),R1_boundedz(k)];
                        
                        if     (Phi_3d_contour(i+1,j,k) ==1 ||  Phi_3d_contour (i-1,j,k)==1) && ( Phi_3d_contour(i,j+1,k)==1 || Phi_3d_contour(i,j-1,k)==1)
                            
                            if Phi0_3d(i,j,k+1)==0 && Phi0_3d(i,j,k-1)==1
                                X_el(L,4:6)=[0,0,1];
                            else
                                X_el(L,4:6)=[0,0,-1];
                            end
                            
                        elseif (Phi_3d_contour(i+1,j,k) ==1 ||  Phi_3d_contour (i-1,j,k)==1) && ( Phi_3d_contour(i,j,k+1)==1 || Phi_3d_contour(i,j,k-1)==1)
                            
                            if Phi0_3d(i,j+1,k)==0 && Phi0_3d(i,j-1,k)==1
                                X_el(L,4:6)=[0,1,0];
                            else
                                X_el(L,4:6)=[0,-1,0];
                            end
                            
                        elseif (Phi_3d_contour(i,j+1,k) ==1 ||  Phi_3d_contour (i,j-1,k)==1) && ( Phi_3d_contour(i,j,k+1)==1 || Phi_3d_contour(i,j,k-1)==1)
                            
                             if Phi0_3d(i+1,j,k)==0 && Phi0_3d(i-1,j,k)==1
                                X_el(L,4:6)=[1,0,0];
                            else
                                X_el(L,4:6)=[-1,0,0];
                            end
                        
                        end
                    
                    
                    
                        
        elseif Phi_3d_contour(i,j,k)==-1 &&( (Phi_3d_contour(i+1,j+1,k)~=0 && Phi_3d_contour(i-1,j-1,k)~=0 && Phi_3d_contour(i+1,j-1,k)~=0 && Phi_3d_contour(i-1,j+1,k)~=0) || (Phi_3d_contour(i+1,j,k+1)~=0 && Phi_3d_contour(i-1,j,k-1)~=0 && Phi_3d_contour(i+1,j,k-1)~=0 && Phi_3d_contour(i-1,j,k+1)~=0 )|| (Phi_3d_contour(i,j+1,k+1)~=0 && Phi_3d_contour(i,j-1,k-1)~=0 && Phi_3d_contour(i,j-1,k+1)~=0 && Phi_3d_contour(i,j+1,k-1)~=0) )             
            
                    L=L+1;
        
                        X_el(L,1:3)=[R1_boundedx(i),R1_boundedy(j),R1_boundedz(k)];
                        X_el(L,7)=0;
                        
                        if     (Phi_3d_contour(i+1,j,k) ==-1 ||  Phi_3d_contour (i-1,j,k)==-1) && ( Phi_3d_contour(i,j+1,k)==-1 || Phi_3d_contour(i,j-1,k)==-1)
                            
                            if dot((X_el(L,1:3)-X_obs_3d'),[0,0,1]) > dot((X_el(L,1:3)-X_obs_3d'),[0,0,-1])
                                X_el(L,4:6)=[0,0,1];
                            else
                                X_el(L,4:6)=[0,0,-1];
                            end
                            
                        elseif (Phi_3d_contour(i+1,j,k) ==-1 ||  Phi_3d_contour (i-1,j,k)==-1) && ( Phi_3d_contour(i,j,k+1)==-1 || Phi_3d_contour(i,j,k-1)==-1)
                            
                           if dot((X_el(L,1:3)-X_obs_3d'),[0,1,0]) > dot((X_el(L,1:3)-X_obs_3d'),[0,-1,0])
                                X_el(L,4:6)=[0,1,0];
                            else
                                X_el(L,4:6)=[0,-1,0];
                            end
                            
                        elseif (Phi_3d_contour(i,j+1,k) ==-1 ||  Phi_3d_contour (i,j-1,k)==-1) && ( Phi_3d_contour(i,j,k+1)==-1 || Phi_3d_contour(i,j,k-1)==-1)
                            
                            if dot((X_el(L,1:3)-X_obs_3d'),[1,0,0]) > dot((X_el(L,1:3)-X_obs_3d'),[-1,0,0])
                                X_el(L,4:6)=[1,0,0];
                            else
                                X_el(L,4:6)=[-1,0,0];
                            end
                        
                        end
                        
                        
                        
                    L_goal=L_goal+1;
                    
                    X_el_goal(L_goal,1:3)=[R1_boundedx(i),R1_boundedy(j),R1_boundedz(k)];
                    
                    
                   
        
        end
        
        
         
            
            
       
        
        end    
     end
end


figure()
hold on
plot3 (Phi_surf_obj(:,1),Phi_surf_obj(:,2),Phi_surf_obj(:,3),'*r')
hold on
plot3(X_el(:,1),X_el(:,2),X_el(:,3),'*k')
hold on
plot3(X_obs_3d(1),X_obs_3d(2),X_obs_3d(3),'*g')
grid on
axis equal
title('All elements')


figure()
hold on
quiver3 (X_el(:,1),X_el(:,2),X_el(:,3),X_el(:,4),X_el(:,5),X_el(:,6))
hold on
plot3(X_el(:,1),X_el(:,2),X_el(:,3),'*k')
hold on
plot3(X_obs_3d(1),X_obs_3d(2),X_obs_3d(3),'*g')
grid on
axis equal
title('All elements')


figure()
 hold on
 plot3(X_el_obj(:,1),X_el_obj(:,2),X_el_obj(:,3),'*r')
hold on
plot3(X_el_goal(:,1),X_el_goal(:,2),X_el_goal(:,3),'*b')
grid on
axis equal
title('goal and object elements')


figure()
 hold on
 plot(X_el_obj(:,1),X_el_obj(:,2),'*r')
hold on
plot(X_el_goal(:,1),X_el_goal(:,2),'*b')
hold on
grid on
axis equal
title('goal and object elements: top view')

%% 
Solve_Laplace=0;

%% Laplace equation discrete solution
if Solve_Laplace==1

Phi_3d=ones(N,M,P);
Phi_3d_new=ones(N,M,P);

 epsilon=1;
 count=0;
 
 displaY_count=0;
 
while epsilon >= 0.1 && count <=500

 for i=2:N-1
     for j=2:M-1
         for k=2:P-1;
             
        if Phi0_3d(i,j,k)==1
            
            Phi_3d(i,j,k)=1;
            
        elseif Phi0_3d(i,j,k)==-1
            
            Phi_3d(i,j,k)=0;
            
        end
           
           
           
           
  Phi_3d_new(i,j,k)=1/6*(Phi_3d(i-1,j,k)+Phi_3d(i+1,j,k)+Phi_3d(i,j-1,k)+Phi_3d(i,j+1,k)+Phi_3d(i,j,k-1)+Phi_3d(i,j,k+1));
  
  
     
         end    
     end
 end
 
 
count=count+1;

if count-displaY_count >= 10
    
    displaY_count=displaY_count+10
    epsilon
end

epsilon=max(max(max(abs(Phi_3d_new-Phi_3d))));


Phi_3d=Phi_3d_new;

end

Phi_2d_xy=Phi_3d(:,:,k_obs);
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
contour(R1_boundedx,R1_boundedz,Phi_2d_xz')
colormap summer
hold on
plot(y_obs,z_obs,'*b')
title('Laplace potential contour x-z')
grid on 


%% ALP method-->laplace guidance

%parameters and global variables

setGlobalstep_grid(step_grid)
setGlobalPhi_3d(Phi_3d)

%step_time to reduce control action
step_time=1;
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

Shape_Velocity=0

if Shape_Velocity==0

   alpha_v=norm(Rho_0_3d-X_obs_3d)/(1*t_end);
   setGlobalalpha_v(alpha_v );
   
  [Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_3d_step_time(dt,t_span,X_0_3d);

else
    
    F_sat=1e-3;
    setGlobalf_sat(F_sat)
    M_cubesat=5;
    setGlobalm_cubesat(M_cubesat)

  [Rho_an,Ni_an,DV_v,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_3d_step_time_Shape_Pro(dt,t_span,X_0_3d);

    
end


figure()
 hold on
plot3 (Phi_vol_obj(:,1),Phi_vol_obj(:,2),Phi_vol_obj(:,3),'*r')
hold on
plot3 (x_obs,y_obs,z_obs,'*k')
hold on
plot3(Rho_an(:,1),Rho_an(:,2),Rho_an(:,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(end,1),Rho_an(end,2),Rho_an(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('collision Avoidance with discrete gradient & analytical trajectory propagation')
xlabel('x')
ylabel('y')
grid on  
axis equal
end