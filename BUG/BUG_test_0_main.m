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
y_obs=-60+3
z_obs=0

X_goal=[x_obs;y_obs;z_obs];

lenght_ISS=50; %--> along x
width_ISS=25; %--> along y
height_ISS=15;%--> along z


%[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,width_ISS,height_ISS);

%% Initial Conditions & Integration time 

Rho_0_2d=[-40;40];
Rho_0_3d=[Rho_0_2d(1);Rho_0_2d(2);0];

tau=1/4*T_orb;


X_obs_3d=[x_obs;y_obs;z_obs];

Ct=cos(w*tau);
St=sin(w*tau);

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

Rho_0=Rho_0_3d;
Ni_0=zeros(3,1);

%% Free Drift !! UNforced Problem!!

t_end=tau;

options = odeset('RelTol',1e-13,'AbsTol', 1e-12);
[t_CW_FD,x_CW_FD]= ode113('Chol_Wilt_Hill',t_end,X_0_3d,options);%(plot at the end of the script)

[CW_FD,~]=size(t_CW_FD);
zebra_FD=zeros(CW_FD,1);

%% 2d Discrete Potential C+I shaped object 
     % Comparison with ALP method
     

L_max=norm(Rho_0_3d);  

X_max=4*abs(Rho_0_3d(1));
Y_max=4*abs(Rho_0_3d(2));

step_grid=0.25;

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
Phi_2d_contour=zeros(N,M);

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
          
            end
            
            
           %contour  
           if y>=0
                
                  if abs(x)==lenght_ISS   && abs(y) <= height_ISS  
                
                      Phi_2d_contour(i,j)=1;
                
                  elseif abs(x)<=lenght_ISS && abs(x) <=3/2*width_ISS &&  abs(y) == height_ISS 
                      
                      Phi_2d_contour(i,j)=1;
                
                  elseif abs(x) ==lenght_ISS && abs(y)<=3/2*width_ISS && abs(y) >= height_ISS || abs(x) ==3/2*width_ISS  && abs(y)<=3/2*width_ISS && abs(y) >= height_ISS
                     
                  Phi_2d_contour(i,j)=1;
                  
                  elseif abs(x) <=lenght_ISS  &&  abs(x) >=3/2*width_ISS  && abs(y)==3/2*width_ISS
                     
                  Phi_2d_contour(i,j)=1;
                  
                  
                  elseif abs(x) <=lenght_ISS && abs(x)>=height_ISS/2  && y==0 
                     
                  Phi_2d_contour(i,j)=1;
                  
                  
                      
                  end
                
            elseif y<0
                
                
                if abs(x)==height_ISS/2 &&  abs(y) <= 2*width_ISS
                
                Phi_2d_contour(i,j)=1;
                end
                
                if abs(x)<=height_ISS/2  &&  abs(y) == 2*width_ISS
                
                Phi_2d_contour(i,j)=1;
                end
                
               
            end
            
            if i==N || i==1 || j==M || j==1
           
           Phi_2d_contour(i,j)=1;
            
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

figure()
hold on
mesh(R1_boundedx,R1_boundedy,Phi_2d_contour')
title('bounduary conditions contour')
xlabel('x')
ylabel('y')
grid on  

[grad_Phi0_2d_x,grad_Phi0_2d_y]=gradient(Phi0_2d);
figure()
quiver(R1_boundedx,R1_boundedy,-grad_Phi0_2d_y',-grad_Phi0_2d_x')
hold on
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
title(' Quiver on bounduary conditions')
xlabel('x')
ylabel('y')
grid on  


% Phi_2d=1*ones(M,N)+Phi_BOX;
% Phi_2d_new=1*ones(M,N)+Phi_BOX;


Phi_2d=1*ones(M,N);
Phi_2d_new=1*ones(M,N);
epsilon=1;
count=0;
%% SELCET ALGORITHM

select_Alg=1
                  
                  %0 --> Bug-0 algorithm with with fixed tau_m 
                  %1 --> real time Bug-0 algorithm ~ unknown obstacle
                  %2 --> Bug-2 filter for all possible tau junction on the original m-line
                  %3 --> [filter Bug-2] before Bug-0
                 




%% ~BUG-0 algorith 1-single m-line at once && @fixed tau

if select_Alg==0
    
tau_bug_0=1*tau

tau_bug_trasl= tau_bug_0/32;

flag_bug_0=1; count_bug_0=0;

%initial m-line
dt_m_bug_0=10;
dt_bug_0=1




[Ni_0_plus,~]=Compute_ci(tau_bug_0,Rho_0,Ni_0,X_goal);
[Rho_m,Ni_m]=Propagate_ci(dt_m_bug_0,tau_bug_0,Rho_0,Ni_0_plus);

figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_m(:,1),Rho_m(:,2),Rho_m(:,3),'ok',Rho_m(1,1),Rho_m(1,2),Rho_m(1,3),'og',Rho_m(end,1),Rho_m(end,2),Rho_m(end,3),'or')
title('Initial (0) m-line ')
grid on 


flag_coll_mm=1;

R_gnc=10;

[Rho_th,~]=Propagate_ci(dt_bug_0,dt_m_bug_0,Rho_0,Ni_0_plus);
R_gnc_safety=norm(Rho_th(end,:)-Rho_0')

Rho_TRAJ=Rho_0';
Ni_TRAJ=Ni_0';

l=0;
TOF=0;

while flag_bug_0==1 && count_bug_0 <=150
    
    
 [Ni_0_plus,~]=Compute_ci(tau_bug_0,Rho_0,Ni_0,X_goal);
 [Rho_m,Ni_m]=Propagate_ci(dt_m_bug_0,tau_bug_0,Rho_0,Ni_0_plus);


[flag_bug_0]=Evaluate_collision(Rho_m,Phi0_2d,step_grid);
    
    if flag_bug_0==1
         %[Rho_m_star,O_k,K]=Tangent_bug(Rho_m,Phi_2d_contour,step_grid,R_gnc);
         %[Rho_m_star,O_k,K]=Tangent_bug_New(Rho_m,Rho_0,Ni_0,tau_bug_0,Phi_2d_contour,step_grid,R_gnc);
         %[Rho_m_star,O_k,K]=Tangent_bug_2Pt(Rho_m,Phi_2d_contour,step_grid,R_gnc);
         [Rho_m_star,Ni_m_star,O_k,K]=Tangent_bug_2Pt_velocity(Rho_m,Ni_m,Phi_2d_contour,step_grid,R_gnc);
          

          %[Rho_0_new]=Heuristic(Rho_m_star,K,O_k,R_gnc);
          %[Rho_0_new]=Heuristic_max(Rho_m_star,K,O_k,R_gnc);
          %[Rho_0_new]=Real_Heuristic(Rho_m_star,K,O_k,X_goal,R_gnc);
          %[Rho_0_new]=Heuristic_tangent(Rho_m_star,K,O_k,R_gnc);
          %[Rho_0_new]=Heuristic_tangent_smart(Rho_m_star,K,O_k,R_gnc,X_goal);
          %[Rho_0_new]=Heuristic_perpendicular(Rho_m_star,K,O_k,R_gnc,X_goal);
          %[Rho_0_new]=Heuristic_tangent_smart_New(Rho_m_star,K,O_k,R_gnc,X_goal,Phi0_2d,step_grid);
          %[Rho_0_new,n_t]=Heuristic_tangent_smart_ortodox(Rho_m_star,Rho_0,K,O_k,R_gnc,X_goal,Phi0_2d,step_grid);
          %[Rho_0_new,n_t]=Heuristic_tangent_smart_ortodox_2Pt(Rho_m_star,Rho_0,K,O_k,R_gnc,X_goal,Phi0_2d,step_grid);
           [Rho_0_new,n_t]=Heuristic_tangent_smart_ortodox_2Pt_velocity(Rho_m_star,Rho_0,Ni_0,K,O_k,R_gnc,X_goal,Phi0_2d,step_grid);
          
      



         
    else
      
        Rho_0_new=X_goal;
        
    end
    
    
  if flag_bug_0==1
      
[Ni_0_plus,~]=Compute_ci(tau_bug_trasl,Rho_0,Ni_0,Rho_0_new);
[Rho_traj,Ni_traj]=Propagate_ci(dt_bug_0,tau_bug_trasl,Rho_0,Ni_0_plus);

TOF=TOF+tau_bug_trasl;
  else
      
[Ni_0_plus,~]=Compute_ci(tau_bug_0,Rho_0,Ni_0,Rho_0_new);
[Rho_traj,Ni_traj]=Propagate_ci(dt_bug_0,tau_bug_0,Rho_0,Ni_0_plus);

TOF=TOF+tau_bug_0;
      
  end
  l=l+1;
  DV_v_BUG_0(l,:)=Ni_traj(1,:)-Ni_TRAJ(end,:);
  
  

Rho_TRAJ=[Rho_TRAJ;Rho_traj];
Ni_TRAJ=[Ni_TRAJ;Ni_traj];


% Rho_0_new=Rho_0;% !! FAKE !! only to test.
%  Ni_0_new= Ni_0 ;

Rho_0=Rho_TRAJ(end,:)';
Ni_0=Ni_TRAJ(end,:)';




count_bug_0=count_bug_0+1; % to be eliminated when speed will be DECENT
end

DV_tot_BUG_0=0;
for i=1:l
    
DV_tot_BUG_0=DV_tot_BUG_0+norm(DV_v_BUG_0(i,:));

end

if K>3
figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_m(:,1),Rho_m(:,2),Rho_m(:,3),'ok',Rho_m(1,1),Rho_m(1,2),Rho_m(1,3),'og',Rho_m(end,1),Rho_m(end,2),Rho_m(end,3),'or')
plot3(O_k(1:K,1),O_k(1:K,2),O_k(1:K,3),'--r',O_k(1,1),O_k(1,2),O_k(1,3),'ok',O_k(2,1),O_k(2,2),O_k(2,3),'ob',O_k(K-1,1),O_k(K-1,2),O_k(K-1,3),'ok',O_k(K,1),O_k(K,2),O_k(K,3),'ob',Rho_m_star(1),Rho_m_star(2),Rho_m_star(3),'*r',Rho_0_new(1),Rho_0_new(2),Rho_0_new(3),'+k')% !! Rho_0_new_Heu to be replaced with Rho_0_new !!
title('Algorithm meaning')
grid on 
else
figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_m(:,1),Rho_m(:,2),Rho_m(:,3),'ok',Rho_m(1,1),Rho_m(1,2),Rho_m(1,3),'og',Rho_m(end,1),Rho_m(end,2),Rho_m(end,3),'or')
title('Algorithm meaning')
grid on   
end

figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_TRAJ(:,1),Rho_TRAJ(:,2),Rho_TRAJ(:,3),'--b',Rho_TRAJ(1,1),Rho_TRAJ(1,2),Rho_TRAJ(1,3),'og',Rho_TRAJ(end,1),Rho_TRAJ(end,2),Rho_TRAJ(end,3),'or',x_obs,y_obs,0,'*b')
title('trajectory ')
grid on 

%% ~BUG-0 Pseudo- real time algorithm 
elseif select_Alg==1

tau_bug_rt=1*tau;

tau_bug_trasl= tau_bug_rt/4;

R_near=2; 

Delta_R= norm(X_goal-Rho_0);

if Delta_R <= R_near
    
flag_bug_rt=0; 

else
    
flag_bug_rt=1; 
end
    

count_bug_rt=0;

%initial m-line
dt_m_bug_rt=10;
dt_bug_rt=2

R_gnc=15;

R_gnc_sens=15;

Rho_TRAJ=Rho_0';
Ni_TRAJ=Ni_0';

flag_control=0;

l=0;

TOF_rt=0;

count_new_goal=0;

while flag_bug_rt==1 && count_bug_rt <=200
    
    
   % initial conditions selection
   if flag_control==1
       
 [Ni_0_plus,~]=Compute_ci(tau_bug_rt,Rho_0,Ni_0,X_goal);
 
 l=l+1;
 DV_v_rt(l,:)=Ni_0_plus-Ni_0;
 
 
   else
 
 Ni_0_plus=Ni_0;
       
   end
   
   
 [Rho_traj,Ni_traj]=Propagate_ci(dt_bug_rt,dt_m_bug_rt,Rho_0,Ni_0_plus);

 TOF_rt=TOF_rt+dt_m_bug_rt;


 

       [~,O_k,K]=Tangent_bug_2Pt(Rho_TRAJ(end,:),Phi_2d_contour,step_grid,R_gnc_sens);
       [Rho_0_new,n_t,flag_control]=Heuristic_tangent_smart_ortodox_2Pt_rt(Rho_TRAJ(end,:)',Rho_TRAJ(end,:)',Ni_0,K,O_k,R_gnc,X_goal,Phi0_2d,step_grid);
       
       

    if flag_control==1
        
   
 
    [Ni_0_plus,~]=Compute_ci(tau_bug_trasl,Rho_0,Ni_0,Rho_0_new);
    [Rho_traj,Ni_traj]=Propagate_ci(dt_bug_rt,tau_bug_trasl,Rho_0,Ni_0_plus);
    
    disp('new goal')
    
    count_new_goal=count_new_goal+1;
    
     l=l+1;
     DV_v_rt(l,:)=Ni_0_plus-Ni_0;
    
     TOF_rt=TOF_rt+tau_bug_trasl;
 
   end
 
  

Rho_TRAJ=[Rho_TRAJ;Rho_traj];
Ni_TRAJ=[Ni_TRAJ;Ni_traj];


[N_TRAJ,~]= size(Rho_TRAJ);

for i_TRAJ=1:N_TRAJ
    
Delta_R= norm(X_goal-Rho_TRAJ (i_TRAJ,:)');

if Delta_R <= R_near
    
flag_bug_rt=0; 

i_TRAJ_goal=i_TRAJ;

end

end

Rho_0=Rho_TRAJ(end,:)';
Ni_0=Ni_TRAJ(end,:)';

    
count_bug_rt=count_bug_rt+1; % to be eliminated when speed will be DECENT
end

% figure()
% contour(R1_boundedx, R1_boundedy,Phi0_2d')
% hold on
% plot3(Rho_m(:,1),Rho_m(:,2),Rho_m(:,3),'ok',Rho_m(1,1),Rho_m(1,2),Rho_m(1,3),'og',Rho_m(end,1),Rho_m(end,2),Rho_m(end,3),'or')
% plot3(O_k(1:K,1),O_k(1:K,2),O_k(1:K,3),'--r',Rho_m_star(1),Rho_m_star(2),Rho_m_star(3),'*r',Rho_0_new(1),Rho_0_new(2),Rho_0_new(3),'+k')% !! Rho_0_new_Heu to be replaced with Rho_0_new !!
% title('Algorithm meaning')
% grid on 

DV_tot_rt=0;
for i=1:l
    
DV_tot_rt=DV_tot_rt+norm(DV_v_rt(i,:));

end
if flag_bug_rt==0
figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_TRAJ(1:i_TRAJ_goal,1),Rho_TRAJ(1:i_TRAJ_goal,2),Rho_TRAJ(1:i_TRAJ_goal,3),'--b',Rho_TRAJ(1,1),Rho_TRAJ(1,2),Rho_TRAJ(1,3),'og',Rho_TRAJ(i_TRAJ_goal,1),Rho_TRAJ(i_TRAJ_goal,2),Rho_TRAJ(i_TRAJ_goal,3),'or',x_obs,y_obs,0,'*b')
title('trajectory ')
grid on 
else
figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_TRAJ(:,1),Rho_TRAJ(:,2),Rho_TRAJ(:,3),'--b',Rho_TRAJ(1,1),Rho_TRAJ(1,2),Rho_TRAJ(1,3),'og',Rho_TRAJ(end,1),Rho_TRAJ(end,2),Rho_TRAJ(end,3),'or',x_obs,y_obs,0,'*b')
title('trajectory ')
grid on 
    
end
%% ~BUG-2 Algorithm @ variable tau

elseif select_Alg==2

d_tau=200
tau_span=T_orb/2:d_tau:1*T_orb;
[~,N_tau]=size(tau_span);

dt_m=200
dt=10;

%m-line computation
[Ni_0_plus,~]=Compute_ci(tau,Rho_0,Ni_0,X_goal);
[Rho_m,Ni_m]=Propagate_ci(dt_m,tau,Rho_0,Ni_0_plus);

figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_m(:,1),Rho_m(:,2),Rho_m(:,3),'ok',Rho_m(1,1),Rho_m(1,2),Rho_m(1,3),'og',Rho_m(end,1),Rho_m(end,2),Rho_m(end,3),'or',x_obs,y_obs,0,'*b')
title('m-line & all possible trajectories')
grid on 

[N_pos_m,~]=size(Rho_m);

i_m=N_pos_m; flag_coll_m=1;%--> i_m=N_pos_m this also impose to test all the tof's to reach the goal point starting from Rho_0

i0_m=1;

i_tau_n=1; flag_coll_n=1;
 
i_win=0;



%while i_m >=1 && flag_coll_m==1 
for i_m=2:1:N_pos_m-1
        
        Rho_M=Rho_m(i_m,:)';
        Ni_M=Ni_m(i_m,:)';
        
        
          
         %while flag_coll_n==1 && i_tau_n <=N_tau 
         
         for i_tau_n=1:N_tau
             
             tau_n=tau_span(i_tau_n);
             
        [Ni_0_plus_n,~]=Compute_ci(tau_n,Rho_0,Ni_0,Rho_M);%Compute_ci(tau,Rho_0,Ni_0_minus,X_goal)

        
        [Rho_n_before,Ni_n_before]=Propagate_ci(dt,tau_n,Rho_0,Ni_0_plus_n);%--> trajectory BEFORE reincounter with  m-line
        [Rho_n_after,Ni_n_after]=Propagate_ci(dt,(tau+dt_m-i_m*dt_m),Rho_M,Ni_M);%--> trajectory AFTER reincounter with m-line
         
        Rho_n=[Rho_n_before;Rho_n_after];
        [flag_coll_n]=Evaluate_collision(Rho_n,Phi0_2d,step_grid);
        
        if flag_coll_n==0
            
            flag_coll_m=0;
            i_win=i_win+1;
            Rho_n_win=Rho_n;
            tau_n_win=tau_n;
            DV_win=norm((Ni_0_plus-Ni_0))+norm(Ni_n_after(1,:)-Ni_n_before(end,:));
            
            
        end
        plot3(Rho_n(:,1),Rho_n(:,2),Rho_n(:,3),'--b',Rho_n(1,1),Rho_n(1,2),Rho_n(1,3),'og',Rho_n(end,1),Rho_n(end,2),Rho_n(end,3),'or',x_obs,y_obs,0,'*b')



        %i_tau_n=i_tau_n+1;
         end

 plot3(Rho_n(:,1),Rho_n(:,2),Rho_n(:,3),'--b',Rho_n(1,1),Rho_n(1,2),Rho_n(1,3),'og',Rho_n(end,1),Rho_n(end,2),Rho_n(end,3),'or',x_obs,y_obs,0,'*b')

    
%i_m=i_m-1;  



if flag_coll_m==0
  
  
figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_n_win(:,1),Rho_n_win(:,2),Rho_n_win(:,3),'--b',Rho_n_win(1,1),Rho_n_win(1,2),Rho_n_win(1,3),'og',Rho_n_win(end,1),Rho_n_win(end,2),Rho_n_win(end,3),'or',x_obs,y_obs,0,'*b')
hold on
plot3(Rho_m(:,1),Rho_m(:,2),Rho_m(:,3),'ok',Rho_m(1,1),Rho_m(1,2),Rho_m(1,3),'og',Rho_m(end,1),Rho_m(end,2),Rho_m(end,3))
title('collision Avoidance pseudo BuG method & m-line ')
grid on 

disp('fine')
    
end 

end

elseif select_Alg==3
    
    

d_tau=500
tau_span=d_tau:d_tau:1*T_orb;
[~,N_tau]=size(tau_span)

tau_bug_0=tau;

tau_bug_trasl= tau_bug_0/8;



dt_m_bug_0=10;
dt_m=dt_m_bug_0*10
dt_bug_0=1;
dt=1;

[Ni_0_plus,~]=Compute_ci(tau,Rho_0,Ni_0,X_goal);
[Rho_m,Ni_m]=Propagate_ci(dt_m,tau,Rho_0,Ni_0_plus);


figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_m(:,1),Rho_m(:,2),Rho_m(:,3),'ok',Rho_m(1,1),Rho_m(1,2),Rho_m(1,3),'og',Rho_m(end,1),Rho_m(end,2),Rho_m(end,3),'or')
title('Initial (0) m-line ')
grid on 



flag_coll_mm=1;

R_gnc=15;

R_gnc_theoretical=0;

Rho_TRAJ=Rho_0';
Ni_TRAJ=Ni_0';


count_bug_3=0;
flag_bug_3=1;

while flag_bug_3==1 && count_bug_3 <=40
    
[Ni_0_plus,~]=Compute_ci(tau,Rho_0,Ni_0,X_goal);
[Rho_m,Ni_m]=Propagate_ci(dt_m,tau,Rho_0,Ni_0_plus);

  [N_m,~]=size(Rho_m);
  
  i_m=N_m;
  flag_coll_m=1;


  while i_m >=1 && flag_coll_m==1 
        
        Rho_M=Rho_m(i_m,:)';
        Ni_M=Ni_m(i_m,:)';
        
        flag_coll_n=1;
        i_tau_n=N_tau;
        
         while flag_coll_n==1 && i_tau_n >=1 
             
             tau_n=tau_span(i_tau_n);
             
        [Ni_0_plus_n,~]=Compute_ci(tau_n,Rho_0,Ni_0,Rho_M);
        
        [Rho_n_before,Ni_n_before]=Propagate_ci(dt,tau_n,Rho_0,Ni_0_plus_n);%--> trajectory BEFORE reincounter with  m-line
        [Rho_n_after,Ni_n_after]=Propagate_ci(dt,(tau+dt_m-i_m*dt_m),Rho_M,Ni_M);%--> trajectory AFTER reincounter with m-line
         
        Rho_n=[Rho_n_before;Rho_n_after];
        Ni_n=[Ni_n_before;Ni_n_after];
        [flag_coll_n]=Evaluate_collision(Rho_n,Phi0_2d,step_grid);
        
        if flag_coll_n==0
            
            flag_coll_m=0;
            flag_bug_3=0;
            Rho_traj=Rho_n;
            Ni_traj=Ni_n;
            tau_n_win=tau_n;
            DV_win=norm((Ni_0_plus-Ni_0))+norm(Ni_n_after(1,:)-Ni_n_before(end,:));
        end
        
        i_tau_n=i_tau_n-1;
         end
  i_m=i_m-1;      
  end
         
     if flag_coll_m==1
            
            [Ni_0_plus,~]=Compute_ci(tau_bug_0,Rho_0,Ni_0,X_goal);
            [Rho_m_bug_0,Ni_m_bug_0]=Propagate_ci(dt_bug_0,tau_bug_0,Rho_0,Ni_0_plus); %m-0 line
 
          %[Rho_m_star,O_k,K]=Tangent_bug(Rho_m_bug_0,Phi_2d_contour,step_grid,R_gnc);
           [Rho_m_star,O_k,K]=Tangent_bug_2Pt(Rho_m,Phi_2d_contour,step_grid,R_gnc);
         
         %[Rho_0_new]=Heuristic_tangent_smart_New(Rho_m_star,K,O_k,R_gnc,X_goal,Phi0_2d,step_grid);
         %[Rho_0_new]=Heuristic_tangent_smart_ortodox(Rho_m_star,Rho_0,K,O_k,R_gnc,X_goal,Phi0_2d,step_grid);
         [Rho_0_new,n_t]=Heuristic_tangent_smart_ortodox_2Pt(Rho_m_star,Rho_0,K,O_k,R_gnc,X_goal,Phi0_2d,step_grid);



          [Ni_0_plus,~]=Compute_ci(tau_bug_trasl,Rho_0,Ni_0,Rho_0_new);
          [Rho_traj,Ni_traj]=Propagate_ci(dt_bug_0,tau_bug_trasl,Rho_0,Ni_0_plus);
          
     end
     
     
     
     

Rho_TRAJ=[Rho_TRAJ;Rho_traj]; 
Ni_TRAJ=[ Ni_TRAJ;Ni_traj];
  
 Rho_0=Rho_TRAJ(end,:)';
 Ni_0=Ni_TRAJ(end,:)';
  
  
  
  
count_bug_3=count_bug_3+1;
end  



figure()
contour(R1_boundedx, R1_boundedy,Phi0_2d')
hold on
plot3(Rho_TRAJ(:,1),Rho_TRAJ(:,2),Rho_TRAJ(:,3),'--b',Rho_TRAJ(1,1),Rho_TRAJ(1,2),Rho_TRAJ(1,3),'og',Rho_TRAJ(end,1),Rho_TRAJ(end,2),Rho_TRAJ(end,3),'or',x_obs,y_obs,0,'*b')
title('trajectory ')
grid on 



%% close
end
