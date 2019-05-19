clc
clear all
close all

%% ISS Orbit ( target)

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


%% Desired Obsrvation Point & ISS Ellipsoid
x_obs=100+3
y_obs=200-3
z_obs=0



lenght_ISS=100; %--> along y
width_ISS=50; %--> alng z
height_ISS=30;


[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,lenght_ISS,height_ISS,width_ISS);





%% Laplace 2D Potential 

X_obs_2d=[x_obs;y_obs];%-> observation point given in LVLH r.f.

m1=lenght_ISS;
m2=width_ISS;


M_2d=[m1, 0;
      0, m2];

flagx_obs=0;
flagy_obs=0;

L_max=3200;
step_grid=10;
R1_bounded=-L_max:step_grid:(L_max-step_grid);

R1_boundedx=R1_bounded;
R1_boundedy=R1_bounded;

[~,N]=size(R1_boundedx);
[~,M]=size(R1_boundedy);

phi_2d=zeros(N,M);
phi_harm_2d=zeros(N,M);

setGlobaltotal_stepx(N);
setGlobaltotal_stepy(M); 


for i=1:N
    
    for j=1:M
        
       
            
            x=R1_boundedx(1,i);
            y=R1_boundedy(1,j);
           
            X_2d= [x;y];
            
           phi_2d(i,j)=1/2*(X_2d-X_obs_2d)'*M_2d*(X_2d-X_obs_2d);
           phi_harm_2d(i,j)=exp((X_2d-X_obs_2d)'*M*(X_2d-X_obs_2d));
          
            
            if y >= y_obs+step_grid && flagy_obs==0 
                
                j_obs=j-2;
                flagy_obs=1;
            end
            
            
            if x >=x_obs+step_grid && flagx_obs==0
                
                i_obs=i-2;
                flagx_obs=1;
            end
            
            
   
    end
end


%test to verify that the generated potential has it's own minima in X_obs

   x_obs_test= R1_bounded(i_obs)

   y_obs_test= R1_bounded(j_obs)


%check for minimum of laplace potential [COINCIDENT WITH (x_obs, y_obs)]
[phi_min,j_min] = min((min(phi_2d)));
[~,i_min]=(min(phi_2d));
i_min=i_min(1)
j_min
i_obs
j_obs

min_phi_2d=phi_2d(j_obs,i_obs)%<- expected to be equal to 0

step_gridx=step_grid;
step_gridy=step_grid;
[px,py] = gradient(phi_2d,step_gridy,step_gridx);

%test for discrete gradient formula: ODE discrete potential
x=x_obs;
y=y_obs;

total_stepx=N;
total_stepy=M;

if x >=0
i_bottom=total_stepx/2+1+floor(x/step_grid);%<--need to use ceil and floor since grad must be computed at the knots.
i_up=total_stepx/2+1+ceil(x/step_grid);

xd=(x-(i_bottom-N/2-1)*step_grid);



else     
i_up=ceil((total_stepx/2*step_grid-abs(x))/step_grid);
i_bottom=floor((total_stepx/2*step_grid-abs(x))/step_grid);

xd=x+(total_stepx/2-i_bottom)*step_grid;



end


if y >=0
j_bottom=total_stepx/2+1+floor(y/step_grid);
j_up=total_stepx/2+1+ceil(y/step_grid);


yd=(y-(j_bottom-M/2-1)*step_grid);


    
else     
j_up=ceil((total_stepy/2*step_grid-abs(y))/step_grid);
j_bottom=floor((total_stepy/2*step_grid-abs(y))/step_grid);

yd=y+(total_stepy/2-j_bottom)*step_grid;



end


Fx=-phi_2d(i_bottom,j_bottom)*(1-yd)  + phi_2d(i_up,j_bottom)*(1-yd) - phi_2d(i_bottom,j_up)*yd     + phi_2d(i_up,j_up)*yd;
Fy=-phi_2d(i_bottom,j_bottom)*(1-xd)  - phi_2d(i_up,j_bottom)*xd     + phi_2d(i_bottom,j_up)*(1-xd) + phi_2d(i_up,j_up)*xd;

i_bottom
i_up
j_bottom
j_up

yd
xd

%phi_c= phi_2d(i_bottom,j_bottom)*(1-xd)*(1-yd) + phi_2d(i_up,j_bottom)*xd*(1-yd) + phi_2d(i_bottom,j_up)*(1-xd)*yd + phi_2d(i_up,j_up)*xd*yd;
%f(x,y)= f(0,0) (1-x)(1-y) + f(1,0) x(1-y) + f(0,1) (1-x)y + f(1,1) xy. 



grad_mod_MY=norm([Fx,Fy])
grad_dir_MY=[Fx/grad_mod_MY;Fy/grad_mod_MY ]

Px=px(i_obs+1,j_obs+1)
Py=py(i_obs+1,j_obs+1)

grad_mod_MAT=norm([Px,Py])
grad_dir_MAT=[Px/grad_mod_MAT;Py/grad_mod_MAT ]


figure()
surf(x_ISS, y_ISS, z_ISS,'facecolor','g')
hold on
surf(R1_bounded,R1_bounded,phi_2d)
grid on  

figure()
contour(R1_bounded,R1_bounded,phi_2d)
hold on
surf(x_ISS, y_ISS,z_ISS,'facecolor','g')
hold on
plot3(x_obs,y_obs,0,'*b')
hold on
quiver(R1_bounded,R1_bounded,px,py)
hold off


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
Ni_0_2d=Ni_0_2d_CT-Ni_0_2d_RAND

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

%% GNC by Laplace Aritificial Potential

setGlobalstep_grid(step_grid);
setGlobalphi_2d(phi_2d );
setGlobalx_obs_2d(X_obs_2d);



% alpha IsSue
alpha_v=norm(Rho_0_2d-X_obs_2d)/t_end
% alpha_v=norm(Ni_0_2d)*1e-4;
setGlobalalpha_v(alpha_v );
t_end_Lap=t_end+0.1*T_orb


%integration
options = odeset('RelTol',1e-13,'AbsTol', 1e-3);
[t_CW,x_CW]= ode113('Chol_Wilt_Hill_Laplace_2d',t_end_Lap,X_0_2d,options);%<--[]

[CW,~]=size(t_CW);
zebra=zeros(CW,1);


%plot
figure()
surf(x_ISS, y_ISS,zeros(size(x_ISS)),'facecolor','g')
hold on
plot3(x_CW(:,1),x_CW(:,2),zebra,'--b',x_CW(1,1),x_CW(1,2),zebra,'og',x_CW(end,1),x_CW(end,2),zebra,'or',x_obs,y_obs,0,'*b')
grid on 





