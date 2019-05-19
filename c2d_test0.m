clear all
%close all

%% ISS Orbit random chaser.

%orbit
Mu=(3.9856e+14);
r_Earth=6380*10^3;

h_orb=400*10^3;
i_orb=51.65*pi/180;
om_orb=10*pi/180;
OM_orb=20*pi/180;
theta0_orb=0;
a_orb=r_Earth+h_orb;
e_orb=0;


T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)



omega=sqrt(Mu/((r_Earth+h_orb)^3));
w=omega;
setGlobalW(w);

d_TC=200;


ROT_plane0=[cos(om_orb)*cos(OM_orb)-sin(om_orb)*cos(i_orb)*sin(OM_orb) -sin(om_orb)*cos(OM_orb)-cos(om_orb)*cos(i_orb)*sin(OM_orb) sin(i_orb)*sin(OM_orb)
            cos(om_orb)*sin(OM_orb)+sin(om_orb)*cos(i_orb)*cos(OM_orb) -sin(om_orb)*sin(OM_orb)+cos(om_orb)*cos(i_orb)*cos(OM_orb) -sin(i_orb)*cos(OM_orb)
              sin(om_orb)*sin(i_orb)                                                 cos(om_orb)*sin(i_orb)                                 cos(i_orb)]  ;

ROT_theta0=[ cos(theta0_orb), sin(theta0_orb),0;
           -sin(theta0_orb), cos(theta0_orb),0;
                 0       ,    0       ,1];
             
ROT_0=ROT_plane0*ROT_theta0;



%target:
[R_0t,V_0t]=OrbPar2RV_MY(a_orb,e_orb,i_orb,OM_orb,om_orb,theta0_orb,Mu);
w_v=cross(R_0t,V_0t)/(norm(R_0t)^2);

%chaser:
R_0c=(R_0t+rand(3,1)*0.25*1e2); 
V_0c=(V_0t-rand(3,1)*0.5*1e1);

% %target:REP
% R_0t=[(r_Earth+h_orb);0;0];
% V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];
% R_T=r_Earth+h_orb;
% V_T=norm(V_0t);
% % 
% %chaser:REP
% R_0c=(R_0t-[d_TC;0;0]); %-> pert. only in plane: circular orbit
% R_C=R_T-d_TC
% omega_c=sqrt(Mu/(norm(R_0c)^3));
% V_0c=[0;sqrt(Mu/norm(R_0c));0];
% V_0c=[0;omega*R_C;0];%-> to set zero velocity between chaser and target
% R_C=R_T-d_TC;
% V_C=norm(V_0c);


%Relative position:
Rrel_I0=R_0c-R_0t;
Vrel_I0=V_0c-V_0t;

Rho_0=ROT_0'*Rrel_I0;
Ni_0=ROT_0'*(Vrel_I0-cross(w_v,Rrel_I0));

%% Linearization & Digitalization 
%actuators:
m=5;
f_sat=2/5*1e-3*m;% ref1: a_sat=4*1e-6N && m_cubesat=3Kg 
                         % ref2[PLASMA]: f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3 http://pepl.engin.umich.edu/thrusters/CAT.html
                         % ref3[COMMERCIAL]: http://www.busek.com/cubesatprop__main.htm
                         % ref4[NANOSAT]:f_sat=100*1e-6 N &&  Imp_duration= 2*1e-3 s && Isp= 50:100 s http://www.cubesatshop.com/index.php?page=shop.product_details&flypage=flypage.tpl&product_id=74&vmcchk=1&option=com_virtuemart&Itemid=65f_sat=2*1e-3N/ && m_cubesat=5Kg--> a_sat=2/5*1e-3
                         
%sys:


A=[0,    0,   0,    1,   0,   0;
   0,    0,   0,    0,   1,   0;
   0,    0,   0,    0,   0,   1;
   3*w^2 0,   0,    0,  2*w,  0;
   0,    0,   0,  -2*w,  0,   0;
   0,    0,  -w^2,  0,   0,   0];

B=[0,    0,    0;
   0,    0,    0;
   0,    0,    0;
   1/m,  0,    0;
   0,    1/m,  0;
   0,    0,  1/m];

C=eye(6);

D=0;

sys=ss(A,B,C,D)

%digital-system:
Tc=1e-1;

[sysd]=c2d(sys,Tc,'zoh');%--> zero-hold derivation
[Ad,Bd,~,~]=ssdata(sysd) 

%% Comparison:  integrated CW  <--> x_k+1=Ad_k*x_k+Bd_k*u_k: 


%FREE-DRIFT:

%continous:
t_end_CW=1/0.5*T_orb;
x_0=[Rho_0;Ni_0];

options = odeset('RelTol',1e-13,'AbsTol', 1e-30);
[t_CW,x_CW]= ode113('Chol_Wilt_Hill',t_end_CW,x_0,options);%<--[integration of C-W eq. in Hill reference frame]

Rho_CW=[x_CW(:,1),x_CW(:,2),x_CW(:,3)];

%discrete:
T_d=0:Tc:t_end_CW;
[~,N_d]=size(T_d);

x_d=zeros(N_d,6);
x_d(1,:)=x_0';

for i=2:N_d
    x_d(i,:)=(Ad*x_d(i-1,:)')';
end

figure()
plot3(Rho_CW(:,1),Rho_CW(:,2),Rho_CW(:,3),'k',Rho_CW(1,1),Rho_CW(1,2),Rho_CW(1,3),'og',Rho_CW(end,1),Rho_CW(end,2),Rho_CW(end,3),'or',x_d(:,1),x_d(:,2),x_d(:,3),'--b')
title('Free Drift: (--b)<--> digital &&  (k)<-->continous ')
hold on
grid on

%FORCED-PROBLEM: R-bar sliding

%continous:

Ni_0_new=[-Ni_0(1)/2;0;0];

x_obs=+100;

V_slide=norm(Ni_0_new);
R_0=Rho_0(1);
Z_0=Rho_0(3);
t_slide=abs(norm(Rho_0(1)-x_obs)/V_slide);


w=omega;
setGlobalW(w);
setGlobalv_slide(V_slide )
setGlobalr_0(R_0 )
setGlobalz_0(Z_0 )

x_0_slide=[Rho_0;Ni_0_new];


options = odeset('RelTol',1e-13,'AbsTol', 1e-10);
[t_int,x_CW_slide]= ode113('Chol_Wilt_Hill_Rbar_slide_Zcontrol',t_slide,x_0_slide,options);%<--[R-bar continous thrust of C-W eq. in Hill reference frame]

Rho_CW_slide=[x_CW_slide(:,1),x_CW_slide(:,2),x_CW_slide(:,3)];


T_d_slide=0:Tc:t_slide;
[~,N_d_slide]=size(T_d_slide);

x_d_slide=zeros(N_d_slide,6);
x_d_slide(1,:)=x_0_slide';

u_d_slide=zeros(N_d_slide,3);
u_d_slide(1,:)=m*[-3*w^2*R_0,-2*w*V_slide,w^2*Z_0];
for i=2:N_d_slide
    t=T_d_slide(i);
    u_d_slide(i,:)=m*[-3*w^2*(R_0-V_slide*t),-2*w*V_slide, w^2*Z_0];
    x_d_slide(i,:)=(Ad*x_d_slide(i-1,:)')'+(Bd*u_d_slide(i-1,:)')';
end


figure()
plot3(Rho_CW_slide(:,1),Rho_CW_slide(:,2),Rho_CW_slide(:,3),'k',Rho_CW_slide(1,1),Rho_CW_slide(1,2),Rho_CW_slide(1,3),'og',Rho_CW_slide(end,1),Rho_CW_slide(end,2),Rho_CW_slide(end,3),'or',x_d_slide(:,1),x_d_slide(:,2),x_d_slide(:,3),'--b')
title('R-BAR slide: (--b)<--> digital  && k<-->continous ')
hold on
grid on

%% Deviation from linearity [C-W model]

 X_0t=[R_0t;V_0t];
  X_0c=[R_0c;V_0c];
   X_0t_c=[X_0t;X_0c];

options = odeset('RelTol',1e-13,'AbsTol', 1e-30);
[t_int,Xt_c]= ode113('r2BP_t_c',t_end_CW,X_0t_c,options);

R_t=[Xt_c(:,1),Xt_c(:,2),Xt_c(:,3)];
R_c=[Xt_c(:,7),Xt_c(:,8),Xt_c(:,9)];
V_t=[Xt_c(:,4),Xt_c(:,5),Xt_c(:,6)];
V_c=[Xt_c(:,10),Xt_c(:,11),Xt_c(:,12)];

Rrel_I=R_c-R_t; %--> free motion
Vrel_I=V_c-V_t;

% LVLH--> I for digital solution:>> ?main" translation

Rrel_I_d=zeros(N_d,3);
Vrel_I_d=zeros(N_d,3);

for i=1:N_d
    t=T_d(i);
    theta=theta0_orb-w*t;%!! TAKE CARE FOR THIS SIGN !! 
    ROT_theta=[ cos(theta), sin(theta),    0;
               -sin(theta), cos(theta),    0;
                 0       ,    0       ,    1];
    ROT=ROT_plane0*ROT_theta;
    
    Rho_d=x_d(i,1:3);
    Ni_d=x_d(i,4:6);
    
    Rrel_I_d(i,:)=(ROT*Rho_d')';
    Vrel_I_d(i,:)=(ROT*(Ni_d'+cross([0;0;omega],Rho_d')))';%!! TAKE CARE FOR THIS CROSS-PRODUCT !!
    
end

figure()
plot3(Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'--b',Rrel_I_d(:,1),Rrel_I_d(:,2),Rrel_I_d(:,3),'--r')
title('(--b) Rrel_I ( kepler_integration) VS (--r) Rrel_I_d (d: rotated digital solution) IN r.f.')
hold on
grid on 

figure()
plot3(Vrel_I(:,1),Vrel_I(:,2),Vrel_I(:,3),'--b',Vrel_I_d(:,1),Vrel_I_d(:,2),Vrel_I_d(:,3),'--r')
title('(--b) Vrel_I  (kepler-integration) VS (--r) (d: rotated digital solution) IN r.f.')
hold on
grid on 


% I--> LVLH for digital solution:>> LQR controlled motion
[N_int,~]=size(t_int)
Rrel_LVLH_int=zeros(N_int,3);
Vrel_LVLH_int=zeros(N_int,3);

i_v=[1;0;0];
j_v=[0;1;0]

for i=1:N_int
    
    
    t=t_int(i);
    
    theta_clock=theta0_orb-w*t%!! TAKE CARE FOR THIS SIGN !!
    
    r_t=R_t(i,:)';
    v_t=V_t(i,:)';
    
    %[~,~,~,~,~,theta]=RV2OrbPar_My(r_t,v_t,Mu);
    RT=norm(R_t(i,:));
    theta=-acos(dot(i_v,(R_t(i,:)/RT)));
     if dot(j_v,(R_t(i,:)/RT)) <0 
         
         theta=2*pi-theta
         
     
     end

    
    ROT_theta=[ cos(theta), sin(theta),    0;
               -sin(theta), cos(theta),    0;
                 0       ,    0       ,    1];


    ROT=ROT_plane0*ROT_theta;
    
    Rrel_LVLH_int(i,:)=(ROT'*(R_c(i,:)-R_t(i,:))')';
    
   
   
   Vrel_LVLH_int(i,:)=(ROT'*((V_c(i,:)-V_t(i,:))'-cross(w_v,(R_c(i,:)-R_t(i,:))')))';%!! TAKE CARE FOR THIS CROSS-PRODUCT !!
   
    
end

figure()
plot3(Rrel_LVLH_int(:,1),Rrel_LVLH_int(:,2),Rrel_LVLH_int(:,3),'--b',x_d(:,1),x_d(:,2),x_d(:,3),'--r')
title('(--b) Rrel_{LVLH{_int}} ( kepler-integration) VS (--r) Rho_d (d: rotated digital solution) LVLH r.f.')
hold on
grid on 

figure()
plot3(Vrel_LVLH_int(:,1),Vrel_LVLH_int(:,2),Vrel_LVLH_int(:,3),'--b',x_d(:,4),x_d(:,5),x_d(:,6),'--r')
title('(--b) Vrel_{LVLH_{int}} (kepler-integration) VS (--r) Ni_d (d: rotated digital solution) LVLH r.f.')
hold on
grid on 

