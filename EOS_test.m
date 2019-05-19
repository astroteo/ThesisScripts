%% EOS (Ellips Of Safety) test 
close all
clear all
%orbit setup.

Mu=(3.9856e+14);
r_Earth=6380*10^3;
h_orb=250*10^3;

T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)

t_end_CW=1/0.5*T_orb
d_TC=200
omega=sqrt(Mu/((r_Earth+h_orb)^3))

%target
R_0t=[(r_Earth+h_orb);0;0];
V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];
R_T=r_Earth+h_orb;
V_T=norm(V_0t);



%% EOS tranfer ( analytical)

x_obs=0
y_obs=50
z_obs=10

X_obs=[x_obs,y_obs,z_obs];

%push-->S1
Rho_0=[-100;0;0];
Ni_0=[-0.1;0;0];%<-- PUSH
DV_PUSH=norm(Ni_0);
w=omega;
t=0;

t_starP_S1=(1/w)*abs(atan(-Ni_0(1)/(3*Rho_0(1)*w+2*Ni_0(2))))
t_starP_S1_APP=2*pi/(w*4)
N_stepP_S1=1000
timeP_S1=linspace(0,t_starP_S1,N_stepP_S1);
% timeP_S1=linspace(0,t_starP_S1_APP,N_stepP_S1);
xP_S1=zeros(N_stepP_S1,3);

for i=1:N_stepP_S1
t=timeP_S1(i);
xP_S1(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*t)+Ni_0(1)/w*sin(w*t)+4*Rho_0(1)+2/w*Ni_0(2);
xP_S1(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*t)-2/w*Ni_0(1)*cos(w*t)+(6*w*Rho_0(1)+3*Ni_0(2))*t-Rho_0(2)+2/w*Ni_0(1));
xP_S1(i,3)=Rho_0(3)*cos(w*t)+Ni_0(3)/w*sin(w*t);
end

x_dotP_S1(i,1)=(3*Rho_0(1)*w+2*Ni_0(2))*sin(w*t)+Ni_0(1)*cos(w*t);
x_dotP_S1(i,2)=-((-4*Ni_0(2)-6*w*Rho_0(1))*cos(w*t)+2*Ni_0(1)*sin(w*t)+6*w*Rho_0(1)+3*Ni_0(2) );
x_dotP_S1(i,3)=Ni_0(3)*cos(w*t)-Rho_0(3)*w*sin(w*t);

Ni_S1_minus=[x_dotP_S1(end,1);x_dotP_S1(end,2);x_dotP_S1(end,3)]
Rho_S1=[xP_S1(end,1);xP_S1(end,2);xP_S1(end,3)];
    

timeP_S1_PROP=linspace(0,200*t_starP_S1,N_stepP_S1);
xP_S1_PROP=zeros(N_stepP_S1,3);

for i=1:N_stepP_S1
t=timeP_S1_PROP(i);
xP_S1_PROP(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*t)+Ni_0(1)/w*sin(w*t)+4*Rho_0(1)+2/w*Ni_0(2);
xP_S1_PROP(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*t)-2/w*Ni_0(1)*cos(w*t)+(6*w*Rho_0(1)+3*Ni_0(2))*t-Rho_0(2)+2/w*Ni_0(1));
xP_S1_PROP(i,3)=Rho_0(3)*cos(w*t)+Ni_0(3)/w*sin(w*t);
end

%S1-->S2
tau=(2*pi/4)/w%;<-- EOS parameters
a=200
z_max=100


y_S2=(y_obs)+(a/z_max)*(z_max+z_obs);
x_S2=0;
z_S2=z_max;

        

x_S1=Rho_S1(1);
y_S1=Rho_S1(2);
z_S1=Rho_S1(3);


Ct=cos(w*tau);
St=sin(w*tau);



%[MATHEMATICA SOLUTION (matrix)]
det=(4*St)/(w^3) -(8*Ct*St)/(w^3) +(4*Ct^2*St)/(w^3)+(4*St^3)/(w^3) -(3*St^2*tau)/(w^2);

N_tau_inv=1/det*[(4*St^2)/(w^2)-(3*St*tau)/w,     -((2*St)/(w^2))+(2*Ct*St)/(w^2),                        0;
                 (2*St)/(w^2)-(2*Ct*St)/(w^2),              St^2/(w^2),                                   0;
                           0,                                  0,              4/(w^2)-(8*Ct)/(w^2)+(4*Ct^2)/(w^2)+(4*St^2)/(w^2)-(3*St*tau)/w];
                   
M_tau=[-3*Ct+4,        0,   0;
        6*St-6*w*tau,  1,   0;                   
               0    ,  0,  Ct];
Rho_S2_design=[x_S2;
               y_S2;
               z_S2];                 

           Ni_S1_plus=N_tau_inv*(Rho_S2_design-M_tau*Rho_S1);                 
           DV_S1=norm(Ni_S1_plus-Ni_S1_minus);
           
Ni_0=Ni_S1_plus;
Rho_0=[x_S1,y_S1,z_S1];
N_stepS1_S2=ceil(N_stepP_S1/10*(tau/t_starP_S1));
timeS1_S2=linspace(0,1*tau,N_stepS1_S2);
xS1_S2=zeros(N_stepS1_S2,3);
          
xS1_S2_MATH=zeros(N_stepS1_S2,3);           
for i=1:N_stepS1_S2
t=timeS1_S2(i);
xS1_S2_MATH(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*t)+Ni_0(1)/w*sin(w*t)+4*Rho_0(1)+2/w*Ni_0(2);
xS1_S2_MATH(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*t)-2/w*Ni_0(1)*cos(w*t)+(6*w*Rho_0(1)+3*Ni_0(2))*t-Rho_0(2)+2/w*Ni_0(1));
xS1_S2_MATH(i,3)=Rho_0(3)*cos(w*t)+Ni_0(3)/w*sin(w*t);
end        
 

x_dotS1_S2_MATH(i,1)=(3*Rho_0(1)*w+2*Ni_0(2))*sin(w*t)+Ni_0(1)*cos(w*t);
x_dotS1_S2_MATH(i,2)=-((-4*Ni_0(2)-6*w*Rho_0(1))*cos(w*t)+2*Ni_0(1)*sin(w*t)+6*w*Rho_0(1)+3*Ni_0(2) );
x_dotS1_S2_MATH(i,3)=Ni_0(3)*cos(w*t)-Rho_0(3)*sin(w*t);

Ni_S2_MATH=[x_dotS1_S2_MATH(end,1);x_dotS1_S2_MATH(end,2);x_dotS1_S2_MATH(end,3)]
Rho_S2_MATH=[xS1_S2_MATH(end,1);xS1_S2_MATH(end,2);xS1_S2_MATH(end,3)]

timeS1_S2_PROP=linspace(0,16*tau,N_stepS1_S2); %<--orbit propagation for 1 period
xS1_S2_PROP=zeros(N_stepS1_S2,3);
for i=1:N_stepS1_S2
tt=timeS1_S2_PROP(i);
xS1_S2_PROP(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*tt)+Ni_0(1)/w*sin(w*tt)+4*Rho_0(1)+2/w*Ni_0(2);
xS1_S2_PROP(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*tt)-2/w*Ni_0(1)*cos(w*tt)+(6*w*Rho_0(1)+3*Ni_0(2))*tt-Rho_0(2)+2/w*Ni_0(1));
xS1_S2_PROP(i,3)=Rho_0(3)*cos(w*tt)+Ni_0(3)/w*sin(w*tt);
end        
 

%S2-->S3 [ trajectory propagation ]
time0_S2_S3=t;
N_stepS2_S3=1000*N_stepS1_S2;
timeS2_S3=linspace(0,2*tau,N_stepS2_S3);%<-- 1/2 of a period to reach the opposite point

xS2_S3_MATH=zeros(N_stepS2_S3,3);
flag_y=0;
flag_z=1;
flag_S3=0;

i=1;
while i <= N_stepS2_S3 && flag_S3==0 
    
t=time0_S2_S3+timeS2_S3(i);
xS2_S3_MATH(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*t)+Ni_0(1)/w*sin(w*t)+4*Rho_0(1)+2/w*Ni_0(2);
xS2_S3_MATH(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*t)-2/w*Ni_0(1)*cos(w*t)+(6*w*Rho_0(1)+3*Ni_0(2))*t-Rho_0(2)+2/w*Ni_0(1));
xS2_S3_MATH(i,3)=Rho_0(3)*cos(w*t)+Ni_0(3)/w*sin(w*t);

   if xS2_S3_MATH(i,2) <= y_obs && flag_y==0;
       
       i_special_y=i;
       flag_y=1;
   end
   
   
    if xS2_S3_MATH(i,3) <= z_obs && flag_z==0;
       
       i_special_z=i;
       flag_z=1;
    end
       
   if flag_z >0  && flag_y > 0
        
        flag_S3=1;
        i_special_S3=i;
        
    end
    
  i=i+1;     
end     

if flag_S3==0;
    
    i_special_S3=N_stepS2_S3;
end

x_dotS2_S3_MATH(i,1)=(3*Rho_0(1)*w+2*Ni_0(2))*sin(w*t)+Ni_0(1)*cos(w*t);
x_dotS2_S3_MATH(i,2)=-((-4*Ni_0(2)-6*w*Rho_0(1))*cos(w*t)+2*Ni_0(1)*sin(w*t)+6*w*Rho_0(1)+3*Ni_0(2) );
x_dotS2_S3_MATH(i,3)=Ni_0(3)*cos(w*t)-Rho_0(3)*sin(w*t);
Ni_S3_minus=[x_dotS2_S3_MATH(end,1);x_dotS2_S3_MATH(end,2);x_dotS2_S3_MATH(end,3)]
Rho_S3=[xS2_S3_MATH(i_special_S3,1);xS2_S3_MATH(i_special_S3,2);xS2_S3_MATH(i_special_S3,3)]

% S3-->OBS [R-bar sliding]
V_slide=Ni_S3_minus(1)/2;
Ni_S3_plus=[-V_slide;0;0]
DV_S3=norm(Ni_S3_minus-Ni_S3_plus);%<-- Brake at S3

V_slide=norm(Ni_S3_plus);
R_0=Rho_S3(1);
Z_0=Rho_S3(3);
Dt_slide=abs(norm(Rho_S3(1)-x_obs)/V_slide);


w=omega;
setGlobalW(w);
setGlobalv_slide(V_slide )
setGlobalr_0(R_0 )
setGlobalz_0(Z_0 )


options = odeset('RelTol',1e-13,'AbsTol', 1e-10);%<--[integration of C-W eq. in Hill reference frame]
[t_int,x_vS3_OBS]= ode113('Chol_Wilt_Hill_Rbar_slide_Zcontrol',Dt_slide, [Rho_S3(1);Rho_S3(2);Rho_S3(3);Ni_S3_plus(1);Ni_S3_plus(2);Ni_S3_plus(3)],options);
%[t_int,x_vS3_OBS]= ode113('Chol_Wilt_Hill_Rbar_slide',Dt_slide, [Rho_S3(1);Rho_S3(2);Rho_S3(3);Ni_S3_plus(1);Ni_S3_plus(2);Ni_S3_plus(3)],options);

xS3_OBS=[x_vS3_OBS(:,1),x_vS3_OBS(:,2),x_vS3_OBS(:,3)];
vS3_OBS=[x_vS3_OBS(:,4),x_vS3_OBS(:,5),x_vS3_OBS(:,6)];

DV_final=norm(vS3_OBS(end,:));%<-- Brake at S3




%results:
disp('PUSH:')
DV_PUSH
disp('S1, impulse: insertion into EOS ')
DV_S1
disp('S3,impulse: brake to reach v_slide*x')
DV_S3
disp('S3 to OBS z control') %<-- continous thrust phase
DV_z_control=3*w^2*Z_0*Dt_slide
disp('S3 to OBS y control') 
DV_y_control=2*w*V_slide*Dt_slide;
disp('S3 to OBS x control') 
DV_x_control=abs(3*w^2*(R_0+(1/2)*V_slide*t)*t);
disp('FINAL BRAKE (zero final velocity alon "x"')
DV_final
DV_total=DV_final+DV_x_control+DV_y_control+DV_z_control+DV_S3+DV_S1+DV_PUSH

 disp(' Design final postion:')
x_obs=X_obs(1)
y_obs=X_obs(2)
z_obs=X_obs(3)

 disp('Real final postion:')
x_OBS=xS3_OBS(end,1)
y_OBS=xS3_OBS(end,2)
z_OBS=xS3_OBS(end,3)



 %plots:
[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,100,50,30);
figure()
plot3(xP_S1(:,1),xP_S1(:,2),xP_S1(:,3),'--k',xP_S1(1,1),xP_S1(1,2),xP_S1(1,3),'ok',xP_S1(end,1),xP_S1(end,2),xP_S1(end,3),'*k')
title('PUSH-->S1')
hold on
surf(x_ISS, y_ISS, z_ISS,'facecolor','g')
axis equal
grid on

figure()
plot3(xP_S1_PROP(:,1),xP_S1_PROP(:,2),xP_S1_PROP(:,3),'--k',xP_S1_PROP(1,1),xP_S1_PROP(1,2),xP_S1_PROP(1,3),'ok',xP_S1_PROP(end,1),xP_S1_PROP(end,2),xP_S1_PROP(end,3),'*k')
title('PUSH-->S1: PROPAGATION')
axis equal
grid on



figure()
plot3(xP_S1(:,1),xP_S1(:,2),xP_S1(:,3),'--k',xP_S1(1,1),xP_S1(1,2),xP_S1(1,3),'ok',xP_S1(end,1),xP_S1(end,2),xP_S1(end,3),'*k',xS1_S2_MATH(:,1),xS1_S2_MATH(:,2),xS1_S2_MATH(:,3),'--b',xS1_S2_MATH(1,1),xS1_S2_MATH(1,2),xS1_S2_MATH(1,3),'ob',xS1_S2_MATH(end,1),xS1_S2_MATH(end,2),xS1_S2_MATH(end,3),'*b' )
title('PUSH-->S1-->S2')
hold on
surf(x_ISS, y_ISS, z_ISS,'facecolor','g')
axis equal
hold on
grid on


figure()
plot3(xP_S1(:,1),xP_S1(:,2),xP_S1(:,3),'--k',xP_S1(1,1),xP_S1(1,2),xP_S1(1,3),'ok',xP_S1(end,1),xP_S1(end,2),xP_S1(end,3),'*k',xS1_S2_MATH(:,1),xS1_S2_MATH(:,2),xS1_S2_MATH(:,3),'--b',xS1_S2_MATH(1,1),xS1_S2_MATH(1,2),xS1_S2_MATH(1,3),'ob',xS1_S2_MATH(end,1),xS1_S2_MATH(end,2),xS1_S2_MATH(end,3),'*b',xS2_S3_MATH(1:i_special_S3,1),xS2_S3_MATH(1:i_special_S3,2),xS2_S3_MATH(1:i_special_S3,3),'b',xS2_S3_MATH(1,1),xS2_S3_MATH(1,2),xS2_S3_MATH(1,3),'ob',xS2_S3_MATH(i_special_S3,1),xS2_S3_MATH(i_special_S3,2),xS2_S3_MATH(i_special_S3,3),'*b' )
title('PUSH-->S2-->S3')
hold on
surf(x_ISS, y_ISS, z_ISS,'facecolor','g')
axis equal
hold on
grid on

figure()
plot(xP_S1(:,2),xP_S1(:,3),'--k',xP_S1(1,2),xP_S1(1,3),'ok',xP_S1(end,2),xP_S1(end,3),'*k',xS1_S2_MATH(:,2),xS1_S2_MATH(:,3),'--b',xS1_S2_MATH(1,2),xS1_S2_MATH(1,3),'ob',xS1_S2_MATH(end,2),xS1_S2_MATH(end,3),'*b',xS2_S3_MATH(1:i_special_S3,2),xS2_S3_MATH(1:i_special_S3,3),'b',xS2_S3_MATH(1,2),xS2_S3_MATH(1,3),'ob',xS2_S3_MATH(i_special_S3,2),xS2_S3_MATH(i_special_S3,3),'*b' )
title('y-z plane')
axis equal
hold on
grid on

figure()
plot(xP_S1(:,1),xP_S1(:,2),'--k',xP_S1(1,1),xP_S1(1,2),'ok',xP_S1(end,1),xP_S1(end,2),'*k',xS1_S2_MATH(:,1),xS1_S2_MATH(:,2),'--b',xS1_S2_MATH(1,1),xS1_S2_MATH(1,2),'ob',xS1_S2_MATH(end,1),xS1_S2_MATH(end,2),'*b',xS2_S3_MATH(1:i_special_S3,1),xS2_S3_MATH(1:i_special_S3,2),'b',xS2_S3_MATH(1,1),xS2_S3_MATH(1,2),'ob',xS2_S3_MATH(i_special_S3,1),xS2_S3_MATH(i_special_S3,2),'*b',xS3_OBS(:,1),xS3_OBS(:,2),'r' )
title('x-y plane')
axis equal
hold on
grid on

figure()
plot3(xP_S1(:,1),xP_S1(:,2),xP_S1(:,3),'--k',xP_S1(1,1),xP_S1(1,2),xP_S1(1,3),'ok',xP_S1(end,1),xP_S1(end,2),xP_S1(end,3),'*k',xS1_S2_PROP(:,1),xS1_S2_PROP(:,2),xS1_S2_PROP(:,3),'--b',xS1_S2_PROP(1,1),xS1_S2_PROP(1,2),xS1_S2_PROP(1,3),'ob',xS1_S2_PROP(end,1),xS1_S2_PROP(end,2),xS1_S2_PROP(end,3),'*b')
title('S2-->S3:(EOS:PROPAGATION)')
hold on
surf(x_ISS, y_ISS, z_ISS,'facecolor','g')
axis equal
hold on
grid on

figure()
plot3(xP_S1(:,1),xP_S1(:,2),xP_S1(:,3),'--k',xP_S1(1,1),xP_S1(1,2),xP_S1(1,3),'ok',xP_S1(end,1),xP_S1(end,2),xP_S1(end,3),'*k',xS1_S2_MATH(:,1),xS1_S2_MATH(:,2),xS1_S2_MATH(:,3),'--b',xS1_S2_MATH(1,1),xS1_S2_MATH(1,2),xS1_S2_MATH(1,3),'ob',xS1_S2_MATH(end,1),xS1_S2_MATH(end,2),xS1_S2_MATH(end,3),'*b',xS2_S3_MATH(1:i_special_S3,1),xS2_S3_MATH(1:i_special_S3,2),xS2_S3_MATH(1:i_special_S3,3),'b',xS2_S3_MATH(1,1),xS2_S3_MATH(1,2),xS2_S3_MATH(1,3),'ob',xS2_S3_MATH(i_special_S3,1),xS2_S3_MATH(i_special_S3,2),xS2_S3_MATH(i_special_S3,3),'*b',xS3_OBS(:,1),xS3_OBS(:,2),xS3_OBS(:,3),'r' )
title('GLOBAL transfer')
hold on
surf(x_ISS, y_ISS, z_ISS,'facecolor','g')
axis equal
hold on
grid on




%% TRASH



% S1--S2 [MY SOLUTION]
% F2=(x_S1-x_S2)/(Ct-1);
% F1=3*w*tau*St+8*(Ct-1);
% 
% Ni_S1_plus(1,1)=-w/F1*(Ct-1)*(2*(y_S1-y_S2)-4*F2*St-3*w*tau*(F2-x_S1));
% Ni_S1_plus(2,1)=-w/2*(3*x_S1+F2-Ni_S1_plus(1,1)/w*St/(Ct-1));
% Ni_S1_plus(3,1)=w/St*(z_S2-z_S1*Ct);
% 
% DV_1=norm(Ni_S1_plus-Ni_S1_minus);
% 
% 
% 
% Ni_0=Ni_S1_plus;
% Rho_0=[x_S1,y_S1,z_S1];
% N_stepS1_S2=ceil(N_stepP_S1*(tau/t_starP_S1));
% timeS1_S2=linspace(0,tau,N_stepS1_S2);
% xS1_S2=zeros(N_stepS1_S2,3);
% 
% for i=1:N_stepS1_S2
% t=timeS1_S2(i);
% xS1_S2(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*t)+Ni_0(1)/w*sin(w*t)+4*Rho_0(1)+2/w*Ni_0(2);
% xS1_S2(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*t)-2/w*Ni_0(1)*cos(w*t)+(6*w*Rho_0(1)+3*Ni_0(2))*t-Rho_0(2)+2/w*Ni_0(1));
% xS1_S2(i,3)=Rho_0(3)*cos(w*t)-Ni_0(3)*w*sin(w*t);
% end
% 
% x_dotS1_S2(i,1)=(3*Rho_0(1)*w+2*Ni_0(2))*sin(w*t)+Ni_0(1)*cos(w*t);
% x_dotS1_S2(i,2)=-((-4*Ni_0(2)-6*w*Rho_0(1))*cos(w*t)+2*Ni_0(1)*sin(w*t)+6*w*Rho_0(1)+3*Ni_0(2) );
% x_dotS1_S2(i,3)=Ni_0(3)*cos(w*t)-Rho_0(3)*sin(w*t);
% 
% 
% Ni_S2=[x_dotS1_S2(end,1);x_dotS1_S2(end,2);x_dotS1_S2(end,3)]
% Rho_S2=[xS1_S2(end,1);xS1_S2(end,2);xS1_S2(end,3)]
% 
% figure()
% plot3(xP_S1(:,1),xP_S1(:,2),xP_S1(:,3),'--k',xP_S1(1,1),xP_S1(1,2),xP_S1(1,3),'ok',xP_S1(end,1),xP_S1(end,2),xP_S1(end,3),'*k',xS1_S2(:,1),xS1_S2(:,2),xS1_S2(:,3),'--b',xS1_S2(1,1),xS1_S2(1,2),xS1_S2(1,3),'ob',xS1_S2(end,1),xS1_S2(end,2),xS1_S2(end,3),'*b' )
% hold on
% grid on



% while i <= N_stepS1_S2 && flag_S3==0 
% t=time0_S2_S3+timeS2_S3(i);
% xS2_S3_MATH(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*t)+Ni_0(1)/w*sin(w*t)+4*Rho_0(1)+2/w*Ni_0(2);
% xS2_S3_MATH(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*t)-2/w*Ni_0(1)*cos(w*t)+(6*w*Rho_0(1)+3*Ni_0(2))*t-Rho_0(2)+2/w*Ni_0(1));
% xS2_S3_MATH(i,3)=Rho_0(3)*cos(w*t)+Ni_0(3)/w*sin(w*t);
% 
%    if xS2_S3_MATH(i,2) <= y_obs && flag_y==0;
%        
%        i_special_y=i;
%        flag_y=1;
%    end
%    
%    
%     if xS2_S3_MATH(i,3) <= z_obs && flag_z==0;
%        
%        i_special_z=i;
%        flag_z=1;
%     end
%        
%    if flag_z >0  || flag_y > 0
%         
%         flag_S3=1;
%     end
%     
%        
% end  

% for i=1:N_stepS2_S3
% t=time0_S2_S3+timeS2_S3(i);
% xS2_S3_MATH(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*t)+Ni_0(1)/w*sin(w*t)+4*Rho_0(1)+2/w*Ni_0(2);
% xS2_S3_MATH(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*t)-2/w*Ni_0(1)*cos(w*t)+(6*w*Rho_0(1)+3*Ni_0(2))*t-Rho_0(2)+2/w*Ni_0(1));
% xS2_S3_MATH(i,3)=Rho_0(3)*cos(w*t)+Ni_0(3)/w*sin(w*t);
% 
% end        
