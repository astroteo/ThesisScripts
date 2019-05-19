%[ TEST 0] trajectory tracking 
clear all
close all

%% 1) relative motion can be ignored from orbital POV

Mu=(3.9856e+14);
r_Earth=6380*10^3;
h_orb=400*10^3;

R_0t=[(r_Earth+h_orb);0;0];
V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];

[~,~,i_0,OM_0,om_0,theta_0]=RV2OrbPar_My(R_0t,V_0t,Mu);
[R_0_LVLH,V_0_LVLH]=In2LVLH(R_0t,V_0t,Mu);
X_0t=[R_0t;V_0t];
T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)

% R_0c=(R_0t+[rand(1,1);0;0]*1e2);
% V_0c=(V_0t+[0;rand(1,1);0]*1e1);

R_0c=(R_0t+rand(3,1)*1e0);
V_0c=(V_0t+rand(3,1)*1e-2);



Rrel_0_I=R_0c-R_0t

X_0c=[R_0c;V_0c]; %-> chaser [SPHERE] is assumed to move up to 1 dm/s in difference from the target at a maximum distnce of 10m
                                                        
N=1e4;
t_end=1/2*T_orb
tspan=linspace(0,t_end,N);

X_0t_c=[X_0t;X_0c]

options = odeset('RelTol',1e-14,'AbsTol', 1e-15);
[t_t,Xt_c]= ode113('r2BP_t_c',t_end,X_0t_c,options);

 
R_t=[Xt_c(:,1),Xt_c(:,2),Xt_c(:,3)];
R_c=[Xt_c(:,7),Xt_c(:,8),Xt_c(:,9)];


V_t=[Xt_c(:,4),Xt_c(:,5),Xt_c(:,6)];
V_c=[Xt_c(:,10),Xt_c(:,11),Xt_c(:,12)];

[a_t,em_t,i_t,OM_t,om_t,theta0_t]=RV2OrbPar_My(X_0t(1:3),X_0t(4:6),Mu)
[a_c,em_c,i_c,OM_c,om_c,theta0_c]=RV2OrbPar_My(X_0c(1:3),X_0c(4:6),Mu)

[x_E, y_E, z_E] = ellipsoid(0,0,0,r_Earth,r_Earth,r_Earth);


figure()
plot3(R_t(:,1),R_t(:,2),R_t(:,3))
hold on
plot3(R_c(:,1),R_c(:,2),R_c(:,3),'r')
hold on
grid on

figure()
plot3(R_t(:,1),R_t(:,2),R_t(:,3))
hold on
plot3(R_c(:,1),R_c(:,2),R_c(:,3),'r')
hold on
surf(x_E, y_E, z_E,'facecolor','g')
grid on

%chaser position rf:{T,x,y,z}
Rrel_I=R_c-R_t;
[N,~]=size(Rrel_I);
Vrel_I=V_c-V_t;
Rrel_abs=zeros(N,1);
Vrel_abs=zeros(N,1);


for i=1:N
    Rrel_abs(i)= norm(Rrel_I(i,:));
    Vrel_abs(i)=norm(Vrel_I(i,:));
end

Rrel_I_max=max(Rrel_abs);



figure()
plot3((Rrel_I(:,1)),Rrel_I(:,2),Rrel_I(:,3))
hold on 
grid on

%target position rf:{T,r,t,z} (LVLH)
Rrel_R=zeros(N,3);
for i=1:N
[a_t,~,i_t,OM_t,om_t,theta_t]=RV2OrbPar_My(R_t,V_t,Mu);
[Rrel_R(i,:),~]=In2LVLH_rel(Rrel_I(i,:),Vrel_I(i,:),a_t,i_t,OM_t,om_t,theta_t,Mu);
end


Rrel_R_max=max(Rrel_abs);

figure()
plot3((Rrel_R(:,1)),Rrel_R(:,2),Rrel_R(:,3))
hold on 
grid on




figure()
plot3(R_t(:,1),R_t(:,2),R_t(:,3),'--b',R_c(:,1),R_c(:,2),R_c(:,3),'--r')
hold on
grid on

figure()
plot3(Rrel_R(:,1),Rrel_R(:,2),Rrel_R(:,3),'b',Rrel_R(1,1),Rrel_R(1,2),Rrel_R(1,3),'ob',Rrel_R(end,1),Rrel_R(end,2),Rrel_R(end,3),'or')
title('relative position RF: centered on t (LVLH) || exact integration')
hold on
grid on

figure()
plot3(Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'k',Rrel_I(1,1),Rrel_I(1,2),Rrel_I(1,3),'ob',Rrel_I(end,1),Rrel_I(end,2),Rrel_I(end,3),'or')
title('relative position RF: centered on t (I) || exact integration')
hold on
grid on

%% comparison with Cholessy-wilthshare model
%--> initial condition must be given according to the assumption of
% In-plane motion
% Circular orbits

t_end_CW=1/2*T_orb;


%exact integration
R_0c=(R_0t-[rand(1,1)*100;0;0]); %-> pert. only in plane: circular orbit
V_0c=[0;sqrt(Mu/norm(R_0c));0];




[a_c,ec_c,i_c,OM_c,om_c,teta0_c]=RV2OrbPar_My(R_0c,V_0c,Mu)
options = odeset('RelTol',1e-13,'AbsTol', 1e-13);
[t_int,Xt_c]= ode113('r2BP_t_c',t_end_CW,X_0t_c,options);

R_t=[Xt_c(:,1),Xt_c(:,2),Xt_c(:,3)];
R_c=[Xt_c(:,7),Xt_c(:,8),Xt_c(:,9)];
V_t=[Xt_c(:,4),Xt_c(:,5),Xt_c(:,6)];
V_c=[Xt_c(:,10),Xt_c(:,11),Xt_c(:,12)];


Rrel_I=R_c-R_t;
[N,~]=size(Rrel_I);

[N,~]=size(Rrel_I);
Vrel_I=V_c-V_t;


Rrel_R=zeros(N,3);
for i=1:N
  [Rrel_R(i,:),~]=In2LVLH(R_t(i,:),V_t(i,:),Mu);
end



 %Cholessy-Wilthshare
Rho_0=R_0c-R_0t;
Ni_0=V_0c-V_0t;
x_0=[Rho_0,Ni_0];

omega=sqrt(Mu/(norm(R_0t)^3));

setGlobalW(omega);
 
[t_cw,x]= ode113('Chol_Wilt',t_end_CW,x_0,options);

Rho=[x(:,1),x(:,2),x(:,3)];


figure()
plot3(R_t(:,1),R_t(:,2),R_t(:,3),'--b',R_c(:,1),R_c(:,2),R_c(:,3),'--r')
hold on
grid on

figure()
plot3(Rho(:,1),Rho(:,2),Rho(:,3),'-k')
title('CW solution')
hold on 
grid on

figure()
plot3((Rrel_I(:,1)),Rrel_I(:,2),Rrel_I(:,3))
title('comparison C-W, exact integration in Inertial RF')
hold on 
grid on

figure()
plot3((Rrel_I(:,1)),Rrel_I(:,2),Rrel_I(:,3),Rho(:,1),Rho(:,2),Rho(:,3),'-k')
title('comparison C-W, exact integration in Inertial RF')
hold on 
grid on

figure()
plot3((Rrel_R(:,1)),Rrel_R(:,2),Rrel_R(:,3),Rho(:,1),Rho(:,2),Rho(:,3),'-k')
title('comparison C-W, exact integration in in rotating RF (LVLH)')
hold on 
grid on


% BRutal TEst via KEPlerian PARameters & ANAlitical solution [JUST TO UNDERSTAND CW-RF]


[a_t,em_t,i_t,OM_t,om_t,theta0_t]=RV2OrbPar_My(R_0t,V_0t,Mu)
[a_c,em_c,i_c,OM_c,om_c,theta0_c]=RV2OrbPar_My(R_0c,V_0c,Mu)



R_t_AN=Plotorb(a_t,em_t,i_t,OM_t,om_t,theta0_t,t_cw,Mu);
R_c_AN=Plotorb(a_c,em_c,i,OM_c,om_c,theta0_c,t_cw,Mu);

R_t_AN=Plotorb(a_t,0,0,0,0,0,t_cw,Mu);
R_c_AN=Plotorb(a_c,0,0,0,0,0,t_cw,Mu);


figure()
plot3((R_t_AN(:,1)),R_t_AN(:,2),R_t_AN(:,3),'-b',R_c_AN(:,1),R_c_AN(:,2),R_c_AN(:,3),'-r')
title('KEPlerian ORBits CHAser & TARget')
hold on
grid on

LOCC=R_t_AN+Rho;
figure()
plot3((R_t_AN(:,1)),R_t_AN(:,2),R_t_AN(:,3),'-b',R_c_AN(:,1),R_c_AN(:,2),R_c_AN(:,3),'-r',LOCC(:,1),LOCC(:,2),LOCC(:,3),'k')
title('Condjungnet CW')
hold on
grid on












%% active feedback control:linear ch-wi -> A  & [ideal {ax,ay,az}] & steady-state tracking

w=2*pi/T_orb;

%% active feedback control:linear ch-wi -> A  & [ideal {ax,ay,az}] & reference tracking
