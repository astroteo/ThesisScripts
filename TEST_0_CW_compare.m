%% comparison with Cholessy-wilthshare model
%--> initial condition must be given according to the assumption of
% In-plane motion
clear all
close all

% Circular orbits (in-plane motion)
Mu=(3.9856e+14);
r_Earth=6380*10^3;
h_orb=200*10^3;


R_0t=[(r_Earth+h_orb);0;0];
V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];
[a_t,e_t,i_t,OM_t,om_t,theta_t]=RV2OrbPar_My(R_0t,V_0t,Mu)

R_0c=(R_0t-[100 ;0;0]); %-> pert. only in plane: circular orbit
V_0c=[0;sqrt(Mu/norm(R_0c));0];
[a_c,e_c,i_c,OM_c,om_c,teta0_c]=RV2OrbPar_My(R_0c,V_0c,Mu)

T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)
AN_or_INT=1 %1-> integration 0-> analytical orbits

%exact integration
X_0t=[R_0t;V_0t];
X_0c=[R_0c;V_0c];

X_0t_c=[X_0t;X_0c];
t_end_CW=1/15*T_orb

if AN_or_INT==1;
options = odeset('RelTol',1e-13,'AbsTol', 1e-13);
[t_int,Xt_c]= ode113('r2BP_t_c',t_end_CW,X_0t_c,options);

R_t=[Xt_c(:,1),Xt_c(:,2),Xt_c(:,3)];
R_c=[Xt_c(:,7),Xt_c(:,8),Xt_c(:,9)];
V_t=[Xt_c(:,4),Xt_c(:,5),Xt_c(:,6)];
V_c=[Xt_c(:,10),Xt_c(:,11),Xt_c(:,12)];


Rrel_I=R_c-R_t;
Vrel_I=V_c-V_t;
[N,~]=size(Rrel_I);

else
    
 N=10000;
time_span=linspace(0,t_end_CW,N);
t_int=time_span;
RT=norm(R_0t);
RC=norm(R_0c);
R_t=zeros(N,3);
V_t=zeros(N,3);
R_c=zeros(N,3);
V_c=zeros(N,3);
Rrel_I=zeros(N,3);
Vrel_I=zeros(N,3);
for k=1:N
theta=2*pi/T_orb*time_span(k);
R_t(k,:)=RT*[cos(theta),sin(theta),0];
V_t(k,:)=2*pi/T_orb*RT*[-sin(theta),cos(theta),0];
R_c(k,:)=RC*[cos(theta),sin(theta),0];
V_c(k,:)=2*pi/T_orb*RC*[-sin(theta),cos(theta),0];
Rrel_I(k,:)=R_c(k,:)-R_t(k,:);
Vrel_I(k,:)=V_c(k,:)-V_t(k,:);
end

end



    
Rrel_R=zeros(N,3);
for i=1:N
    [a_t,e_t,i_t,OM_t,om_t,theta_t]=RV2OrbPar_My(R_t(i,:),V_t(i,:),Mu);
    [Rrel_R(i,:),~]=In2LVLH_Par(Rrel_I(i,:),Vrel_I(i,:),a_t,e_t,i_t,OM_t,om_t,theta_t,Mu);
    
    %[Rrel_R(i,:),~]=In2LVLH(R_t(i,:),R_t(i,:),Rrel_I(i,:),Vrel_I(i,:),Mu);
end

  
    
    



%Cholessy-Wilthshare
Rho_0=R_0c-R_0t;
Ni_0=V_0c-V_0t;
Rho_0=[Rho_0(2);Rho_0(1);Rho_0(3)];
Ni_0=[Ni_0(2);Ni_0(1);Ni_0(3)];

x_0=[Rho_0,Ni_0];

omega=sqrt(Mu/(norm(R_0t)^3));

setGlobalW(omega);
 
options = odeset('RelTol',1e-13,'AbsTol', 1e-13);
[t_cw,x]= ode113('Chol_Wilt',t_end_CW,x_0,options);

Rho=[x(:,2),x(:,1),x(:,3)];


%comparison with analytical solution of CW-equations
[M,~]=size(t_cw);
x_AN=zeros(M,3);

x_CW_norm=zeros(M,1);
for i=1:M
    t=t_cw(i);
    
x_AN(i,1)=(4*Ni_0(1)/omega-6*omega*Rho_0(2))*sin(omega*t)-2/omega*Ni_0(2)*cos(omega*t)+(6*omega*Rho_0(2)-3*Ni_0(1))*t+Rho_0(1)+2*Ni_0(2)/omega;
x_AN(i,2)=(2*Ni_0(1)/omega-3*Rho_0(2))*cos(omega*t)+Ni_0(2)/omega*sin(omega*t)+4*Rho_0(2)-2*Ni_0(1)/omega;
x_AN(i,3)=Rho_0(3)*cos(omega*t)-Ni_0(3)*omega*sin(omega*t);

x_CW_norm(i)=norm([x_AN(1),x_AN(2),x_AN(3)]);
    end

x_AN=[x_AN(:,2),x_AN(:,1),x_AN(:,3)];

%brutal analytic rotation: In->LVLH
Rrel_BRUT=zeros(N,3);
Rrel_BRUT_norm=zeros(N,1);
for j=1:N
    t=t_int(j);
    theta=-omega*t;
ROT_Theta=[cos(theta), sin(theta),  0;
           -sin(theta), cos(theta),  0;
             0       ,     0     ,  1];
 Rrel_BRUT(j,:)=ROT_Theta'* Rrel_I(j,:)';
% theta_r=acos(dot(Rrel_I(j,:),[1,0,0])/norm(Rrel_I(j,:)))
% Rrel_mod=norm(Rrel_I(j,:))
% Rrel_BRUT(j,1)=cos(theta_r-theta)*Rrel_mod;
% Rrel_BRUT(j,2)=sin(theta_r-theta)*Rrel_mod;
%  Rrel_BRUT(j,3)=Rrel_I(j,3)
Rrel_BRUT_norm(j)=norm(Rrel_I(j,:));
end

%brutal analytic rotation: LVLH-->In
Rrel_BRUT_I=zeros(M,3);
Rrel_BRUT_I_norm=zeros(M,1);
for j=1:M
    t=t_int(j);
    theta=-omega*t;
ROT_Theta=[cos(theta), sin(theta),  0;
           -sin(theta), cos(theta),  0;
             0       ,     0     ,  1];
 Rrel_BRUT_I(j,:)=ROT_Theta* Rho(j,:)';
% theta_r=acos(dot(Rrel_I(j,:),[1,0,0])/norm(Rrel_I(j,:)))
% Rrel_mod=norm(Rrel_I(j,:))
% Rrel_BRUT(j,1)=cos(theta_r-theta)*Rrel_mod;
% Rrel_BRUT(j,2)=sin(theta_r-theta)*Rrel_mod;
%  Rrel_BRUT(j,3)=Rrel_I(j,3)
Rrel_BRUT_I_norm(j)=norm(Rrel_BRUT_I(j,:));
end

figure()
plot3(Rho(:,1),Rho(:,2),Rho(:,3),'b',x_AN(:,1),x_AN(:,2),x_AN(:,3),'--k', Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or',x_AN(1,1),x_AN(1,2),x_AN(1,3),'og',x_AN(end,1),x_AN(end,2),x_AN(end,3),'or')
title('analytical vs integrated CW-equations')
hold on
grid on 

figure()
plot3(R_t(:,1),R_t(:,2),R_t(:,3),'--b',R_c(:,1),R_c(:,2),R_c(:,3),'--r')
hold on
grid on

figure()
plot3(Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3))
title('exact integration in Inertial RF')
hold on 
grid on

figure()
plot3(Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),Rho(:,1),Rho(:,2),Rho(:,3),'-k')
title('comparison ! INT C-W ! , exact integration in Inertial RF (IN)')
hold on 
grid on

figure()
plot3(Rrel_R(:,1),Rrel_R(:,2),Rrel_R(:,3),'b',Rho(:,1),Rho(:,2),Rho(:,3),'-k')
title('comparison !INT C-W !, exact integration in in rotating RF (LVLH)')
hold on 
grid on


figure()
plot3(Rrel_BRUT(:,1),Rrel_BRUT(:,2),Rrel_BRUT(:,3),'b',Rho(:,1),Rho(:,2),Rho(:,3),'--k')
title('comparison  !INT C-W!, exact integration in in rotating RF(LVLH) !BRUTAL ROT!')
hold on 
grid on

figure()
plot3(Rrel_BRUT(:,1),Rrel_BRUT(:,2),Rrel_BRUT(:,3),'b',x_AN(:,1),x_AN(:,2),x_AN(:,3),'--k')
title('comparison !ANALYTICAL C-W!, exact integration in in rotating RF(LVLH) !BRUTAL ROT!')
hold on 
grid on
figure()
plot3(Rrel_BRUT_I(:,1),Rrel_BRUT_I(:,2),Rrel_BRUT_I(:,3),'--k',Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'b')
title('comparison !ANALYTICAL C-W!, exact integration in in inerital RF(IN) !BRUTAL ROT!')
hold on 
grid on


figure()
plot3(Rrel_BRUT(:,1),Rrel_BRUT(:,2),Rrel_BRUT(:,3),'b')
title('Rrel rotating RF(LVLH) !BRUTAL ROT!')
hold on 
grid on

figure()
plot3(Rrel_BRUT(:,1),Rrel_BRUT(:,2),Rrel_BRUT(:,3),'b',Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'--r')
title('UNDESTAND ROTATION: Rrel_I (r) & Rrel_R(b) RF(LVLH) !BRUTAL ROT!')
hold on 
grid on
%% order of magnitude test

%c.i. 
Delta_r=[1 10 100]
Delta_v=omega*Delta_r
R_0t=[(r_Earth+h_orb);0;0];
V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];


R_0c=[(r_Earth+h_orb+1);0;0];
V_0c=[0;sqrt(Mu/(r_Earth+h_orb+1));0];

DV=norm(V_0c-V_0t)

R_0c=[(r_Earth+h_orb-100);0;0];
V_0c=[0;sqrt(Mu/(r_Earth+h_orb-100));0];

DV=norm(V_0c-V_0t)
