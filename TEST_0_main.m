%[ TEST 0] trajectory tracking 
clear all
close all

%% 1)verify that relative motion can be ignored from orbital POV

Mu=3.9856e+14;
r_Earth=6380*10^3;
h_orb=400*10^3;

R_0t=[(r_Earth+h_orb);0;0];
V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];

[~,~,i_0,OM_0,om_0,theta_0]=RV2OrbPar_My(R_0t,V_0t,Mu);
[R_0_LVLH,V_0_LVLH]=In2LVLH(R_0t,V_0t,i_0,OM_0,om_0,theta_0);

X_0t=[R_0t;V_0t];
T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)

% R_0c=(R_0t+[rand(1,1);0;0]*1e1);
% V_0c=(V_0t+[0;rand(1,1);0]*1e-1);

R_0c=(R_0t+rand(3,1)*1e1);
V_0c=(V_0t+rand(3,1)*1e-1);



Rrel_0_I=R_0c-R_0t

X_0c=[R_0c;V_0c]; %-> chaser [SPHERE] is assumed to move up to 1 dm/s in difference from the target at a maximum distnce of 10m
                                                        
N=1e4;
t_end=1*T_orb
tspan=linspace(0,t_end,N);


% [t,X_t]= ode113('r2BP',tspan,X_0t);
% [t,X_c]= ode113('r2BP',tspan,X_0c);
% options = odeset('RelTol',1e-9,'AbsTol', 1e-9);
% [t,X_t]= ode113('r2BP',10*T_orb,X_0t,options);
% [t,X_c]= ode113('r2BP',10*T_orb,X_0c,options);
% 
% R_t=[X_t(:,1),X_t(:,2),X_t(:,3)];
% R_c=[X_c(:,1),X_c(:,2),X_c(:,3)];
% V_c=[X_c(:,4),X_c(:,5),X_c(:,6)];

% [a_t,em_t,i_t,OM_t,om_t,theta0_t]=RV2OrbPar_My(X_0t(1:3),X_0t(4:6),Mu)
% [a_c,em_c,i_c,OM_c,om_c,theta0_c]=RV2OrbPar_My(X_0c(1:3),X_0c(4:6),Mu)
% 
% R_t=zeros(N,3);
% V_t=zeros(N,3);
% 
% R_c=zeros(N,3);
% V_c=zeros(N,3);
% 
% for i=1:N
%     t=tspan(i);
%     
%     [theta_t]=time2theta(a_t,em_t,theta0_t,t,Mu);
%     [R_t(i,:),V_t(i,:)]=ParOrb2RV(a_t,em_t,i_t,OM_t,om_t,theta_t,Mu);
%     
%     [theta_c]=time2theta(a_c,em_c,theta0_c,t,Mu);
%     [R_c(i,:),V_c(i,:)]=ParOrb2RV(a_c,em_c,i_c,OM_c,om_c,theta_c,Mu);
%     
%     
% end
% 
% 
% [x_E, y_E, z_E] = ellipsoid(0,0,0,r_Earth,r_Earth,r_Earth);
% 
% 
% figure()
% plot3(R_t(:,1),R_t(:,2),R_t(:,3))
% hold on
% plot3(R_c(:,1),R_c(:,2),R_c(:,3))
% hold on
% grid on
% 
% figure()
% plot3(R_t(:,1),R_t(:,2),R_t(:,3))
% hold on
% plot3(R_c(:,1),R_c(:,2),R_c(:,3))
% hold on
% surf(x_E, y_E, z_E,'facecolor','g')
% grid on
% 
% %chaser position rf:{T,x,y,z}
% Rrel_I=R_c-R_t;
% Vrel_I=V_c;
% Rrel_abs=zeros(N,1);
% 
% for i=1:N
%     Rrel_abs(i)= norm(Rrel_I(i,:));
% end
% 
% Rrel_I_max=max(Rrel_abs);
% 
% figure()
% plot3((Rrel_I(:,1)),Rrel_I(:,2),Rrel_I(:,3))
% hold on 
% grid on
% 
% %target position rf:{T,r,t,z} (LVLH)
% Rrel_R=zeros(N,3);
% for i=1:N
%     [~,~,i_t,OM_t,om_t,theta_t]=RV2OrbPar_My(R_t(i,:),V_t(i,:),Mu);
%     [Rrel_R(i,:),~]=In2LVLH(Rrel_I(i,:),Vrel_I(i,:),i_t,OM_t,om_t,theta_t);
% end
% 
% 
% Rrel_R_max=max(Rrel_abs);
% 
% figure()
% plot3((Rrel_R(:,1)),Rrel_R(:,2),Rrel_R(:,3))
% hold on 
% grid on



%simulink output

%out=sim('TEST_0_sim');
%out=sim('TEST_0_sim_NO_3D');
out=('ubutest')



[~,~,M]=size(R_t_sim.signals.values);

R_t_int=zeros(M,3);
R_c_int=zeros(M,3);
Rrel_I_int=zeros(M,3);
Rrel_R_int=zeros(M,3);

for i=1:M
    R_t_int(i,:)=R_t_sim.signals.values(:,1,i);
    R_c_int(i,:)=R_c_sim.signals.values(:,1,i);
    Rrel_I_int(i,:)=Rrel_I_sim.signals.values(:,1,i);
    Rrel_R_int(i,:)=Rrel_R_sim.signals.values(i,:);
    
end

figure()
plot3(R_t_int(:,1),R_t_int(:,2),R_t_int(:,3),'--b',R_c_int(:,1),R_c_int(:,2),R_c_int(:,3),'--r')
hold on
grid on

figure()
plot3(Rrel_R_int(:,1),Rrel_R_int(:,2),Rrel_R_int(:,3),'b',Rrel_R_int(1,1),Rrel_R_int(1,2),Rrel_R_int(1,3),'ob',Rrel_R_int(end,1),Rrel_R_int(end,2),Rrel_R_int(end,3),'or')
hold on
grid on

figure()
plot3(Rrel_I_int(:,1),Rrel_I_int(:,2),Rrel_I_int(:,3),'k',Rrel_I_int(1,1),Rrel_I_int(1,2),Rrel_I_int(1,3),'ob',Rrel_I_int(end,1),Rrel_I_int(end,2),Rrel_I_int(end,3),'or')
hold on
grid on


%% active feedback control:linear ch-wi -> A  & [ideal {ax,ay,az}] & steady-state tracking

w=2*pi/T_orb;

%% active feedback control:linear ch-wi -> A  & [ideal {ax,ay,az}] & reference tracking
