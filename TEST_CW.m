%% cholessy-wilthshare test
close all
clear all
%orbit
Mu=(3.9856e+14);
r_Earth=6380*10^3;
h_orb=200*10^3;
T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)
%target
R_0t=[(r_Earth+h_orb);0;0];
V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];

t_end_CW=1/4*T_orb

omega=sqrt(Mu/(norm(R_0t)^3));
omega=omega; %?????? + OR - ?????%
d_0=[1 10 50 100 200 ]
[~,T_d]=size(d_0);

%% model's hypotesis.--> ONLY LINEARITY LOSS

for k=1:T_d
% d_TC=d_0(k);
d_TC=d_0(k);
%chaser
R_0c=(R_0t+[d_TC;0;0]); %-> pert. only in plane: circular orbit
V_0c=[0;sqrt(Mu/norm(R_0c));0];


%exact integration
X_0t=[R_0t;V_0t];
X_0c=[R_0c;V_0c];

X_0t_c=[X_0t;X_0c];

options = odeset('RelTol',1e-13,'AbsTol', 1e-20);
[t_int,Xt_c]= ode113('r2BP_t_c',t_end_CW,X_0t_c,options);


R_t=[Xt_c(:,1),Xt_c(:,2),Xt_c(:,3)];
R_c=[Xt_c(:,7),Xt_c(:,8),Xt_c(:,9)];
V_t=[Xt_c(:,4),Xt_c(:,5),Xt_c(:,6)];
V_c=[Xt_c(:,10),Xt_c(:,11),Xt_c(:,12)];

Rrel_I=R_c-R_t;
Vrel_I=V_c-V_t;

[M,~]=size(t_int);
Rrel_Inorm=zeros(M,1);
Rrel_LVLH=zeros(M,3);
Vrel_LVLH_IN=zeros(M,3);
Vrel_LVLH_ROT=zeros(M,3);
Vrel_LVLH_R=zeros(M,3);

for j=1:M
Rrel_Inorm(j)=norm(Rrel_I(j,:));
t=t_int(j);
theta=omega*t;
ROT=[cos(theta), sin(theta),0; 
    -sin(theta),cos(theta), 0; 
      0        ,     0    , 1];
  
Rrel_LVLH(j,:)=(ROT'*Rrel_I(j,:)')';

Vrel_LVLH_IN(j,:)=Vrel_I(j,:)-cross([0;0;omega],Rrel_I(j,:));
Vrel_LVLH_R(j,:)=(ROT*Vrel_LVLH_IN(j,:)')';
Vrel_LVLH_ROT(j,:)=(ROT*Vrel_I(j,:)')'-(cross([0;0;omega],Rrel_LVLH(j,:)'))';



end

max_dist_T_C_INT=max(Rrel_Inorm);


%CW-integration.
Rho_0=R_0c-R_0t;
Rho_0=[Rho_0(2);Rho_0(1);Rho_0(3)];

Ni_0=((V_0c-V_0t)-cross([0;0;omega],(R_0c-R_0t)));
Ni_0_corr=[Ni_0(2);Ni_0(1);Ni_0(3)];
x_0=[Rho_0;Ni_0_corr];




setGlobalW(omega);
 
options = odeset('RelTol',1e-13,'AbsTol', 1e-20);
[t_cw,x]= ode113('Chol_Wilt',t_end_CW,x_0,options);

Rho=[x(:,2),x(:,1),x(:,3)];
 Ni=[x(:,5),x(:,4),x(:,6)];

[N,~]=size(t_cw);
Rrel_CWnorm=zeros(M,1);


for i=1:N
Rrel_CWnorm(i)=norm(Rho(i,:));

end

max_dist_T_C_CW=max(Rrel_CWnorm);

model_error(k)= abs(max_dist_T_C_CW-max_dist_T_C_INT);

end



model_error

figure()
plot3(R_t(:,1),R_t(:,2),R_t(:,3),'--b',R_c(:,1),R_c(:,2),R_c(:,3),'b',R_c(1,1),R_c(1,2),R_c(1,3),'og',R_c(end,1),R_c(end,2),R_c(end,3),'or',R_t(1,1),R_t(1,2),R_t(1,3),'og',R_t(end,1),R_t(end,2),R_t(end,3),'or')
hold on
grid on

figure()
plot3(Rho(:,1),Rho(:,2),Rho(:,3),'--k',Rrel_LVLH(:,1),Rrel_LVLH(:,2),Rrel_LVLH(:,3),'b', Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or',Rrel_LVLH(1,1),Rrel_LVLH(1,2),Rrel_LVLH(1,3),'og',Rrel_LVLH(end,1),Rrel_LVLH(end,2),Rrel_LVLH(end,3),'or')
title('!CW-INT! vs !Rrel_{LVLH}-INT!')
hold on
grid on 

figure()
plot3(Rho(:,1),Rho(:,2),Rho(:,3),'--k',Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'b', Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or',Rrel_I(1,1),Rrel_I(1,2),Rrel_I(1,3),'og',Rrel_I(end,1),Rrel_I(end,2),Rrel_I(end,3),'or')
title('!CW-INT! vs !Rrel_{IN}-INT!')
hold on
grid on

figure()
plot3(Ni(:,1),Ni(:,2),Ni(:,3),'--k',Vrel_LVLH_IN(:,1),Vrel_LVLH_IN(:,2),Vrel_LVLH_IN(:,3),'r')
title('!CW-INT! vs !Vrel_I-INT! ')
hold on
grid on 

figure()
plot3(Ni(:,1),Ni(:,2),Ni(:,3),'--k',Vrel_LVLH_ROT(:,1),Vrel_LVLH_ROT(:,2),Vrel_LVLH_ROT(:,3),'r')
title('!CW-INT! vs !Vrel_{LVLH-ROT}-INT! ')
hold on
grid on 

figure()
plot3(Ni(:,1),Ni(:,2),Ni(:,3),'--k',Vrel_LVLH_R(:,1),Vrel_LVLH_R(:,2),Vrel_LVLH_R(:,3),'r')
title('!CW-INT! vs !Vrel_{LVLH-R}-INT! ')
hold on
grid on 
omega_test_ortdox=omega

% %% random initial condition --> APPROXIMATION
% for k=1:T_d
% d_TC=d_0(k);
% %target
% R_0t=[(r_Earth+h_orb);0;0];
% V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];
% 
% %chaser
% R_0c=(R_0t+d_TC*rand(3,1)); %-> pert. only in plane: circular orbit
% V_0c=sqrt(Mu/norm(R_0c))*[0;1;0]-(sqrt(Mu/norm(R_0t))-sqrt(Mu/norm(R_0c)))*rand(3,1);
% 
% % R_0c=(R_0t-rand(3,1)*1e0);%->  pert. random
% % V_0c=(V_0t-rand(3,1)*1e2);
% 
% %exact integration
% X_0t=[R_0t;V_0t];
% X_0c=[R_0c;V_0c];
% 
% X_0t_c=[X_0t;X_0c];
% 
% options = odeset('RelTol',1e-13,'AbsTol', 1e-13);
% [t_int,Xt_c]= ode113('r2BP_t_c',t_end_CW,X_0t_c,options);
% 
% 
% R_t=[Xt_c(:,1),Xt_c(:,2),Xt_c(:,3)];
% R_c=[Xt_c(:,7),Xt_c(:,8),Xt_c(:,9)];
% V_t=[Xt_c(:,4),Xt_c(:,5),Xt_c(:,6)];
% V_c=[Xt_c(:,10),Xt_c(:,11),Xt_c(:,12)];
% 
% Rrel_I=R_c-R_t;
% Vrel_I=V_c-V_t;
% 
% [M,~]=size(t_int);
% Rrel_Inorm=zeros(M,1);
% 
% 
% for j=1:M
% Rrel_Inorm(j)=norm(Rrel_I(j,:));
% 
% end
% 
% max_dist_T_C_INT=max(Rrel_Inorm);
% 
% 
% %CW-integration.
% Rho_0=R_0c-R_0t;
% Rho_0=[Rho_0(2);Rho_0(1);Rho_0(3)];
% Ni_0=V_0c-V_0t;
% Ni_0=[Ni_0(2);Ni_0(1);Ni_0(3)];
% x_0=[Rho_0,Ni_0];
% 
% 
% options = odeset('RelTol',1e-13,'AbsTol', 1e-15);
% [t_cw,x]= ode113('Chol_Wilt',t_end_CW,x_0,options);
% 
% Rho=[x(:,2),x(:,1),x(:,3)];
% Ni=[x(:,5),x(:,4),x(:,6)];
% 
% [N,~]=size(t_cw);
% Rrel_CWnorm=zeros(M,1);
% 
% 
% for i=1:N
% Rrel_CWnorm(i)=norm(Rho(i,:));
% 
% end
% 
% max_dist_T_C_RND=max(Rrel_CWnorm);
% 
% model_error_rnd(k)= abs(max_dist_T_C_RND-max_dist_T_C_CW);
% 
% end
% 
% 
% figure()
% plot3(Rho(:,1),Rho(:,2),Rho(:,3),'--k',Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'b', Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or',Rrel_I(1,1),Rrel_I(1,2),Rrel_I(1,3),'og',Rrel_I(end,1),Rrel_I(end,2),Rrel_I(end,3),'or')
% title('!CW-INT! vs !Rrel_{IN}-INT!')
% hold on
% grid on
% 
