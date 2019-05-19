clear all
close all

%CW-integrated
Mu=(3.9856e+14);
r_Earth=6380*10^3;
h_orb=200*10^3;
T_orb=2*pi*sqrt((r_Earth+h_orb)^3/Mu)
%target
R_T=r_Earth+h_orb;
R_0t=[R_T;0;0];
V_0t=[0;sqrt(Mu/(r_Earth+h_orb));0];


omega=sqrt(Mu/R_T^3);
omega=omega; %?????? + OR - ?????%
t_end_CW=1/3*T_orb

V_T=omega*R_T;
%chaser
d_TC=60000;

R_C=R_T-d_TC;
omega_c=sqrt(Mu/(R_C^3));
V_C=omega_c*R_C;
R_0c=[R_C;0;0]; %-> pert. only in plane: circular orbit
V_0c=[0;sqrt(Mu/norm(R_0c));0];

%CW-integrated.
Rho_0=R_0c-R_0t;
Rho_0=[Rho_0(2);Rho_0(1);Rho_0(3)]

Ni_0=((V_0c-V_0t)-cross([0;0;omega],[R_C-R_T;0;0]));
Ni_0=[Ni_0(2);Ni_0(1);Ni_0(3)];

x_0=[Rho_0;Ni_0];

setGlobalW(omega);

 options = odeset('RelTol',1e-13,'AbsTol', 1e-20);
[t_cw,x]= ode113('Chol_Wilt',t_end_CW,x_0,options);

Rho=[x(:,1),x(:,2),x(:,3)];
 Ni=[x(:,4),x(:,5),x(:,6)];
 
 
%C_W analytical

[N,~]=size(t_cw);
 x_AN=zeros(N,3);
% B=Ni_0(2);
% A=3*omega*Rho_0(2)-2*Ni_0(1);
% C=4*Rho_0(2)-2*Ni_0(1)/omega;
% D=6*Rho_0(2)*omega-3*Ni_0(1)/omega;
% E=Rho_0(1)+2*Ni_0(2)/omega;
for i=1:N
    t=t_cw(i);
% x_AN(i,1)=-2*A/omega*sin(omega*t)-2*B/omega+D*t+E;
% x_AN(i,2)=-A/omega*cos(omega*t)+B/omega*sin(omega*t)+C;
% x_AN(i,3)=Rho_0(3)*cos(omega*t)-Ni_0(3)*omega*sin(omega*t);

x_AN(i,1)=(4*Ni_0(1)/omega-6*Rho_0(2))*sin(omega*t)-2/omega*Ni_0(2)*cos(omega*t)+(6*omega*Rho_0(2)-3*Ni_0(1))*t+Rho_0(1)+2*Ni_0(2)/omega;
x_AN(i,2)=(2*Ni_0(1)/omega-3*Rho_0(2))*cos(omega*t)+Ni_0(2)/omega*sin(omega*t)+4*Rho_0(2)-2*Ni_0(1)/omega;
x_AN(i,3)=Rho_0(3)*cos(omega*t)-Ni_0(3)*omega*sin(omega*t);

end
x_AN=[x_AN(:,1),x_AN(:,2),x_AN(:,3)];

figure()
plot3(Rho(:,1),Rho(:,2),Rho(:,3),'b',x_AN(:,1),x_AN(:,2),x_AN(:,3),'--k',Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or',x_AN(1,1),x_AN(1,2),x_AN(1,3),'og',x_AN(end,1),x_AN(end,2),x_AN(end,3),'or')
title('!CW-INT! vs !CW-AN!')
hold on
grid on 

%analytical orbits
R_t=zeros(N,3);
R_c=zeros(N,3);
V_t=zeros(N,3);
V_c=zeros(N,3);

for i=1:N
    t=t_cw(i);
    theta=omega*t;
    theta_c=omega_c*t;
    R_t(i,:)=[R_T*cos(theta),R_T*sin(theta),0];
    R_c(i,:)=[R_C*cos(theta_c),R_C*sin(theta_c),0];
    V_t(i,:)=[V_T*(-sin(theta)),V_T*cos(theta),0];
    V_c(i,:)=[V_C*(-sin(theta_c)),V_C*cos(theta_c),0];
    
end
figure()
plot3(R_t(:,1),R_t(:,2),R_t(:,3),'b',R_c(:,1),R_c(:,2),R_c(:,3),'r',R_t(1,1),R_t(1,2),R_t(1,3),'og',R_t(end,1),R_t(end,2),R_t(end,3),'or',R_c(1,1),R_c(1,2),R_c(1,3),'og',R_c(end,1),R_c(end,2),R_c(end,3),'or')
title('!R_t_{in}! vs !R_c_{in}!')
hold on
grid on 

Rrel_I=R_c-R_t;
Vrel_I=V_c-V_t;

Rrel_LVLH=zeros(N,3);
Rrel_I_rec=zeros(N,3);
for i=1:N
    t=t_cw(i);
    theta=omega*t;
    ROT=[sin(theta), -cos(theta),0;
           cos(theta),sin(theta),0;
             0       ,    0     ,1];
        Rrel_LVLH(i,:)= (ROT*Rrel_I(i,:)')';
        Rrel_I_rec(i,:)=(ROT'*Rho(i,:)')';
end

figure()
plot3(Rho(:,1),Rho(:,2),Rho(:,3),'--k',Rrel_LVLH(:,1),Rrel_LVLH(:,2),Rrel_LVLH(:,3),'r',Rrel_LVLH(1,1),Rrel_LVLH(1,2),Rrel_LVLH(1,3),'og',Rrel_LVLH(end,1),Rrel_LVLH(end,2),Rrel_LVLH(end,3),'or',Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or')
title('!CW-INT! vs !CIRC-diff! in LVLH ref')
hold on
grid on 

figure()
plot3(Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'b',Rrel_I_rec(:,1),Rrel_I_rec(:,2),Rrel_I_rec(:,3),'--k',Rrel_I_rec(1,1),Rrel_I_rec(1,2),Rrel_I_rec(1,3),'og',Rrel_I_rec(end,1),Rrel_I_rec(end,2),Rrel_I_rec(end,3),'or',Rrel_I(1,1),Rrel_I(1,2),Rrel_I(1,3),'og',Rrel_I(end,1),Rrel_I(end,2),Rrel_I(end,3),'or')
title('!CW-INT! vs !CIRC-diff!in IN ref')
hold on
grid on 

figure()
plot3(Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'b',Rho(:,1),Rho(:,2),Rho(:,3),'--k',Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or',Rrel_I(1,1),Rrel_I(1,2),Rrel_I(1,3),'og',Rrel_I(end,1),Rrel_I(end,2),Rrel_I(end,3),'or')
title('!CW-INT! vs !CIRC-diff!in IN ref')
hold on
grid on 
