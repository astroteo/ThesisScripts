
close all
clear all
%orbit
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

%chaser
R_0c=(R_0t-[d_TC;0;0]); %-> pert. only in plane: circular orbit
R_C=R_T-d_TC
omega_c=sqrt(Mu/(norm(R_0c)^3));
V_0c=[0;sqrt(Mu/norm(R_0c));0];
V_0c=[0;omega*R_C;0];%-> to set zero velocity between chaser and target
R_C=R_T-d_TC;
V_C=norm(V_0c);
%Cholessy-Wiltshare integration
Rho_0=R_0c-R_0t;


w=omega;
setGlobalW(w);

Ni_0=(V_0c-V_0t)-cross([0;0;w],(R_0c-R_0t));
X_0t=[R_0t;V_0t];
X_0c=[R_0c;V_0c];

X_0t_c=[X_0t;X_0c];

options = odeset('RelTol',1e-13,'AbsTol', 1e-40);
[t_int,Xt_c]= ode113('r2BP_t_c',t_end_CW,X_0t_c,options);


R_t=[Xt_c(:,1),Xt_c(:,2),Xt_c(:,3)];
R_c=[Xt_c(:,7),Xt_c(:,8),Xt_c(:,9)];
V_t=[Xt_c(:,4),Xt_c(:,5),Xt_c(:,6)];
V_c=[Xt_c(:,10),Xt_c(:,11),Xt_c(:,12)];

%analytical Cholessy-Wiltshare
[N,~]=size(t_int);
x_AN=zeros(N,3);
x_dot_AN=zeros(N,3);

Rho_0=R_0c-R_0t;


w=omega;
setGlobalW(w);

Ni_0=(V_0c-V_0t)-cross([0;0;w],(R_0c-R_0t));



for i=1:N
    t=t_int(i);
   
    x_AN(i,1)=(-2/w*Ni_0(2)-3*Rho_0(1))*cos(w*t)+Ni_0(1)/w*sin(w*t)+4*Rho_0(1)+2/w*Ni_0(2);
    x_AN(i,2)=-1*((-4/w*Ni_0(2)-6*Rho_0(1))*sin(w*t)-2/w*Ni_0(1)*cos(w*t)+(6*w*Rho_0(1)+3*Ni_0(2))*t-Rho_0(2)+2/w*Ni_0(1));
    x_AN(i,3)=Rho_0(3)*cos(w*t)+Ni_0(3)/w*sin(w*t);
    
    x_dot_AN(i,1)=(3*Rho_0(1)*w+2*Ni_0(2))*sin(w*t)+Ni_0(1)*cos(w*t);
    x_dot_AN(i,2)=-((-4*Ni_0(2)-6*w*Rho_0(1))*cos(w*t)+2*Ni_0(1)*sin(w*t)+6*w*Rho_0(1)+3*Ni_0(2) );
    x_dot_AN(i,3)=Ni_0(3)*cos(w*t)-Rho_0(3)*w*sin(w*t);
end

Rho=x_AN;
Ni=x_dot_AN;
t_cw=t_int;

x_0=[Rho_0;Ni_0];
options = odeset('RelTol',1e-13,'AbsTol', 1e-20);
[t_CW,x]= ode113('Chol_Wilt_Hill',t_end_CW,x_0,options);

RHO=[x(:,1),x(:,2),x(:,3)];
 NI=[x(:,4),x(:,5),x(:,6)];


%relative positions & velocities
Rrel_I=R_c-R_t;
Vrel_I=V_c-V_t;

Rrel_LVLH=zeros(N,3);
Vrel_LVLH=zeros(N,3);
Rrel_BASTA=zeros(N,3);
Rrel_I_rec=zeros(N,3);
Vrel_I_rec=zeros(N,3);
%LVLH-->IN
for l=1:N
    t=t_cw(l);
    theta=omega*t;
    ROT=[cos(theta), sin(theta),0;
         -sin(theta), cos(theta),0;
             0       ,    0     ,1];
        
        Rrel_I_rec(l,:)=(ROT'*Rho(l,:)')';
        Vrel_I_rec(l,:)=(ROT'*(Ni(l,:)'+cross([0;0;omega],Rho(l,:)')))';
       
        
end
%IN-->LVLH
for j=1:N
    t=t_cw(j);
    theta=omega*t;
    ROT=[ cos(theta), sin(theta),0;
           -sin(theta), cos(theta),0;
                 0       ,    0       ,1];
             
        Rrel_LVLH(j,:)=(ROT*Rrel_I(j,:)')';
        Vrel_LVLH(j,:)=(ROT*(Vrel_I(j,:)'-cross([0;0;omega],Rrel_I(j,:)')))';
        
end

figure()
plot3(Rho(:,1),Rho(:,2),Rho(:,3),'--r',Rrel_LVLH(:,1),Rrel_LVLH(:,2),Rrel_LVLH(:,3),'--b',Rrel_LVLH(1,1),Rrel_LVLH(1,2),Rrel_LVLH(1,3),'og',Rrel_LVLH(end,1),Rrel_LVLH(end,2),Rrel_LVLH(end,3),'or',Rho(1,1),Rho(1,2),Rho(1,3),'og',Rho(end,1),Rho(end,2),Rho(end,3),'or')
title('!CW-INT! vs !CIRC-diff! in LVLH ref')
hold on
grid on 
figure()
plot3(Ni(:,1),Ni(:,2),Ni(:,3),'--r',Vrel_LVLH(:,1),Vrel_LVLH(:,2),Vrel_LVLH(:,3),'--b')
title('VELOCITY COMPARISON:!CW-INT! vs !CIRC-diff! in LVLH ref')
hold on
grid on 


figure()
plot3(Rrel_I(:,1),Rrel_I(:,2),Rrel_I(:,3),'--b',Rrel_I_rec(:,1),Rrel_I_rec(:,2),Rrel_I_rec(:,3),'--r',Rrel_I(1,1),Rrel_I(1,2),Rrel_I(1,3),'og',Rrel_I(end,1),Rrel_I(end,2),Rrel_I(end,3),'or',Rrel_I_rec(1,1),Rrel_I_rec(1,2),Rrel_I_rec(1,3),'og',Rrel_I_rec(end,1),Rrel_I_rec(end,2),Rrel_I_rec(end,3),'or')
title('!CW-INT! vs !CIRC-diff! in IN ref')
hold on
grid on 
figure()
plot3(Vrel_I(:,1),Vrel_I(:,2),Vrel_I(:,3),'--b',Vrel_I_rec(:,1),Vrel_I_rec(:,2),Vrel_I_rec(:,3),'--r')
title('VELOCITY COMPARISON:!CW-INT! vs !CIRC-diff! in IN ref')
hold on
grid on 

figure()
plot3(Rho(:,1),Rho(:,2),Rho(:,3),'b',RHO(:,1),RHO(:,2),RHO(:,3),'--k')
title('!CW-INT! vs !CW-AN!')
hold on
grid on 

figure()
plot3(Ni(:,1),Ni(:,2),Ni(:,3),'r',NI(:,1),NI(:,2),NI(:,3),'--k')
title('VELOCITY COMPARISON:!CW-INT! vs !CW-AN! ')
hold on
grid on 


