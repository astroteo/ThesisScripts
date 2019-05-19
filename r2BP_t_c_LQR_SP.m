function   [X_dot_t_c,t]=r2BP_t_c_LQR_SP(t,X_t_c)

%FT--> final translation to the desired x_obs

global k_lqr x_obs ROTplane_0, w_v, w ;



Mu=3.9856e+14;

R_t=X_t_c(1:3);
V_t=X_t_c(4:6);

R_c=X_t_c(7:9);
V_c=X_t_c(10:12);



%target
R_abs_t=sqrt(R_t(1)^2+R_t(2)^2+R_t(3)^2);

R_dot_t=zeros(3,1);

R_dot_t(1)=V_t(1);
R_dot_t(2)=V_t(2);
R_dot_t(3)=V_t(3);

V_dot_t(1)=-(Mu/(R_abs_t^3)) *R_t(1);
V_dot_t(2)=-(Mu/(R_abs_t^3)) *R_t(2);
V_dot_t(3)=-(Mu/(R_abs_t^3)) *R_t(3);


    
    theta=-w*t;%!! TAKE CARE FOR THIS SIGN !!
    ROT_theta=[ cos(theta), sin(theta),    0;
               -sin(theta), cos(theta),    0;
                 0       ,    0       ,    1];

    ROT=ROTplane_0*ROT_theta;
   
   
Rho=ROT'*(R_c-R_t);



Ni=ROT'*((V_c-V_t)-cross(omega_t,(R_c-R_t)));
Ni_obs=cross(omega_t,(x_obs-Rho));


%chaser

x=[Rho(1)-x_obs(1);Ni(1)-Ni_obs(1);Rho(2)-x_obs(2);Ni(2)-Ni_obs(2);Rho(3)-x_obs(3);Ni(3)-Ni_obs(3)];

a_LQR=ROT*(k_lqr*x);

a_LQR_SAT=a_LQR;


R_abs_c=sqrt(R_c(1)^2+R_c(2)^2+R_c(3)^2);

R_dot_c=zeros(3,1);


R_dot_c(1)=V_c(1);
R_dot_c(2)=V_c(2);
R_dot_c(3)=V_c(3);

V_dot_c(1)=-(Mu/(R_abs_c^3)) *R_c(1)+a_LQR_SAT(1);
V_dot_c(2)=-(Mu/(R_abs_c^3)) *R_c(2)+a_LQR_SAT(2);
V_dot_c(3)=-(Mu/(R_abs_c^3)) *R_c(3)+a_LQR_SAT(3);

X_dot_t_c=[R_dot_t(1); R_dot_t(2);R_dot_t(3);V_dot_t(1);V_dot_t(2);V_dot_t(3);R_dot_c(1); R_dot_c(2);R_dot_c(3);V_dot_c(1);V_dot_c(2);V_dot_c(3)];



return