function   [X_dot_t_c]=r2BP_t_c_LQR_sat_SK(t,X_t_c)

%FT--> final translation to the desired x_obs

global k_lqr f_sat m_cubesat x_obs  kp kv ;

contol_method=1; %<-- 1 to use LQR, 0 to use PD

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

[~,~,i_t,OM_t,om_t,theta_t]=RV2OrbPar_My(R_t,V_t,Mu);


ROT_plane=[cos(om_t)*cos(OM_t)-sin(om_t)*cos(i_t)*sin(OM_t) -sin(om_t)*cos(OM_t)-cos(om_t)*cos(i_t)*sin(OM_t) sin(i_t)*sin(OM_t)
    cos(om_t)*sin(OM_t)+sin(om_t)*cos(i_t)*cos(OM_t) -sin(om_t)*sin(OM_t)+cos(om_t)*cos(i_t)*cos(OM_t) -sin(i_t)*cos(OM_t)
    sin(om_t)*sin(i_t) cos(om_t)*sin(i_t) cos(i_t)]  ;

ROT_theta=[ cos(theta_t), sin(theta_t),0;
           -sin(theta_t), cos(theta_t),0;
                 0       ,    0       ,1];
             
ROT=ROT_plane*ROT_theta;
Rho=ROT'*(R_c-R_t);


omega_t=cross(R_t,V_t)/(norm(R_t)^2);

Ni=ROT'*((V_c-V_t)-cross(omega_t,(R_c-R_t)));

Ni_obs=[0; -2*norm(omega_t)*x_obs(1); 0 ];

%chaser
x=[Rho(1)-x_obs(1);Ni(1)-Ni_obs(1);Rho(2)-x_obs(2);Ni(2)-Ni_obs(2);Rho(3)-x_obs(3);Ni(3)--Ni_obs(3)];

if contol_method == 1
a_LQR=ROT*(k_lqr*x);

a_LQR_SAT=a_LQR;
% for j=1:3
%     
%     if  a_LQR_SAT(j) >  (f_sat/m_cubesat)
%         
%         a_LQR_SAT(j)=f_sat/m_cubesat;
%    
%     elseif  a_LQR_SAT(j) < -(f_sat/m_cubesat)
%         
%         a_LQR_SAT(j)=-f_sat/m_cubesat;
%     end
%         
%     
% end

else
Kv=[kv 0 0;
    0 kv 0;
    0 0  kv];

Kp=[kp 0 0;
    0 kp 0;
    0 0  kp];
    
a_LQR_SAT=ROT*(Kp*(Rho-x_obs)+Kv*Ni);    
    
end


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