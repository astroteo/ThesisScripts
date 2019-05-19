function   [X_dot_t_c]=r2BP_t_c_PD_sat_traj1(t,X_t_c)


global kp kv f_sat m_cubesat y_SP0 v_slide l_bar;




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
rho=ROT'*(R_c-R_t);


omega_t=cross(R_t,V_t)/(norm(R_t)^2);

Ni=ROT'*((V_c-V_t)-cross(omega_t,(R_c-R_t)));

%trajectory

y_SP=y_SP0+v_slide*t;
    
    if y_SP >= 0  && y_SP < l_bar  
    v_slide_sp=v_slide;
    y_sp=y_sp+v_slide_sp*t;
   %disp('a')
    
   elseif y_SP > 0 && y_SP >= l_bar
    v_slide_sp=-v_slide;
    y_sp=y_sp+v_slide_sp*t;
    %disp('b')
   elseif y_SP <0 && y_SP >= (-l_bar)
    v_slide_sp=-v_slide;
    y_sp=y_sp+v_slide_sp*t;
    %disp('c')
    
    
   elseif y_SP <=0 && y_SP < (-l_bar)
    v_slide_sp=v_slide;
    y_sp=y_sp+v_slide*t;

    end
    
rho_SP=[0;y_sp;0];
Ni_SP=[0;v_slide_sp;0];
%PD
Kv=[kv 0 0;
    0 kv 0;
    0 0  kv];

Kp=[kp 0 0;
    0 kp 0;
    0 0  kp];


a_PD=ROT*(Kp*(rho-rho_SP)+Kv*(Ni-Ni_SP));
a_PD_SAT=a_PD;

for j=1:3
    
    if  a_PD_SAT(j) >  (f_sat/m_cubesat)
        
        a_PD_SAT(j)=f_sat/m_cubesat;
        
    elseif  a_PD_SAT(j) < -(f_sat/m_cubesat)
        
        a_PD_SAT(j)=-f_sat/m_cubesat;
        
   
    end
 end

F_PD_sat=m_cubesat*a_PD_SAT;

%chaser

R_abs_c=sqrt(R_c(1)^2+R_c(2)^2+R_c(3)^2);
R_dot_c=zeros(3,1);

R_dot_c(1)=V_c(1);
R_dot_c(2)=V_c(2);
R_dot_c(3)=V_c(3);

V_dot_c(1)=-(Mu/(R_abs_c^3)) *R_c(1)+a_PD_SAT(1);
V_dot_c(2)=-(Mu/(R_abs_c^3)) *R_c(2)+a_PD_SAT(2);
V_dot_c(3)=-(Mu/(R_abs_c^3)) *R_c(3)+a_PD_SAT(3);

X_dot_t_c=[R_dot_t(1); R_dot_t(2);R_dot_t(3);V_dot_t(1);V_dot_t(2);V_dot_t(3);R_dot_c(1); R_dot_c(2);R_dot_c(3);V_dot_c(1);V_dot_c(2);V_dot_c(3)];



return