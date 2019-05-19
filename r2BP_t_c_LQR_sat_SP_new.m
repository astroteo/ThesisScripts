function   [X_dot_t_c,t]=r2BP_t_c_LQR_sat_SP_new(t,X_t_c)

%FT--> final translation to the desired x_obs

global k_lqr  ROT_plane W w_v f_sat m_cubesat Set_Point x_obs i_Sp Horizon time_before;

i_v=[1;0;0];
j_v=[0;1;0];

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

%IN-->LVLH
   RT=norm(R_t);
   theta=-acos(dot(i_v,(R_t/RT)));

        if dot(j_v,(R_t/RT)) <0 
           theta=2*pi-theta;
        end

    ROT_theta=[ cos(theta), sin(theta),    0;
               -sin(theta), cos(theta),    0;
                 0       ,    0       ,    1];
                         

    ROT=ROT_plane*ROT_theta;
   
Rho=ROT'*(R_c-R_t);
Ni=ROT'*((V_c-V_t)-cross(w_v,(R_c-R_t)));

[N,~]=size(Set_Point);

if Horizon ==1 %SET-POINT: !pure trajectory! DON't WORK BECAUSE OF REDUCED ACTUATORS!

  for i=1:N
    
      if   Set_Point(i,1) > t
        
        i_up=i-1;
        Ni_Sp=Set_Point(i_up,5:7)';
        Rho_Sp=Set_Point(i_up,2:4)';
        
   
        
      else
        Ni_Sp=cross([0;0;W],(x_obs-Rho));%Ni_obs=cross([0;0;W],(Rho_Sp-Rho));
        Rho_Sp=x_obs;
%         Ni_Sp=zeros(3,1);
       
      end
    
  end

else %SET-POINT: ! Static Horizon !

Near_Rho=1;
Time_2_Change=200;

Rho_Sp=Set_Point(i_Sp,2:4)';
%Ni_Sp=Set_Point(i_Sp,5:7)'+cross([0;0;W],(Rho_Sp-Rho));
Ni_Sp=Set_Point(i_Sp,5:7)';

        
        if  i_Sp <= (N-1) 
            
             if norm(Rho-Rho_Sp)<= Near_Rho || (t-time_before) >= Time_2_Change
                i_Sp=i_Sp+1
                setGlobali_Sp(i_Sp)
                
                time_before=t;
                setGlobaltime_before(time_before);
             end
        else
            Rho_Sp=x_obs;
            Ni_Sp=cross([0;0;W],(x_obs-Rho));
           %Ni_Sp=zeros(3,1);
            
        
        end
    
end

%Actuators saturation

x=[Rho(1)-Rho_Sp(1);Ni(1)-Ni_Sp(1);Rho(2)-Rho_Sp(2);Ni(2)-Ni_Sp(2);Rho(3)-Rho_Sp(3);Ni(3)-Ni_Sp(3)];

a_LQR=ROT*(k_lqr*x);

a_LQR_SAT=a_LQR;

for j=1:3
    
         if  a_LQR(j) >  (f_sat/m_cubesat)
        
             a_LQR_SAT(j)=f_sat/m_cubesat;
             
             %disp('sat')
   
         elseif  a_LQR_SAT(j) < -(f_sat/m_cubesat)
        
                 a_LQR_SAT(j)=-f_sat/m_cubesat;
                 
                 %disp('sat')
         end
        
    
end

%Chaser integration 
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