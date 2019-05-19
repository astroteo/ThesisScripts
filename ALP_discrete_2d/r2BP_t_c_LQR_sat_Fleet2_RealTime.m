function   [X_dot_t_c,t]=r2BP_t_c_LQR_sat_Fleet2_RealTime(t,X_t_c1_c2)

%FT--> final translation to the desired x_obs

global k_lqr  ROT_plane W w_v f_sat m_cubesat x_obs1 x_obs2 i_Sp Horizon   M_2d N_2d N2d_f lambda_v l_max x_obs_2d flag_time step_time time_before time_last_fire flag_fire;

i_v=[1;0;0];
j_v=[0;1;0];

Mu=3.9856e+14;

R_t=X_t_c1_c2(1:3);
V_t=X_t_c1_c2(4:6);

R_c1=X_t_c1_c2(7:9);
V_c1=X_t_c1_c2(10:12);


R_c2=X_t_c1_c2(7:9);
V_c2=X_t_c1_c2(10:12);








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
   
Rho1=ROT'*(R_c1-R_t);
Ni1=ROT'*((V_c1-V_t)-cross(w_v,(R_c1-R_t)));

Rho2=ROT'*(R_c2-R_t);
Ni2=ROT'*((V_c2-V_t)-cross(w_v,(R_c2-R_t)));

[N,~]=size(Set_Point);



%Set Point generation (ALP--> Analitycal gradient)

m1=M_2d(1,1);
m2=M_2d(2,2);

n1=N_2d(1,1);
n2=N_2d(2,2);

x1=Rho1(1);
y1=Rho1(2);

f_harm_x1=lambda_1*(2*n1*(lambda_2)*(x1/l_max)*exp((-lambda_2)*(n1*(x1/l_max)^2+n2*(y1/l_max)^2)))+lambda_1*(2*n1*(lambda_2)*(x1-Rho2(1)/l_max)*exp((-lambda_2)*(n1*(x1-Rho2(1))/l_max)^2+n2*((y1-Rho2(2)/l_max)^2)));
f_harm_y1=lambda_1*(2*n2*(lambda_2)*(y1/l_max)*exp((-lambda_2)*(n1*(x1/l_max)^2+n2*(y1/l_max)^2)))+lambda_1*(2*n2*(lambda_2)*(y/l_max)*exp((-lambda_2)*(n1*(x1/l_max)^2+n2*(y1/l_max)^2)));

%phi_2d(i,j)=1/2*(X_2d/L_max-X_obs_2d/L_max)'*M_2d*(X_2d/L_max-X_obs_2d/L_max);
f_pos_x1=m1*(x1/l_max-x_obs1(1,1)/l_max);
f_pos_y1=m2*(y1/l_max-x_obs1(2,1)/l_max);

fx1=f_harm_x1+f_pos_x1;
fy1=f_harm_y1+f_pos_y1;

f1=[fx1;fy1];



% coorection at discrete time:

if  (t-time_before) <= step_time

  if (t-time_last_fire) <= step_time 
    
      
    flag_time=1;
    setGlobalflag_fire(flag_time);
    
     
  else
    
    flag_time=0; 
    setGlobalflag_time(flag_time);
    
    if flag_time ==1  && flag_fire==0;
    
    flag_fire=1;
    setGlobalflag_fire(flag_fire);
    
    time_last_fire=t;
    setGlobaltime_last_fire(time_last_fire);
    end
  end

else
  
    flag_time=1;
end


if dot(f,v) >=0 

K=2*alpha_v*flag_time*flag_fire; %2*alpha_v ref.KOREA

% if K ~= 0
%     disp ('corr')
% end

V=-K*f/norm(f);
v=V;
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

Near_Rho=0.2;
Time_2_Change=100;

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


%Control Action
x=[Rho(1)-Rho_Sp(1);Ni(1)-Ni_Sp(1);Rho(2)-Rho_Sp(2);Ni(2)-Ni_Sp(2);Rho(3)-Rho_Sp(3);Ni(3)-Ni_Sp(3)];
a_LQR=ROT*(k_lqr*x);
a_LQR_SAT=a_LQR;

%Actuators saturation
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