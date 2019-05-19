function   [Rho,Ni,DV,count_impulse_time,count_impulse_done,i_obs]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_2d_step_time(dt,t_span,X_0,x_obs)

global W  alpha_v  step_grid step_time Phi_2d L_box_x L_box_y  ;

[N_time,~]=size(t_span);

Rho=zeros(N_time,3);
Ni=zeros(N_time,3);

i_obs=1


Rho(1,:)=X_0(1:3)';
Ni(1,:)=X_0(4:6)';

%X_max=4*abs(X_0(1));
%Y_max=4*abs(X_0(2));

X_max=L_box_x;
Y_max=L_box_y;

R1_boundedx=-X_max:step_grid:X_max;
R1_boundedy=-Y_max:step_grid:Y_max;

t=0;
tau=dt;

DV=zeros(N_time,3);

count_impulse_time=0;
count_impulse_done=0;
for i=2:1:N_time
    
t=t+i*dt;

r=Rho(i-1,1:3);
v=Ni(i-1,1:3);

x_t=r(1);
y_t=r(2);

i_bottom=ceil((x_t+X_max)/step_grid);
i_top=i_bottom+1;
        
        
j_bottom=ceil((y_t+Y_max)/step_grid);
j_top=j_bottom+1;
      



fx=(1/step_grid)^2*( -Phi_2d(i_bottom,j_bottom)*(R1_boundedy(j_top)-y_t) + Phi_2d(i_top,j_bottom)*( R1_boundedy(j_top)-y_t)  - Phi_2d(i_bottom,j_top)*(y_t-R1_boundedy(j_bottom))+ Phi_2d(i_top,j_top)*(y_t-R1_boundedy(j_bottom)));
fy=(1/step_grid)^2*( -Phi_2d(i_bottom,j_bottom)*(R1_boundedx(i_top)-x_t) - Phi_2d(i_top,j_bottom)*(x_t-R1_boundedx(i_bottom))+ Phi_2d(i_bottom,j_top)*(R1_boundedx(i_top)-x_t)   + Phi_2d(i_top,j_top)*(x_t-R1_boundedx(i_bottom)));


f=[fx;fy;0];

 if t-floor(t/step_time)*step_time==0
%     
   
 count_impulse_time=count_impulse_time+1;

if dot(f,v) >=0 
    
    %disp('fight_gradient')

 K=1*alpha_v; 

V=-K*f/norm(f);

Ni_old=Ni(i-1,:);
Ni(i-1,:)=[V(1),V(2),V(3)];
DV(i-1,:)=Ni(i-1,:)-Ni_old;

count_impulse_done=count_impulse_done+1;

else
    
Ni(i-1,:)=[v(1),v(2),v(3)];

end

 end  

Ct=cos(W*tau);
St=sin(W*tau);

Rho_Ni_0=[Rho(i-1,:)';Ni(i-1,:)'];

PHI=[4-3*Ct       0      0       St/W          2/W*(1-Ct)          0;
    6*(St-W*tau)  1      0    -2/W*(1-Ct)    1/W*(4*St-3*W*tau)    0;
        0         0      Ct       0                0           1/W*St;
     3*W*St       0      0        Ct             2*St              0;
   -6*W*(1-Ct)    0      0      -2*St           4*Ct-3             0;
        0         0    -W*St      0                0              Ct ];
    
Rho_Ni=PHI*Rho_Ni_0;
    
Rho(i,:)=Rho_Ni(1:3)';
Ni(i,:)=Rho_Ni(4:6)';


% if norm (Rho(i,:)-x_obs') <= 4*step_grid
%     
%     i_obs=i;  
% end

if isnan(Rho)
    disp('Rho_nan')
elseif  isnan(Ni)
     disp('Ni_nan')
end

end

return

