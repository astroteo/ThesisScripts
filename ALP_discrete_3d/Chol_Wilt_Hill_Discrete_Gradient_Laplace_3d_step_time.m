function   [Rho,Ni,DV,count_impulse_time,count_impulse_done]=Chol_Wilt_Hill_Discrete_Gradient_Laplace_3d_step_time(dt,t_span,X_0)

global W  alpha_v  step_grid step_time Phi_3d L_box_x L_box_y L_box_z x_obs tau_or ;

[N_time,~]=size(t_span);

Rho=zeros(N_time,3);
Ni=zeros(N_time,3);

Rho(1,:)=X_0(1:3)';
Ni(1,:)=X_0(4:6)';



X_max=L_box_x;
Y_max=L_box_y;
Z_max=L_box_z;

R1_boundedx=-X_max:step_grid:X_max;
R1_boundedy=-Y_max:step_grid:Y_max;
R1_boundedz=-Z_max:step_grid:Z_max;

t=0;
tau=dt;

DV=zeros(N_time,3);

count_impulse_time=0;
count_impulse_done=0;
for i_t=2:1:N_time
    
t=t+i_t*dt;

r=Rho(i_t-1,1:3);
v=Ni(i_t-1,1:3);

x_t=r(1);
y_t=r(2);
z_t=r(3);

% if abs(x_t)< X_max && abs(y_t) < Y_max && abs(z_t)< Z_max

i_bottom=ceil((x_t+X_max)/step_grid);
i_top=i_bottom+1;


        
j_bottom=ceil((y_t+Y_max)/step_grid);
j_top=j_bottom+1;


k_bottom=ceil((z_t+Z_max)/step_grid);
k_top=k_bottom+1;

% 2D recap:
% [~,PP]=size(R1_boundedz);
% P=ceil(PP/2);

% fx=(1/step_grid)^2*( -Phi_3d(i_bottom,j_bottom,P)*(R1_boundedy(j_top)-y_t) + Phi_3d(i_top,j_bottom,P)*( R1_boundedy(j_top)-y_t)  - Phi_3d(i_bottom,j_top,P)*(y_t-R1_boundedy(j_bottom))+ Phi_3d(i_top,j_top,P)*(y_t-R1_boundedy(j_bottom)));
% fy=(1/step_grid)^2*( -Phi_3d(i_bottom,j_bottom,P)*(R1_boundedx(i_top)-x_t) - Phi_3d(i_top,j_bottom,P)*(x_t-R1_boundedx(i_bottom))+ Phi_3d(i_bottom,j_top,P)*(R1_boundedx(i_top)-x_t)   + Phi_3d(i_top,j_top,P)*(x_t-R1_boundedx(i_bottom)));
% fz=0;

%manually computed gradient
% fx=(1/step_grid)^3*( Phi_3d(i_bottom,j_bottom,k_bottom)*(y_t-R1_boundedy(j_top))*(R1_boundedz(k_top)-z_t) + Phi_3d(i_top,j_bottom,k_bottom)*( R1_boundedy(j_top)-y_t)*(R1_boundedz(k_top)-z_t)    + Phi_3d(i_bottom,j_top,k_bottom)*(R1_boundedy(j_bottom)-y_t)*(R1_boundedz(k_top)-z_t) + Phi_3d(i_top,j_top,k_bottom)*(y_t-R1_boundedy(j_bottom))*(R1_boundedz(k_top)-z_t)   + Phi_3d(i_bottom,j_bottom,k_top)*(y_t-R1_boundedy(j_top))*(z_t-R1_boundedz(k_bottom)) + Phi_3d(i_top,j_bottom,k_top)*(R1_boundedy(j_top)-y_t)*(z_t-R1_boundedz(k_bottom))   + Phi_3d(i_bottom,j_top,k_top)*(R1_boundedy(j_bottom)-y_t)*(z_t-R1_boundedz(k_bottom)) + Phi_3d(i_top,j_top,k_top)*(y_t-R1_boundedy(j_bottom))*(z_t-R1_boundedz(k_bottom)));
% fy=(1/step_grid)^3*( Phi_3d(i_bottom,j_bottom,k_bottom)*(x_t-R1_boundedx(i_top))*(R1_boundedz(k_top)-z_t) + Phi_3d(i_top,j_bottom,k_bottom)*( R1_boundedx(i_bottom)-x_t)*(R1_boundedz(k_top)-z_t) + Phi_3d(i_bottom,j_top,k_bottom)*(R1_boundedx(i_top)-x_t)*(R1_boundedz(k_top)-z_t)    + Phi_3d(i_top,j_top,k_bottom)*(x_t-R1_boundedx(i_bottom))*(R1_boundedz(k_top)-z_t)   + Phi_3d(i_bottom,j_bottom,k_top)*(x_t-R1_boundedx(i_top))*(z_t-R1_boundedz(k_bottom)) + Phi_3d(i_top,j_bottom,k_top)*(R1_boundedx(i_bottom)-x_t)*(z_t-R1_boundedz(k_bottom))+ Phi_3d(i_bottom,j_top,k_top)*(R1_boundedx(i_top)-x_t)*(z_t-R1_boundedz(k_bottom))    + Phi_3d(i_top,j_top,k_top)*(x_t-R1_boundedx(i_bottom))*(z_t-R1_boundedz(k_bottom)));
% fz=(1/step_grid)^3*( Phi_3d(i_bottom,j_bottom,k_bottom)*(R1_boundedy(j_top)-y_t)*(R1_boundedx(i_top)-x_t) + Phi_3d(i_top,j_bottom,k_bottom)*( R1_boundedy(j_top)-y_t)*(x_t-R1_boundedx(i_bottom)) + Phi_3d(i_bottom,j_top,k_bottom)*(y_t-R1_boundedy(j_bottom))*(R1_boundedx(i_top)-x_t) + Phi_3d(i_top,j_top,k_bottom)*(y_t-R1_boundedy(j_bottom))*(x_t-R1_boundedx(i_bottom))+ Phi_3d(i_bottom,j_bottom,k_top)*(R1_boundedy(j_top)-y_t)*(x_t-R1_boundedx(i_top))    + Phi_3d(i_top,j_bottom,k_top)*(R1_boundedy(j_top)-y_t)*(R1_boundedx(i_bottom)-x_t)   + Phi_3d(i_bottom,j_top,k_top)*(y_t-R1_boundedy(j_bottom))*(x_t-R1_boundedx(i_top))    + Phi_3d(i_top,j_top,k_top)*(y_t-R1_boundedy(j_bottom))*(x_t-R1_boundedx(i_bottom)));

% mathematica computed gradient
Xbottom= R1_boundedx(i_bottom);
Xtop=R1_boundedx(i_top);

Ybottom= R1_boundedy(j_bottom);
Ytop=R1_boundedy(j_top);

Zbottom= R1_boundedz(k_bottom);
Ztop=R1_boundedz(k_top);

% fz=(1/step_grid)^3*(Phi_3d(i_top,j_top,k_bottom)*(x_t - Xbottom)*(y_t - Ybottom) + Phi_3d(i_top,j_top,k_top)*(x_t - Xbottom)*(y_t - Ybottom) - Phi_3d(i_bottom,j_top,k_bottom)*(-x_t + Xtop)*(y_t - Ybottom) + Phi_3d(i_bottom,j_top,k_top)*(-x_t + Xtop)*(y_t - Ybottom) - Phi_3d(i_top,j_bottom,k_bottom)*(x_t - Xbottom)*(-y_t + Ytop) + Phi_3d(i_top,j_bottom,k_top)*(x_t - Xbottom)*(-y_t + Ytop) - Phi_3d(i_bottom,j_bottom,k_bottom)*(-x_t + Xtop)*(-y_t + Ytop) + Phi_3d(i_bottom,j_bottom,k_top)*(-x_t + Xtop)*(-y_t + Ytop));
% fx=(1/step_grid)^3*(Phi_3d(i_top,j_top,k_top)*(y_t - Ybottom)*(z_t - Zbottom) - Phi_3d(i_bottom,j_top,k_top)*(y_t - Ybottom)*(z_t - Zbottom) + Phi_3d(i_top,j_bottom,k_top)*(-y_t + Ytop)*(z_t - Zbottom) - Phi_3d(i_bottom,j_bottom,k_top)*(-y_t + Ytop)*(z_t - Zbottom) + Phi_3d(i_top,j_top,k_bottom)*(y_t - Ybottom)*(-z_t + Ztop) - Phi_3d(i_bottom,j_top,k_bottom)*(y_t - Ybottom)*(-z_t + Ztop) - Phi_3d(i_bottom,j_bottom,k_bottom)*(-y_t + Ytop)*(-z_t + Ztop) + Phi_3d(i_top,j_bottom,k_bottom)*(-y_t + Ytop)*(-z_t + Ztop));
% fy=(1/step_grid)^3*(Phi_3d(i_top,j_top,k_top)*(x_t - Xbottom)*(z_t - Zbottom) - Phi_3d(i_top,j_bottom,k_top)*(x_t - Xbottom)*(z_t - Zbottom) + Phi_3d(i_bottom,j_top,k_top)*(-x_t + Xtop)*(z_t - Zbottom) - Phi_3d(i_bottom,j_bottom,k_top)*(-x_t + Xtop)*(z_t - Zbottom) - Phi_3d(i_top,j_bottom,k_bottom)*(x_t - Xbottom)*(-z_t + Ztop) + Phi_3d(i_top,j_top,k_bottom)*(x_t - Xbottom)*(-z_t + Ztop) - Phi_3d(i_bottom,j_bottom,k_bottom)* (-x_t + Xtop)*(-z_t + Ztop) + Phi_3d(i_bottom,j_top,k_bottom)*(-x_t + Xtop)*(-z_t + Ztop));
K=Phi_3d(i_bottom,j_bottom,k_bottom);
Ki=Phi_3d(i_top,j_bottom,k_bottom);
Kj=Phi_3d(i_bottom,j_top,k_bottom);
Kk=Phi_3d(i_bottom,j_bottom,k_top);
Kij=Phi_3d(i_top,j_top,k_bottom);
Kik=Phi_3d(i_top,j_bottom,k_top);
Kjk=Phi_3d(i_bottom,j_top,k_top);
Kijk=Phi_3d(i_top,j_top,k_top);
x=x_t;
y=y_t;
z=z_t;
  %1) no gradient felt
% fz=-Kij*(x - Xbottom)*(y - Ybottom) + Kijk*(x - Xbottom)*(y - Ybottom) -  Kj*(-x + Xtop)*(y - Ybottom) + Kjk*(-x + Xtop)*(y - Ybottom) -  Ki*(x - Xbottom)*(-y + Ytop) + Kik*(x - Xbottom)*(-y + Ytop) - K*(-x + Xtop)*(-y + Ytop) + Kk*(-x + Xtop)*(-y + Ytop);
% fx=Kijk*(y - Ybottom)*(z - Zbottom) - Kjk*(y - Ybottom)*(z - Zbottom) + Kik*(-y + Ytop)*(z - Zbottom) - Kk*(-y + Ytop)*(z - Zbottom) + Kij*(y - Ybottom)*(-z + Ztop) - Kj*(y - Ybottom)*(-z + Ztop) - K*(-y + Ytop)*(-z + Ztop) + Ki*(-y + Ytop)*(-z + Ztop);
% fy=Kijk*(x - Xbottom)*(z - Zbottom) - Kik*(x - Xbottom)*(z - Zbottom) + Kjk*(-x + Xtop)*(z - Zbottom) - Kk*(-x + Xtop)*(z - Zbottom) - Ki*(x - Xbottom)*(-z + Ztop) + Kij*(x - Xbottom)*(-z + Ztop) - K*(-x + Xtop)*(-z + Ztop) + Kj*(-x + Xtop)*(-z + Ztop);

%2) 

fx= (1/step_grid)^3*(Kijk*(y - Ybottom)*(z - Zbottom) - Kjk*(y - Ybottom)*(z - Zbottom) + Kk*(y - Ytop)*(z - Zbottom) + Kik*(-y + Ytop)*(z - Zbottom) + Kj*(y - Ybottom)*(z - Ztop) + Ki*(y - Ytop)*(z - Ztop) + Kij*(y - Ybottom)*(-z + Ztop) - K*(-y + Ytop)*(-z + Ztop));
fy= (1/step_grid)^3*(Kijk*(x - Xbottom)*(z - Zbottom) - Kik*(x - Xbottom)*(z - Zbottom) + Kk*(x - Xtop)*(z - Zbottom) + Kjk*(-x + Xtop)*(z - Zbottom) + Ki*(x - Xbottom)*(z - Ztop) + Kj*(x - Xtop)*(z - Ztop) + Kij*(x - Xbottom)*(-z + Ztop) - K*(-x + Xtop)*(-z + Ztop));
fz= (1/step_grid)^3*(Kijk*(x - Xbottom)*(y - Ybottom) - Kij*(x - Xbottom)*(y - Ybottom) + Kj*(x - Xtop)*(y - Ybottom) + Kjk*(-x + Xtop)*(y - Ybottom) + Ki*(x - Xbottom)*(y - Ytop) + Kk*(x - Xtop)*(y - Ytop) + Kik*(x - Xbottom)*(-y + Ytop) - K*(-x + Xtop)*(-y + Ytop));
f=[fx;fy;fz];

 if t-floor(t/step_time)*step_time==0
%     
   
 count_impulse_time=count_impulse_time+1;

if dot(f,v) >=0 
    
    disp('fight_gradient')

 Kv=1*alpha_v; 

V=-Kv*f/norm(f);

Ni_old=Ni(i_t-1,:);
Ni(i_t-1,:)=[V(1),V(2),V(3)];
DV(i_t-1,:)=Ni(i_t-1,:)-Ni_old;

count_impulse_done=count_impulse_done+1;

else
    
Ni(i_t-1,:)=[v(1),v(2),v(3)];


end

% Ni(i_t-1,:)=[v(1),v(2),v(3)];

 end 
 
% else
%     disp('we are out of the grid')
%     Ni(i_t-1,:)=[v(1),v(2),v(3)];

    
% Ct=cos(W*tau_or);
% St=sin(W*tau_or);    
%     
% det=(4*St)/(W^3) -(8*Ct*St)/(W^3) +(4*Ct^2*St)/(W^3)+(4*St^3)/(W^3) -(3*St^2*tau_or)/(W^2);
% 
% N_tau_inv=1/det*[(4*St^2)/(W^2)-(3*St*tau_or)/W,     -((2*St)/(W^2))+(2*Ct*St)/(W^2),                        0;
%                  (2*St)/(W^2)-(2*Ct*St)/(W^2),              St^2/(W^2),                                   0;
%                            0,                                  0,              4/(W^2)-(8*Ct)/(W^2)+(4*Ct^2)/(W^2)+(4*St^2)/(W^2)-(3*St*tau_or)/W];
%                    
% M_tau=[-3*Ct+4,        0,   0;
%         6*St-6*W*tau_or,  1,   0;                   
%                0    ,  0,  Ct];
% 
% 
% Ni(i_t-1,:)=(N_tau_inv*(x_obs-M_tau*r'))';

    
    
% end

Ct=cos(W*tau);
St=sin(W*tau);

Rho_Ni_0=[Rho(i_t-1,:)';Ni(i_t-1,:)'];

PHI=[4-3*Ct       0      0       St/W          2/W*(1-Ct)          0;
    6*(St-W*tau)  1      0    -2/W*(1-Ct)    1/W*(4*St-3*W*tau)    0;
        0         0      Ct       0                0           1/W*St;
     3*W*St       0      0        Ct             2*St              0;
   -6*W*(1-Ct)    0      0      -2*St           4*Ct-3             0;
        0         0    -W*St      0                0              Ct ];
    
Rho_Ni=PHI*Rho_Ni_0;
    
Rho(i_t,:)=Rho_Ni(1:3)';
Ni(i_t,:)=Rho_Ni(4:6)';

% 
% 
% 
end

return

