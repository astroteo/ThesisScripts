function   [x_dot]=Chol_Wilt_Hill_Laplace_2d_analytic(~,x)

global W  alpha_v M_2d N_2d lambda_v l_max x_obs_2d;


r=x(1:2);
v=x(3:4);


%parameters recovering

m1=M_2d(1,1);
m2=M_2d(2,2);

n1=N_2d(1,1);
n2=N_2d(2,2);

x=r(1);
y=r(2);

x_obs=x_obs_2d(1);
y_obs=x_obs_2d(2);

lambda_1=lambda_v(1);
lambda_2=lambda_v(2);

% gradient analytic computation:

%phi_harm_2d(i,j)=lambda_1*exp(-lambda_2*(X_2d)'*N_2d*(X_2d));
f_harm_x=lambda_1*(2*n1*(-lambda_2)*(x/l_max)*exp((-lambda_2)*(n1*(x/l_max)^2+n2*(y/l_max)^2)));
f_harm_y=lambda_1*(2*n2*(-lambda_2)*(y/l_max)*exp((-lambda_2)*(n1*(x/l_max)^2+n2*(y/l_max)^2)));

%phi_2d(i,j)=1/2*(X_2d/L_max-X_obs_2d/L_max)'*M_2d*(X_2d/L_max-X_obs_2d/L_max);
f_pos_x=m1*(x/l_max-x_obs/l_max);
f_pos_y=m2*(y/l_max-y_obs/l_max);

%d/dx(f(x,y)+g(x,y))=df(x,y)/dx+dg(x,y)/dx 
%d/dy(f(x,y)+g(x,y))=df(x,y)/dy+dg(x,y)/dy
fx=f_harm_x+f_pos_x;
fy=f_harm_y+f_pos_y;

f=[fx;fy];

if dot(f,v) >=0 

K=1*alpha_v; %2*alpha_v ref.KOREA

V=-K*f/norm(f);
v=V;
%disp('corr');
 
end


  if abs(r(1)) > l_max || abs(r(2)) > l_max
      disp('we are out of reticulate')
      x
      y
  end
    
r_dot(1)=v(1);
r_dot(2)=v(2);


v_dot(1)=3*(W^2)*r(1)+2*W*v(2);
v_dot(2)=-2*W*v(1);

	
x_dot=[r_dot(1);r_dot(2);v_dot(1);v_dot(2)];




return