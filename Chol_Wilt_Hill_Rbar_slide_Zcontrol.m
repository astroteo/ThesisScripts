function   [x_dot]=Chol_Wilt_Hill_Rbar_slide_Zcontrol(t,x)

global W v_slide r_0 z_0;




r=x(1:3);
v=x(4:6);

r_dot(1)=v(1);
r_dot(2)=v(2);
r_dot(3)=v(3);

v_dot(1)=3*(W^2)*r(1)+2*W*v(2)-3*W^2*(r_0-v_slide*t);
v_dot(2)=-2*W*v(1)-2*W*v_slide;
v_dot(3)=-W^2*r(3)+W^2*z_0;
	
x_dot=[r_dot(1);r_dot(2);r_dot(3);v_dot(1);v_dot(2);v_dot(3)];




return