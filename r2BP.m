function   [X_dot]=r2BP(t,X)
	
Mu=3.9856e+14;
R=X(1:3);
V=X(4:6);
R_abs=sqrt(R(1)^2+R(2)^2+R(3)^2);

R_dot(1)=V(1);
R_dot(2)=V(2);
R_dot(3)=V(3);

V_dot(1)=-(Mu/(R_abs^3)) *R(1);
V_dot(2)=-(Mu/(R_abs^3)) *R(2);
V_dot(3)=-(Mu/(R_abs^3)) *R(3);



X_dot=[R_dot(1); R_dot(2);R_dot(3);V_dot(1);V_dot(2);V_dot(3)];

return