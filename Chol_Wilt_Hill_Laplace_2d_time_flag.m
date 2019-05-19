function   [x_dot]=Chol_Wilt_Hill_Laplace_2d_time_flag(~,x)

global W phi_2d step_grid step_time total_stepy  total_stepx  alpha_v x_obs_2d Rho_0 flag_time;

r=x(1:2);
v=x(3:4);


%compute grad(phi_2d)
x=r(1);
y=r(2);

if x >=0
i_bottom=total_stepx/2+floor(x/step_grid);
i_up=i_bottom+1;

xd=(x-i_bottom*step_grid);

    
else     
i_up=ceil((total_stepx/2*step_grid-abs(x))/step_grid);
i_bottom=i_up-1;

xd=x+(total_stepx/2-i_up)*step_grid;

end


if y >=0
j_bottom=total_stepy/2+floor(y/step_grid);
j_up=j_bottom+1;

yd=(y-j_bottom*step_grid);

    
else     
j_up=ceil((total_stepy/2*step_grid-abs(y))/step_grid);
j_bottom=j_up-1;

yd=y+(total_stepy/2-j_up)*step_grid;


    
end




%phi_c= phi_2d(i_bottom,j_bottom)*(1-xd)*(1-yd) + phi_2d(i_up,j_bottom)*xd*(1-yd) + phi_2d(i_bottom,j_up)*(1-xd)*yd + phi_2d(i_up,j_up)*xd*yd; 

fx=-phi_2d(i_bottom,j_bottom)*(1-yd)  + phi_2d(i_up,j_bottom)*(1-yd) - phi_2d(i_bottom,j_up)*yd     + phi_2d(i_up,j_up)*yd;
fy=-phi_2d(i_bottom,j_bottom)*(1-xd)  - phi_2d(i_up,j_bottom)*xd     + phi_2d(i_bottom,j_up)*(1-xd) + phi_2d(i_up,j_up)*xd;

f=[fx;fy];



if (t-floor(t/step_time)*step_time) <= step_time && flag_time==0
    
    flag_time=1;
    setGlobalflag_time(flag_time);
else
    
    flag_time=0;
    setGlobalflag_time(flag_time);

end



if dot(f,v) >=0 



K=2*alpha_v*flag_time;

V=-K*f/norm(f);
v=V;
%disp('corr');
 
end


  if abs(r(1)) > total_stepx*step_grid || abs(r(2)) > total_stepx*step_grid
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