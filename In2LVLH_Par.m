function [R_LVLH,V_LVLH]=In2LVLH_Par(R_In,V_In,a_t,e_t,i_t,om_t,OM_t,theta_t,Mu)

R_In=[R_In(1);R_In(2);R_In(3)];
V_In=[V_In(1);V_In(2);V_In(3)];

ROT_Plane=[cos(om_t)*cos(OM_t)-sin(om_t)*cos(i_t)*sin(OM_t), -sin(om_t)*cos(OM_t)-cos(om_t)*cos(i_t)*sin(OM_t), sin(i_t)*sin(OM_t);
    cos(om_t)*sin(OM_t)+sin(om_t)*cos(i_t)*cos(OM_t), -sin(om_t)*sin(OM_t)+cos(om_t)*cos(i_t)*cos(OM_t), -sin(i_t)*cos(OM_t);
    sin(om_t)*sin(i_t), cos(om_t)*sin(i_t), cos(i_t)]; 

ROT_Theta=[cos(theta_t), sin(theta_t),  0;
          -sin(theta_t), cos(theta_t),  0;
             0       ,     0     ,  1];
         


R_LVLH=ROT_Theta'* (ROT_Plane'* R_In);

[R_t,V_t]=ParOrb2RV(a_t,e_t,i_t,om_t,OM_t,theta_t,Mu);
hv=cross(R_t,V_t)/norm(cross(R_t,V_t));
w=sqrt(Mu/(a_t^3));
%V_LVLH=V_In+w*cross(hv,R_LVLH);
V_LVLH=V_In;%NON E' VERO !!!!
end