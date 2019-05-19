function [R_LVLH,V_LVLH]=In2LVLH_rel(R_In,V_In,a,i,om,OM,theta,Mu)


R_In=[R_In(1);R_In(2);R_In(3)];

ROT_Plane=[cos(om)*cos(OM)-sin(om)*cos(i)*sin(OM), -sin(om)*cos(OM)-cos(om)*cos(i)*sin(OM), sin(i)*sin(OM);
    cos(om)*sin(OM)+sin(om)*cos(i)*cos(OM), -sin(om)*sin(OM)+cos(om)*cos(i)*cos(OM), -sin(i)*cos(OM);
    sin(om)*sin(i), cos(om)*sin(i), cos(i)]; 

ROT_Theta=[cos(theta), sin(theta),  0;
          -sin(theta), cos(theta), 0;
             0       ,     0     ,1 ];
         


R_LVLH=ROT_Theta'* (ROT_Plane'* R_In);

hv=cross(R_In,V_In)/norm(cross(R_In,V_In));
w=sqrt(Mu/(a^3));
V_LVLH=V_In+w*cross(hv,R_LVLH);
