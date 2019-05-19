function [R_In,V_In]=LVLH2In(R_LVLH,V_LVLH,i,OM,om,theta)

Mu=3.9856e+14;

R_LVLH=[R_LVLH(1);R_LVLH(2);R_LVLH(3)];

ROT_Plane=[cos(om)*cos(OM)-sin(om)*cos(i)*sin(OM), -sin(om)*cos(OM)-cos(om)*cos(i)*sin(OM), sin(i)*sin(OM);
    cos(om)*sin(OM)+sin(om)*cos(i)*cos(OM), -sin(om)*sin(OM)+cos(om)*cos(i)*cos(OM), -sin(i)*cos(OM);
    sin(om)*sin(i), cos(om)*sin(i), cos(i)]; 

ROT_Theta=[cos(theta), sin(theta),  0;
          -sin(theta), cos(theta), 0;
             0       ,     0     ,1 ];
         
 

R_In=ROT_Plane* (ROT_Theta* R_LVLH);
V_In=V_LVLH;%NON E' VERO !!!!!
end