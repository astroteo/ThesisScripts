%contour
%             if y>=0
%                 
%                   if abs(x)<=lenght_ISS+step_grid && abs(x)>=lenght_ISS-step_grid  && abs(y) <= height_ISS  
%                 
%                       Phi_2d_contour(i,j)=1;
%                 
%                   elseif abs(x)<=lenght_ISS && abs(x) <=3/2*width_ISS &&  abs(y) <= height_ISS +step_grid && abs(y) >= height_ISS -step_grid
%                       
%                       Phi_2d_contour(i,j)=1;
%                 
%                   elseif (abs(x) <=lenght_ISS+step_grid && abs(x) >=lenght_ISS-step_grid) && abs(y)<=3/2*width_ISS && abs(y) >= height_ISS || ( abs(x) <=3/2*width_ISS+step_grid && abs(x) >= 3/2*width_ISS-step_grid) && abs(y)<=3/2*width_ISS && abs(y) >= height_ISS
%                      
%                   Phi_2d_contour(i,j)=1;
%                   
%                   elseif abs(x) <=lenght_ISS  &&  abs(x) >=3/2*width_ISS  && abs(y)<=3/2*width_ISS+step_grid && abs(y)>=3/2*width_ISS-step_grid
%                      
%                   Phi_2d_contour(i,j)=1;
%                   
%                   
%                   elseif abs(x) <=lenght_ISS && abs(x)>=height_ISS/2  && y<=step_grid 
%                      
%                   Phi_2d_contour(i,j)=1;
%                   
%                   
%                       
%                   end
%                 
%             elseif y<0
%                 
%                 
%                 if abs(x)<=height_ISS/2+step_grid && abs(x)>=height_ISS/2-step_grid  &&  abs(y) <= 2*width_ISS+step_grid
%                 
%                 Phi_2d_contour(i,j)=1;
%                 end
%                 
%                 if abs(x)<=height_ISS/2  &&  abs(y) <= 2*width_ISS+step_grid &&  abs(y) >= 2*width_ISS-step_grid 
%                 
%                 Phi_2d_contour(i,j)=1;
%                 end
%                 
%                
%             end
%             
%             if i==N || i==1 || j==M || j==1
%            
%            Phi_2d_contour(i,j)=1;
%             
%             end