function   [flag_coll]=Evaluate_collision(Rho,Surf,step_grid)

[N_pos,~]=size(Rho);

flag_coll=0;

[N_surf,M_surf,~]=size(Surf);

X_max=floor(N_surf/2)*step_grid;
Y_max=floor(M_surf/2)*step_grid;


   for i=1:N_pos
       
       x_p=Rho(i,1);
       y_p=Rho(i,2);
      
     if abs(x_p) <= X_max && abs(y_p) <= Y_max
       
       if x_p >=0
       i_p=floor(x_p/step_grid)+floor((X_max/step_grid));
       else
       i_p=floor((x_p+floor(X_max/step_grid)*step_grid)/step_grid);
       end

      if y_p >=0
       j_p=floor(y_p/step_grid)+floor((Y_max/step_grid));
      else
       j_p=floor((y_p+floor(Y_max/step_grid)*step_grid)/step_grid);
      end
      
      if i_p==0 || j_p==0
         
           flag_coll=1;
       
     elseif i_p~=0 && j_p~=0 && Surf(i_p,j_p) ~= 0 
           
           flag_coll=1;
       
          
      end
       
     else
         
         flag_coll=1;
     end
       
     
       
       
      
       
   end





return