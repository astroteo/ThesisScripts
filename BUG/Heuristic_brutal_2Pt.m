function [Rho_0_new,n_t]=Heuristic_brutal_2Pt(Rho_m_star,O_k,R_gnc,Surf,step_grid)


%[K,~]=size(O_k);
%d_min=norm(Rho_m_star-O_k(1)')+norm(X_goal-O_k(1));

n_t=(O_k(2,:)'-O_k(1,:)')/norm(O_k(2,:)'-O_k(1,:)');

%Rho_0_new=Rho_m_star+cross([0;0;1],[n_t(1);n_t(2);0])*R_gnc;
  
  
  
  
  
  
         Rho_0_new_t=Rho_m_star+(n_t)*R_gnc;
         [flag_coll]=Evaluate_collision(Rho_0_new_t',Surf,step_grid);
         
         
         if flag_coll==1
             
             disp('ooooops')
             
           Rho_0_new=Rho_m_star-(n_t)*R_gnc;
         else
             Rho_0_new=Rho_0_new_t;
           
         end
         
         
return