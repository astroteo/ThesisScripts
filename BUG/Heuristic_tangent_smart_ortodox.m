function [Rho_0_new,n_t]=Heuristic_tangent_smart_ortodox(Rho_m_star,Rho,K,O_k,R_gnc,X_goal,Surf,step_grid)


%[K,~]=size(O_k);
d_min=norm(Rho_m_star-O_k(1)')+norm(X_goal-O_k(1));

n_t=(O_k(2,:)'-O_k(1,:)')/norm(O_k(2,:)'-O_k(1,:)');

%Rho_0_new=Rho_m_star+cross([0;0;1],[n_t(1);n_t(2);0])*R_gnc;
  Rho_0_new=Rho_m_star+n_t*R_gnc;
  
  if K< 3
      %disp('aksjdbakjdb')
  end

 for k=2:K-1
     o_k=O_k(k,:)';
   if norm(o_k-X_goal) > norm( Rho-X_goal)
         
     d=norm(Rho_m_star-o_k)+norm(X_goal-o_k);
     
     if d <= d_min 
         
         k_min=k;
         o_k_min=O_k(k_min,:)';
         
         if (norm(X_goal-O_k(k_min-1,:)')+norm(Rho-X_goal)) >= (norm(X_goal-O_k(k_min+1,:)')+ norm(Rho-X_goal))
%          if norm(X_goal-O_k(k_min-1,:)') >= norm(X_goal-O_k(k_min+1,:)')
             
           n_t=(o_k_min-O_k(k_min-1,:)')/norm(o_k_min-O_k(k_min-1,:)');
          
          
         %disp('Pseudo-right')
         
         else
         
         n_t=(o_k_min-O_k(k_min-1,:)')/norm(o_k_min-O_k(k_min-1,:)');
         %disp('Pseudo-left')
        
         
         end
         
         Rho_0_new_t=Rho_m_star+(n_t)*R_gnc;
         
         [flag_coll]=Evaluate_collision(Rho_0_new_t',Surf,step_grid);
         
         
         if flag_coll==1
             
           Rho_0_new=Rho_m_star-(n_t)*R_gnc;
         else
             Rho_0_new=Rho_0_new_t;
           
         end
         
         


     end
     
  end
     
 end
 
 


return