function [Rho_0_new]=Heuristic_tangent_smart(Rho_m_star,K,O_k,R_gnc,X_goal)


%[K,~]=size(O_k);
d_min=norm(Rho_m_star-O_k(1)');

n_t=(O_k(2)'-O_k(1)')/norm(O_k(2)'-O_k(1)');


if K >=3

 for k=2:K-1
     o_k=O_k(k,:)';
     d=norm(Rho_m_star-o_k);
     
     if d <= d_min
         
         k_min=k;
         o_k_min=O_k(k_min,:)';
         
         if (norm(X_goal-O_k(k_min-1,:)')+norm(Rho_m_star-X_goal)) >= (norm(X_goal-O_k(k_min+1,:)')+ norm(Rho_m_star-X_goal))
          
          
          n_t=(o_k_min-O_k(k_min-1,:)')/norm(o_k_min-O_k(k_min-1,:)');
         %disp('Pseudo-right')
         
         else
         
         
         %disp('Pseudo-left')
          n_t=(o_k_min-O_k(k_min+1,:)')/norm(o_k_min-O_k(k_min+1,:)'); 
         
         end
         
         
         
         
         

     end
     
  
     
 end
 
Rho_0_new=Rho_m_star+(n_t)*R_gnc;  

else

Rho_0_new=Rho_m_star;
disp('not enough surf point')
end

return