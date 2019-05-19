function [Rho_0_new]=Real_Heuristic(Rho_m_star,K,O_k,X_goal,R_gnc)

%[K,~]=size(O_k);
d_max=norm(Rho_m_star-O_k(1)')+norm(X_goal-O_k(1)');

n_v=(Rho_m_star-O_k(1)')/norm(Rho_m_star-O_k(1)');

 for k=2:K
     o_k=O_k(k,:)';
     d=norm(Rho_m_star-o_k)+norm(X_goal-o_k);

     
     if d >=d_max
         
         k_max=k;
         o_k_max=O_k(k_max,:)';
         n_v=(Rho_m_star-o_k_max)/norm(Rho_m_star-o_k_max);

         
         
  
     end
     
 end
 
 

Rho_0_new=Rho_m_star+(n_v)*R_gnc;

return