function [Rho_0_new]=Heuristic(Rho_m_star,K,O_k,R_gnc)

%[K,~]=size(O_k);
d_min=norm(Rho_m_star-O_k(1)');

n_v=(Rho_m_star-O_k(1)')/norm(Rho_m_star-O_k(1)');

 for k=2:K
     o_k=O_k(k,:)';
     d=norm(Rho_m_star-o_k);
     
     if d <=d_min
         
         k_min=k;
         o_k_min=O_k(k_min,:)';
         n_v=(Rho_m_star-o_k_min)/norm(Rho_m_star-o_k_min);

         
         
  
     end
     
 end
 
 

Rho_0_new=Rho_m_star+(n_v)*R_gnc;

return