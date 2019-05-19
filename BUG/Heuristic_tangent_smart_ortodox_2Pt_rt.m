function[Rho_0_new,n_t,flag_control]=Heuristic_tangent_smart_ortodox_2Pt_rt(Rho_m_star,Rho,Ni,K,O_k,R_gnc,X_goal,Surf,step_grid)

if K<=1
    flag_control=0;
    Rho_0_new=Rho;
    n_t=zeros(3,1);
    
elseif K<3 && K>1

n_t=(O_k(2,:)'-O_k(1,:)')/norm(O_k(2,:)'-O_k(1,:)');

%Rho_0_new=Rho_m_star+cross([0;0;1],[n_t(1);n_t(2);0])*R_gnc;
Rho_0_new=Rho_m_star+n_t*R_gnc;

flag_control=1;


  
else
    
flag_control=1;
  
  o_k_L=O_k(1,:);
  
  o_k_R=O_k(K,:);
  
  
  if (norm(X_goal-o_k_L')+norm(Rho-X_goal)) >= (norm(X_goal-o_k_R')+ norm(Rho-X_goal))
   %if (norm(X_goal-o_k_L')) >= (norm(X_goal-o_k_R'))
     %o_k_min=o_k_L;
     n_t_T=(o_k_L'-O_k(2,:)')/norm(o_k_L'-O_k(2,:)');
     
     Rho_0_new_R=Rho_m_star+(n_t_T)*R_gnc;
     Rho_0_new_L=Rho_m_star-(n_t_T)*R_gnc;
     
     
        if dot((Rho_0_new_R-Rho),Ni) <= dot((Rho_0_new_L-Rho),Ni)
     
            Rho_0_new_t=Rho_0_new_L;
            n_t=-n_t_T;
         
        else
         
            Rho_0_new_t=Rho_0_new_R;
            n_t=n_t_T;
     
       
        end
     
     
  else
      
     %o_k_min=o_k_R;
     n_t_T=(o_k_R'-O_k(K-1,:)')/norm(O_k(K-1,:)'-o_k_R');
     
     Rho_0_new_R=Rho_m_star+(n_t_T)*R_gnc;
     Rho_0_new_L=Rho_m_star-(n_t_T)*R_gnc;
     
       if dot((Rho_0_new_R-Rho),Ni) <= dot((Rho_0_new_L-Rho),Ni)
     
            Rho_0_new_t=Rho_0_new_L;
            n_t=-n_t_T;
         
        else
         
            Rho_0_new_t=Rho_0_new_R;
            n_t=n_t_T;
     
       
        end
     
     
  end
  
         
         [flag_coll]=Evaluate_collision(Rho_0_new_t',Surf,step_grid);
         
         
         if flag_coll==1
             
             disp('ooooops')
             
           Rho_0_new=Rho_m_star-(n_t)*R_gnc;
         else
             Rho_0_new=Rho_0_new_t;
           
         end
         
end        
return




