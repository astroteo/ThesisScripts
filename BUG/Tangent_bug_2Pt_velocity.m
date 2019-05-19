function [Rho_m_star,Ni_m_star,O_k,K]=Tangent_bug_2Pt_velocity(Rho_m,Ni_m,Surf,step_grid,R_gnc)


flag_star=0;
i=1;
k=0;

[N_surf,M_surf]=size(Surf);

O_k=zeros(N_surf,3);

[N_m,~]=size(Rho_m);

Rho_m_star=Rho_m(1,:)';

while flag_star==0 && i <= N_m;
    
    Rho=Rho_m(i,:)';
    Ni=Ni_m(i,:)';
    
    
    
    for i_surf=1:N_surf %-> this cycle is meant to be improved...lokking for closer points first.
        for j_surf=1:M_surf
            
            if Surf(i_surf,j_surf)==1
            x=(i_surf-N_surf/2)*step_grid;
            y=(j_surf-M_surf/2)*step_grid;
            r_surf=[x;y;0];
            
              if abs(norm(Rho-r_surf)-R_gnc) <= step_grid/2 
                  
                  k=k+1;
                  O_k(k,:)=r_surf';%-> containing interception points
                  
                  Rho_m_star=Rho;
                  Ni_m_star=Ni;
                  flag_star=1;
                  
              end
            end
            
        end
    end
    
    
i=i+1;  
   
end

K=k;
return