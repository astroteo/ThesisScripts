
Phi0_3d_plot_plus=1*zeros(M,N);
Phi0_3d_plot_minus=1*zeros(M,N);
Phi0_3d_plot_plus_Iz=1*zeros(M,N);



for i=2:N-1
     for j=2:M-1
         for k=2:P-1;
             
        x=R1_boundedx(i);
        y=R1_boundedy(j);
        z=R1_boundedz(k);
        
        
        
           if z <= 0 && z >= - height_ISS/2 || z > 0 && z<=height_ISS/2
            
            if y>=0 
                
                  if abs(x)<=lenght_ISS && abs(y) <= height_ISS/2
                      
                  Phi0_3d_plot_plus(i,j)=height_ISS/2;
                  Phi0_3d_plot_minus(i,j)=-height_ISS/2;
                
                  Phi_3d(i,j,k)=1;
                
                  elseif abs(x) <=lenght_ISS   && abs(x) >=2*height_ISS && abs(y)<=2*width_ISS
                      
                 Phi0_3d_plot_plus(i,j)=height_ISS/2;
                 Phi0_3d_plot_minus(i,j)=-height_ISS/2;
                     
                  Phi_3d(i,j,k)=1;
                      
                  end
                
            elseif y<0 
                
                
                if abs(x)<=height_ISS/2 && abs(y) <= 1*width_ISS
                    
                Phi0_3d_plot_plus(i,j)=height_ISS/2;
                Phi0_3d_plot_minus(i,j)=-height_ISS/2;
                
                Phi_3d(i,j,k)=1;
                
                
                end
                
               
            end
            
            
            if i==i_obs && j==j_obs && k==k_obs
            
            Phi_3d(i,j,k)=0;
            end
            
           elseif z > height_ISS/2 && z <=2*width_ISS
               
               
               if abs(x)< height_ISS/2 && abs(y) < height_ISS/2
               
                Phi0_3d_plot_plus_Iz(i,j)=z;
                Phi_3d(i,j,k)=1;
               end
               
           
               
        
           end
           
           
           
  Phi_3d_new(i,j,k)=1/6*(Phi_3d(i-1,j,k)+Phi_3d(i+1,j,k)+Phi_3d(i,j-1,k)+Phi_3d(i,j+1,k)+Phi_3d(i,j,k-1)+Phi_3d(i,j,k+1));
  
  
     
     
         end    
     end
 end
 
 
count=count+1;

if count-displaY_count >= 100
    
    displaY_count=displaY_count+100
    epsilon
end

figure()
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus')
colormap summer
hold on
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_minus')
colormap summer
hold on
contour3(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus_Iz',30)
colormap autumn
hold on
plot3(Rho_an(:,1),Rho_an(:,2),Rho_an(:,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(end,1),Rho_an(end,2),Rho_an(end,3),'or',x_obs,y_obs,z_obs,'*b')
title('collision Avoidance with discrete gradient & analytical trajectory propagation')
hold on
axis equal
grid on 


% figure()
% hold on
% mesh(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus')
% colormap summer
% mesh(R1_boundedx,R1_boundedy,Phi0_3d_plot_minus')
% colormap summer
% mesh(R1_boundedx,R1_boundedy,Phi0_3d_plot_plus_Iz')
% colormap summer
% %plot3(Rho_an(:,1),Rho_an(:,2),Rho_an(:,3),'--b',Rho_an(1,1),Rho_an(1,2),Rho_an(1,3),'og',Rho_an(end,1),Rho_an(end,2),Rho_an(end,3),'or',x_obs,y_obs,z_obs,'*b')
% title('modified obstacle')
% axis equal
% grid on 




