clear all
close all

L_max=10;
step_gridx=0.1;
R1_boundedx=-L_max:step_gridx:(L_max-step_gridx);

step_gridy=0.1;
R1_boundedy=-L_max:step_gridy:(L_max-step_gridy);



[~,N]=size(R1_boundedx);
[~,M]=size(R1_boundedy);

phi_harm=zeros(N,M);
phi_scod=zeros(N,M);

grad_scod=zeros(N,M,2);
grad_scod_dir=zeros(N,M,2);

grad_harm=zeros(N,M,2);
grad_harm_dir=zeros(N,M,2);

%potential & gradient discrete.
for i=1:N
    
    for j=1:M
        
       
            
            x=R1_boundedx(1,i)/L_max;
            y=R1_boundedy(1,j)/L_max;
           
         
            
           phi_scod(i,j)=x^2+y^2;
           phi_harm(i,j)=exp(x^2+y^2);
           
           grad_scod(i,j,1)=2*x+y^2;
           grad_scod(i,j,2)=2*y+x^2;
           
           grad_harm(i,j,1)=2*x*exp(x^2+y^2);
           grad_harm(i,j,2)=2*y*exp(x^2+y^2);
           
           
           grad_scod_dir(i,j,1)=grad_scod(i,j,1)/norm([grad_scod(i,j,1);grad_scod(i,j,2)]);
           grad_scod_dir(i,j,2)=grad_scod(i,j,2)/norm([grad_scod(i,j,1);grad_scod(i,j,2)]);
           
           grad_harm_dir(i,j,1)=grad_harm(i,j,1)/norm([grad_harm(i,j,1);grad_harm(i,j,2)]);
           grad_harm_dir(i,j,2)=grad_harm(i,j,2)/norm([grad_harm(i,j,1);grad_harm(i,j,2)]);
           
     end
end

x_t=0
y_t=0

total_stepx=N;
total_stepy=M;


if x_t >=0
i_bottom=total_stepx/2+floor(x_t/step_gridx);
i_up=i_bottom+1;

xd=(x_t-(i_bottom-total_stepx/2)*step_gridx);

    
else     
i_up=ceil((total_stepx/2*step_grid-abs(x_t))/step_gridx);
i_bottom=i_up-1;

xd=x_t+(total_stepx/2-i_up)*step_gridx;

end


if y_t >=0
j_bottom=total_stepy/2+floor(y_t/step_gridy);
j_up=j_bottom+1;

yd=(y_t-(j_bottom-total_stepy/2)*step_gridy);

    
else     
j_up=ceil((total_stepy/2*step_gridy-abs(y_t))/step_gridy);
j_bottom=j_up-1;

yd=y_t+(total_stepy/2-j_up)*step_gridy;
end

scod_or_harm=1;% 1 to evaluate harmonic potential 0 to evaluate scod. potential

if scod_or_harm==1
    disp('evaluating harmonic potential')
   phi_2d=phi_harm;
   grad_AN_x=2*x_t*exp(x_t^2+y_t^2);
   grad_AN_y=2*y_t*exp(x_t^2+y_t^2);
   grad_AN=[grad_AN_x;grad_AN_y];
   grad_dir_AN=[grad_AN_x/norm(grad_AN);grad_AN_y/norm(grad_AN)];
else
    disp('evaluating scodella potential')
    phi_2d=phi_scod;
    
   grad_AN_x= 2*x_t+y_t^2;
   grad_AN_y= 2*y_t+x_t^2;
   
   grad_AN=[grad_AN_x;grad_AN_y];
   grad_dir_AN=[grad_AN_x/norm(grad_AN);grad_AN_y/norm(grad_AN)];
end


fx=-phi_2d(i_bottom,j_bottom)*(1-yd)  + phi_2d(i_up,j_bottom)*(1-yd) - phi_2d(i_bottom,j_up)*yd     + phi_2d(i_up,j_up)*yd;
fy=-phi_2d(i_bottom,j_bottom)*(1-xd)  - phi_2d(i_up,j_bottom)*xd     + phi_2d(i_bottom,j_up)*(1-xd) + phi_2d(i_up,j_up)*xd;

grad_mod_MY=norm([fx,fy]);
grad_dir_MY=[fx/grad_mod_MY;fy/grad_mod_MY ];


[p_harmx,p_harmy] = gradient(phi_harm,step_gridx,step_gridy);
[p_scodx,p_scody] = gradient(phi_scod,step_gridx,step_gridy);

i_up
i_bottom

j_up
j_bottom

xd
yd


grad_dir_MY
grad_dir_AN

figure()
mesh(R1_boundedx,R1_boundedy,phi_harm)
title('harmonic potential')
grid on
hold on

figure()
mesh(R1_boundedx,R1_boundedy,phi_scod)
title('scod potential')
grid on
hold on

figure()
contour(R1_boundedx,R1_boundedy,phi_harm)
title('harmonic potential contour')
hold on
quiver(R1_boundedx,R1_boundedy,p_harmx,p_harmy)
hold off

figure()
contour(R1_boundedx,R1_boundedy,phi_scod)
title('scod potential contour')
hold on
quiver(R1_boundedx,R1_boundedy,p_scodx,p_scody)
hold off


