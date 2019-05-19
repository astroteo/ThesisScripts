[x_ISS, y_ISS, z_ISS] = ellipsoid(0,0,0,50,5,10);
figure()
%plot3(Rrel_LVLH_LQR(:,1),Rrel_LVLH_LQR(:,2),Rrel_LVLH_LQR(:,3),'r')
hold on
surf(x_ISS, y_ISS, z_ISS,'facecolor','g')
axis equal
hold on
grid on