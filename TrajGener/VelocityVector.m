%***************************************
%@Author: Wentao Zhang
%@Date: 2023-3-10
%@E-mail: zwt190315@163.com
%@Reference: ##########
%@Problems: 
%@Description:
%@TODO：##########
%***************************************
clear
close all

YAW_ANG = -20;
ORIEN_VEL = 3;
VERDIT_VEL = 1.2;
UNIT_A = ORIEN_VEL * 0.8;
UNIT_B = VERDIT_VEL * 0.8;
%测试范围的速度
v = [1;1.8];

yaw = deg2rad(YAW_ANG);
% theta = -pi/2:0.1:pi/2;
theta = -pi:0.1:pi;
orien_vel = VERDIT_VEL + cos(theta)*(ORIEN_VEL- VERDIT_VEL);
vel_angle = yaw+theta;
% x_vel = orien_vel.*cos(vel_angle);
% y_vel = orien_vel.*sin(vel_angle);
x_vel = ORIEN_VEL*cos(theta);
y_vel = VERDIT_VEL*sin(theta);

vel(:,1) = x_vel;
vel(:,2) = y_vel;
Rmatrix = [cos(yaw),sin(yaw);-sin(yaw),cos(yaw)];
oval_vel = vel * Rmatrix;

figure
% subplot(1,2,1);
hold on
% plot(x_vel,y_vel,'k-');
% plot(oval_vel(:,1),oval_vel(:,2),'c-');
% plot(vel(:,1,:),vel(:,2),'b-');
vel_length = length(vel_angle);
for idx = 1:vel_length
%     plot([0,VERDIT_VEL*cos(vel_angle(idx))],[0,ORIEN_VEL*sin(vel_angle(idx))],'c-');
%     plot(VERDIT_VEL*cos(vel_angle(idx)),ORIEN_VEL*sin(vel_angle(idx)),'bo');
%     plot(VERDIT_VEL*cos(theta),ORIEN_VEL*sin(theta),'c-');
%     plot(VERDIT_VEL*cos(vel_angle),ORIEN_VEL*sin(vel_angle),'c-');
    plot(oval_vel(idx,1),oval_vel(idx,2),'bo','MarkerFaceColor','b');
    plot([0,oval_vel(idx,1)],[0,oval_vel(idx,2)],'c-');
end
plot(oval_vel(:,1),oval_vel(:,2),'b-');
% plot([0,UNIT_A*cos(yaw)],[0,UNIT_A*sin(yaw)],'r-','LineWidth',2);
% plot([0,UNIT_B*cos(yaw+pi/2)],[0,UNIT_B*sin(yaw+pi/2)],'g-','LineWidth',2);
% title('oval velocity vector');
% xlabel("x velocity");
% ylabel("y velocity");
% grid on
axis equal
axis off
% axis([-5,5,-5,5]);
% axis square 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 速度范围约束测试

v_l2 = sqrt(v(1)^2+v(2)^2);
% plot([0,v(1)],[0,v(2)],'r-','LineWidth',1);
% plot(v(1),v(2),'r.');
% Rmatrix = transpose(Rmatrix);
v_b = Rmatrix*v;
v_bl2 = v_b(1)^2/ORIEN_VEL^2+v_b(2)^2/VERDIT_VEL^2;
v_theta = atan(v(2)/v(1));
if (v_bl2 < 1)
    fprintf('Yes: v is in the oval \n');
else
    fprintf('No: v is‘t in the oval \n');
end


hold off
% figure
% % subplot(1,2,2);
% hold on
% plot(vel_angle,oval_vel(:,1),'r');
% plot(vel_angle,oval_vel(:,2),'b');
% legend('x vel','y vel');
% grid on
% subtitle("linear velocity");

%% 
% plot(vel(1,:),vel(2,:),'b-');
% print -depsc2 ellipsevel.eps