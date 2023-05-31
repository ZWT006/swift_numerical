%***************************************
%@Author: Wentao Zhang
%@Date: 2023-3-23
%@E-mail: zwt190315@163.com
%@Reference: ##########
%@Problems: 经过代码测试,发现使用
%@Description: Bezier Curve calculate
%@TODO：
%***************************************
close all
clear
clc
n=0;
dt=0.01;
T = 1;
ts=3;
%##########################################################################
%计算Bezier曲线
control_points = (ginput()-0.5)*100;
% load('control_points.mat');
path_length = length(control_points);
x_p = control_points(:,1);
y_p = control_points(:,2);
order = path_length-1;
M_c = BerneteinCoeff(order);
% M_c = getM(path_length-1);
x_poly=M_c*x_p;
y_poly=M_c*y_p;

x_poly=flip(x_poly);
y_poly=flip(y_poly);

t=0:dt:T;
n=n+1;
figure(n)
hold on
plot(control_points(:,1),control_points(:,2),'ob')
x_pos = polyval(x_poly,t);
y_pos = polyval(y_poly,t);
plot(x_pos,y_pos,'k-');
title('xy bezier pose')
axis equal
grid on
hold off
%##########################################################################
%计算Bezier xy的变化
n=n+1;
figure(n)
hold on
subplot(1,2,1)
plot(t,x_pos,'r-')
subtitle('x pose')

grid on
subplot(1,2,2)
plot(t,y_pos,'b-')
subtitle('y pose')

grid on
hold off

%##########################################################################
%计算Bezier的一阶导: c'i=n(ci - c(i+1))
n=n+1;
figure(n)
hold on
control_points_d(:,1) = order*(x_p(2:end)-x_p(1:end-1));
control_points_d(:,2) = order*(y_p(2:end)-y_p(1:end-1));
plot(control_points_d(:,1),control_points_d(:,2),'ob')
xd_pos = polyval(polyder(x_poly),t);
yd_pos = polyval(polyder(y_poly),t);
plot(xd_pos,yd_pos,'k--','LineWidth',2);

xd_p = control_points_d(:,1);
yd_p = control_points_d(:,2);
M_c = BerneteinCoeff(order-1);
xd_poly=M_c*xd_p;
yd_poly=M_c*yd_p;
xd_poly=flip(xd_poly);
yd_poly=flip(yd_poly);
xdc_pos = polyval(xd_poly,t);
ydc_pos = polyval(yd_poly,t);
plot(xdc_pos,ydc_pos,'r--');
grid on
title('bezier first-order derivative')

%##########################################################################
%计算具有时间缩放的Bezier
M_c = BerneteinCoeff(order);
coeff = getCoeffCons(1/ts, order, 1);
coeff=coeff';
% sx_ploy=M_c*x_p.*coeff;
% sy_ploy=M_c*y_p.*coeff;
coeff=flip(coeff);
sx_ploy=coeff.*x_poly;
sy_ploy=coeff.*y_poly;
n=n+1;
figure(n)
t=0:dt:ts;
hold on
plot(control_points(:,1),control_points(:,2),'ob')
sx_pos = polyval(sx_ploy,t);
sy_pos = polyval(sy_ploy,t);
plot(sx_pos,sy_pos,'k-');
title('bezier curve time = ts')
axis equal
grid on
hold off

%##########################################################################
%对比时间缩放的速度问题
n=n+1;
figure(n)
subplot(1,2,1)
hold on
plot(control_points_d(:,1)/ts,control_points_d(:,2)/ts,'*r')
sxd_ploy=polyder(sx_ploy);
syd_ploy=polyder(sy_ploy);
sxd_pos = polyval(sxd_ploy,t);
syd_pos = polyval(syd_ploy,t);
plot(sxd_pos,syd_pos,'k-');
grid on
title('t = ts')

subplot(1,2,2)
hold on
plot(control_points_d(:,1),control_points_d(:,2),'ob')
plot(xdc_pos,ydc_pos,'r--','LineWidth',2);
grid on
title('t = 1')

%##########################################################################
%对比时间缩放的加速度问题

%计算Bezier曲线的控制点
control_points_dd(:,1) = (order-1)*(xd_p(2:end)-xd_p(1:end-1));
control_points_dd(:,2) = (order-1)*(yd_p(2:end)-yd_p(1:end-1));
xdd_p = control_points_dd(:,1);
ydd_p = control_points_dd(:,2);

M_c = BerneteinCoeff(order-2);
xdd_poly=M_c*xdd_p;
ydd_poly=M_c*ydd_p;
xdd_poly=flip(xdd_poly);
ydd_poly=flip(ydd_poly);
t=0:dt:T;
xddc_pos = polyval(xdd_poly,t);
yddc_pos = polyval(ydd_poly,t);

% 

t=0:dt:ts;
n=n+1; 
figure(n)
subplot(1,3,1)
hold on
plot(control_points_dd(:,1)/ts^2,control_points_dd(:,2)/ts^2,'*r')
sxdd_ploy=polyder(sxd_ploy);
sydd_ploy=polyder(syd_ploy);
sxdd_pos = polyval(sxdd_ploy,t);
sydd_pos = polyval(sydd_ploy,t);
plot(sxdd_pos,sydd_pos,'k-');
grid on
title('t = ts')

subplot(1,3,2)
hold on
plot(xdd_p,ydd_p,'ob')
plot(xddc_pos,yddc_pos,'r-' );
grid on
title('t = 1')

xdd_p = control_points_dd(:,1)/ts^2;
ydd_p = control_points_dd(:,2)/ts^2;

M_c = BerneteinCoeff(order-2);
xdd_poly=M_c*xdd_p;
ydd_poly=M_c*ydd_p;
xdd_poly=flip(xdd_poly);
ydd_poly=flip(ydd_poly);
t=0:dt:T;
xddc_pos = polyval(xdd_poly,t);
yddc_pos = polyval(ydd_poly,t);

subplot(1,3,3)
hold on
plot(xdd_p,ydd_p,'ob')
plot(xddc_pos,yddc_pos,'b-');
grid on
title('t = 1')