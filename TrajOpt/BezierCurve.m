%***************************************
%@Author: Wentao Zhang
%@Date: 2023-3-23
%@E-mail: zwt190315@163.com
%@Reference: ##########
%@Problems: 
%@Description: Bezier Curve calculate
%@TODO：
%***************************************
close all
clear
clc

dt=0.01;
T = 1;

control_points = (ginput()-0.5)*100;
path_length = length(control_points);
x_p = control_points(:,1);
y_p = control_points(:,2);
M_c = BerneteinCoeff(path_length-1);
% M_c = getM(path_length-1);
x_poly=M_c*x_p;
y_poly=M_c*y_p;

x_poly=flip(x_poly);
y_poly=flip(y_poly);

t=0:dt:T;
hold on
plot(control_points(:,1),control_points(:,2),'ob')
x_pos = polyval(x_poly,t);
y_pos = polyval(y_poly,t);
plot(x_pos,y_pos,'k-');
axis equal
grid on
hold off