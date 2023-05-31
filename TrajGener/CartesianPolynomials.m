%***************************************
%@Author: Wentao Zhang
%@Date: 2023-3-5
%@E-mail: zwt190315@163.com
%@Description: unicycle trajectory generation by using cartesian polynomials
%@Reference: Siciliano, Bruno, Lorenzo Sciavicco, Luigi Villani和Giuseppe Oriolo. 
% Robotics Modelling, Planning and Control. Advanced Textbooks in Control and Signal Processing. 
% London: Springer London, 2009. https://doi.org/10.1007/978-1-84628-642-1.
% Ref[1]: 11.5.2 Flat Outputs
% if x_ds = y_ds = 0 than vs = 0, theta_s and ws not define the orientation
% but can be derived by continuity;
% 说人话就是，特殊点处，不能用定义算theta和ω了，但是可以根据连续性算出来一个
% Ref[2]: 11.5.3 Path Planning: Planning via Cartesian polynomials
%@Problems: 1) how to add constraints in start and goal point, like vt、ωt;
% 2)ki and kf must > 0
%@Tricks: x(s),y(s)定义的多项式好像比较随意呀，可以自己设计其它形式
%***************************************
%% Step1: Initial configuration
clc;clear;close all;
global k b ax ay bx by xi xf yi yf;
%cubic polynomial
% xs = s^3*xf - (s - 1)^3*xi + ax*s^2*(s - 1) + bx*s*(s - 1)^2;
% ys = s^3*yf - (s - 1)^3*yi + ay*s^2*(s - 1) + by*s*(s - 1)^2;
% x_ds = ax*s*(-2 + 3*s) + bx*(1 - 4*s + 3*s^2) + 3*s^2*xf - 3*xi + 6*s*xi - 3*s^2*xi;
% y_ds = ay*s*(-2 + 3*s) + by*(1 - 4*s + 3*s^2) + 3*s^2*yf - 3*yi + 6*s*yi - 3*s^2*yi;
% x_dds = 2*ax*(-1 + s) + 4*bx*(-1 + s) + 4*ax*s + 2*bx*s + 6*s*xf - 6*(-1 + s)*xi;
% y_dds = 2*ay*(-1 + s) + 4*by*(-1 + s) + 4*ay*s + 2*by*s + 6*s*yf - 6*(-1 + s)*yi;
% st = k*t + b;
T = 5;dt=0.02;
k = 1/T;b = 0;
% s_dt = k;
ki=20;kf=20;
%给定起始点

% qi=[8,7,0.3*pi];
% qf=[1,2,0.2*pi];

qi=[1,7,0.3*pi];
qf=[8,2,0.2*pi];

xi=qi(1);yi=qi(2);theta_i=qi(3);
xf=qf(1);yf=qf(2);theta_f=qf(3);
%添加起始点的约束
x_d_start = ki*cos(theta_i);
x_d_end   = kf*cos(theta_f);
y_d_start = ki*sin(theta_i);
y_d_end   = kf*sin(theta_f);
%求解多项式的系数解
ax=x_d_end - 3*xf   ;ay=y_d_end - 3*yf;
bx=x_d_start + 3*xi ;by=y_d_start + 3*yi;
% alpha =[x_d_end - 3*xf   ,y_d_end - 3*yf];
% beta  =[x_d_start + 3*xi ,y_d_start + 3*yi];
%xy的轨迹曲线
%% Step2: calculate function
t = 0:dt:T;
length_t = length(t);
st = k*t+b;
% st = k^2*t.^2+b;
xs = zeros(length_t,1);
ys = zeros(length_t,1);
vt = zeros(length_t,1);
wt = zeros(length_t,1);
thetas = zeros(length_t,1);
thetat = zeros(length_t,1);
for i=1:length_t
    xs(i) = fx_s(st(i));
    ys(i) = fy_s(st(i));
    vt(i) = fv_t(t(i));
    wt(i) = fw_t(t(i));
    thetas(i) = ftheta_s(st(i));
    thetat(i) = ftheta_t(t(i));
end
%% Step3: plot figure
% setup unicycle shape
unicycle_size=0.4;
vertices_unicycle_shape=unicycle_size*[[-0.25;-0.5;1/unicycle_size],...
    [0.7;0;1/unicycle_size],[-0.25;0.5;1/unicycle_size]];
faces_unicycle_shape=[1 2 3];
figure(1)
title("x-y 2D trajectory")
hold on
% plot unicycle initial configuration

M=[cos(theta_i) -sin(theta_i) xi; sin(theta_i) cos(theta_i)  yi;0 0 1]; 
vertices_unicycle_shape_0=(M*vertices_unicycle_shape)';
vertices_unicycle_shape_0=vertices_unicycle_shape_0(:,1:2);
patch('Faces',faces_unicycle_shape,'Vertices',vertices_unicycle_shape_0,...
    'FaceColor','r','EdgeColor','k');


% plot unicycle final configuration


M=[cos(theta_f) -sin(theta_f) xf; sin(theta_f) cos(theta_f)  yf;0 0 1];
vertices_unicycle_shape_f=(M*vertices_unicycle_shape)';
vertices_unicycle_shape_f=vertices_unicycle_shape_f(:,1:2);
patch('Faces',faces_unicycle_shape,'Vertices',vertices_unicycle_shape_f,...
    'FaceColor','b','EdgeColor','k'); %,'EraseMode','none'

text(xi+0.2,yi+0.2,"S");
text(xf+0.2,yf+0.2,"G");
% text(xi,yi,"S");
% text(xf,yf,"G");

plot(xs,ys);
% scatter(xi,yi,'red');
% scatter(xf,yf,'blue');
xlabel('x distance');
ylabel('y distance');
grid on
hold off
axis square
%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% % title("time axis for v and ω");
% subplot(2,1,1)
% plot(t,vt);
% title("time axis for v");
% xlabel('time / t');
% ylabel('velocity / vt');
% grid on
% subplot(2,1,2)
% plot(t,wt);
% subtitle("time axis for ω");
% xlabel('time / t');
% ylabel('velocity / ωt');
% grid on
%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(2,1,1)
title("theta for s")
plot(st,thetas);
title("theta for s")
xlabel('length / s');
ylabel('orientatition / theta');
grid on
subplot(2,1,2)
plot(t,thetat);
subtitle("theta for t")
xlabel('time / t');
ylabel('orientatition / theta');
grid on
%%%%%%%%%%%%%%%%%%%%%%
% figure(4)
% plot(t,st);
% title("s for time");
% xlabel('time / t');
% ylabel('distance / s');
%% PS：数值计算函数
%x_s
% function xs = fx_s(s,xf,xi,ax,bx)
function xs = fx_s(s)
global ax bx xi xf;
    xs = s^3*xf - (s - 1)^3*xi + ax*s^2*(s - 1) + bx*s*(s - 1)^2;
end
%y_s
% function ys = fy_s(s,yf,yi,ay,by)
function ys = fy_s(s)
global  ay  by yi yf;
    ys = s^3*yf - (s - 1)^3*yi + ay*s^2*(s - 1) + by*s*(s - 1)^2;
end
%x_ds
% function x_ds = fx_ds(s,ax,bx,xf,xi)
function x_ds = fx_ds(s)
global ax bx xi xf;
    x_ds = ax*s*(-2 + 3*s) + bx*(1 - 4*s + 3*s^2) + 3*s^2*xf - 3*xi + 6*s*xi - 3*s^2*xi;
end %end x_ds
%y_ds
% function y_ds = fy_ds(s,ay,by,yf,yi)
function y_ds = fy_ds(s)
global  ay by yi yf;
    y_ds = ay*s*(-2 + 3*s) + by*(1 - 4*s + 3*s^2) + 3*s^2*yf - 3*yi + 6*s*yi - 3*s^2*yi;
end %end x_ds
%x_dds
% function x_dds = fx_dds(s,ax,bx,xf,xi)
function x_dds = fx_dds(s)
global ax  bx xi xf;
    x_dds = 2*ax*(-1 + s) + 4*bx*(-1 + s) + 4*ax*s + 2*bx*s + 6*s*xf - 6*(-1 + s)*xi;
end %end x_dds
%y_dds
% function y_dds = fy_dds(s,ay,by,yf,yi)
function y_dds = fy_dds(s)
global ay by yi yf;
    y_dds = 2*ay*(-1 + s) + 4*by*(-1 + s) + 4*ay*s + 2*by*s + 6*s*yf - 6*(-1 + s)*yi;
end %end x_ds
% vs
function vs=fv_s(s)
    vs = sqrt(fx_ds(s)^2+fy_ds(s)^2);
end %end vs
% ws
function ws=fw_s(s)
    ws = (fy_dds(s)*fx_ds(s)-fx_dds(s)*fy_ds(s))/(fx_ds(s)^2+fy_ds(s)^2);
end %end vs
% fs_t
function st=fs_t(t)
global k b;
    st = k*t+b;
end %end st
% fv_t
function vt=fv_t(t)
% global s_dt;
    s_dt = fs_dt(t);
    st = fs_t(t);
    vt = fv_s(st)*s_dt;
end %end fv_t
% fw_t
function wt=fw_t(t)
% global s_dt;
    s_dt = fs_dt(t);
    st = fs_t(t);
    wt = fw_s(st)*s_dt;
end %end fw_t
% ftheta_t
function theta_t=ftheta_t(t)
    % + k*π; k=0,forward motion;k=1,backward motion;
    st = fs_t(t);
    theta_t =atan2(fy_ds(st),fx_ds(st)); 
end %end fw_t
% ftheta_s
function theta_s=ftheta_s(s)
    % + k*π; k=0,forward motion;k=1,backward motion;
    theta_s =atan2(fy_ds(s),fx_ds(s)); 
end %end ftheta_s
% fs_dt
function s_dt=fs_dt(t)
global k;
%     s_dt = 2*k^2*t;
s_dt = k;
end %end st
