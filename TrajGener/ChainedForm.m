%***************************************
%@Author: Wentao Zhang
%@Date: 2023-3-6
%@E-mail: zwt190315@163.com
%@Description: unicycle trajectory generation by using cartesian polynomials
%@Reference: Siciliano, Bruno, Lorenzo Sciavicco, Luigi Villani和Giuseppe Oriolo. 
% Robotics Modelling, Planning and Control. Advanced Textbooks in Control and Signal Processing. 
% London: Springer London, 2009. https://doi.org/10.1007/978-1-84628-642-1.
% 11.5.3 Path Planning: Planning via chained form
%@Problem: 链式法则搞出来的路径比较极端
%***************************************
%% Step1: Initial configuration
clc;clear;close all;
global k b z1f z1i z3f z3i a3 b3 xi xf yi yf;
%cubic polynomial
% 这样定义z1则角速度就就从开始到结束匀速运动
% z1s = z1f*s - (s - 1)*z1i;
% z3s = s^3*z3f - (s - 1)^3*z3i + a3*s^2*(s - 1) + b3*s*(s - 1)^2;
% z2s = z3_ds/z1_ds = (b3*(-1 + s)^2 + 2*a3*(-1 + s)*s + 2*b3*(-1 + s)*s + a3*s^2 + 
% 3*s^2*z3f - 3*(-1 + s)^2*z3i)/(z1f - z1i);
% z1_ds = z1f - z1i;
% z3_ds = a3*s*(-2 + 3*s) + b3*(1 - 4*s + 3*s^2) + 3*s^2*z3f - 3*z3i +
% 6*s*z3i - 3*s^2*z3i;
% z2_ds = (2*a3*(-1 + s) + 4*b3*(-1 + s) + 4*a3*s + 2*b3*s + 6*s*z3f - 
% 6*(-1 + s)*z3i)/(z1f - z1i)
% st = k*t + b;
T = 5;dt=0.02;
k = 1/T;b = 0;
% s_dt = k;
% ki=10;kf=10;
%给定起始点
qi=[1,2,0.2*pi];
qf=[8,7,-0.2*pi];
z1i = qi(3); z1f = qf(3);
z2i = fz2(qi(1),qi(2),qi(3));
z2f = fz2(qf(1),qf(2),qf(3));
z3i = fz3(qi(1),qi(2),qi(3));
z3f = fz3(qf(1),qf(2),qf(3));
%添加起始点的约束

%求解多项式的系数解
a3 = z2f*(z1f - z1i) - 3*z3f;
b3 = z2i*(z1f - z1i) + 3*z3i;
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
    xs(i) = fz2_s(st(i))*cos(ftheta_s(st(i))) + fz3_s(st(i))*sin(ftheta_s(st(i)));
    ys(i) = fz2_s(st(i))*sin(ftheta_s(st(i))) - fz3_s(st(i))*cos(ftheta_s(st(i)));
    vt(i) = fv_t(t(i));
    wt(i) = fw_t(t(i));
    thetas(i) = ftheta_s(st(i));
    thetat(i) = ftheta_t(t(i));
end
%% Step3: plot figure
xi=qi(1);yi=qi(2);theta_i=qi(3);
xf=qf(1);yf=qf(2);theta_f=qf(3);
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

plot(xs,ys);
% scatter(xi,yi,'red');
% scatter(xf,yf,'blue');
xlabel('x distance');
ylabel('y distance');
grid on
hold off
axis square
figure(2)
% title("time axis for v and ω");
subplot(2,1,1)
plot(t,vt);
title("time axis for v");
xlabel('time / t');
ylabel('velocity / vt');
grid on
subplot(2,1,2)
plot(t,wt);
subtitle("time axis for ω");
xlabel('time / t');
ylabel('velocity / ωt');
grid on
figure(3)
subplot(2,1,1)
% title("theta for s")
plot(st,thetas);
title("theta for s ")
xlabel('length / s');
ylabel('orientatition / theta');
grid on
subplot(2,1,2)
plot(t,thetat);
subtitle("theta for t")
xlabel('time / t');
ylabel('orientatition / theta');
grid on
figure(4)
plot(t,st);
title("s for time");
xlabel('time / t');
ylabel('distance / s');
%% PS：数值计算函数
%z1s
function z1s = fz1_s(s)
global z1f z1i;
   z1s = z1f*s - (s - 1)*z1i;
end % end z1(s)
%z3s
function z3s = fz3_s(s)
global z3f z3i a3 b3;
   z3s = s^3*z3f - (s - 1)^3*z3i + a3*s^2*(s - 1) + b3*s*(s - 1)^2;
end % end z3(s)
%z2s
function z2s = fz2_s(s)
global z1f z1i z3f z3i a3 b3;
   z2s = (b3*(-1 + s)^2 + 2*a3*(-1 + s)*s + 2*b3*(-1 + s)*s + a3*s^2 + ...
          3*s^2*z3f - 3*(-1 + s)^2*z3i)/(z1f - z1i);
end % end z2(s)
%fz1ds()
function z1ds = fz1_ds(~)
global z1f z1i;
    z1ds = z1f-z1i;
end % end d_z1(s)/d_s
%fz2ds(s)
function z2ds = fz2_ds(s)
global z1f z1i z3f z3i a3 b3;
    z2ds = (2*a3*(-1 + s) + 4*b3*(-1 + s) + 4*a3*s + 2*b3*s + 6*s*z3f - ...
6*(-1 + s)*z3i)/(z1f - z1i);
end % end d_z2(s)/d_s
%fz3ds(s)
function z3ds = fz3_ds(s)
global z3f z3i a3 b3;
     z3ds = a3*s*(-2 + 3*s) + b3*(1 - 4*s + 3*s^2) + 3*s^2*z3f - 3*z3i +...
6*s*z3i - 3*s^2*z3i;
end % end d_z3(s)/d_s
% z1d(t)
function z1dt = fz1_dt(t)
    st = fs_t(t);
    z1dt = fz1_ds(st)*fs_dt(t);
end
% z2d(t)
function z2dt = fz2_dt(t)
    st = fs_t(t);
    z2dt = fz2_ds(st)*fs_dt(t);
end
% z2d(t)
function z3dt = fz3_dt(t)
    st = fs_t(t);
    z3dt = fz3_ds(st)*fs_dt(t);
end
% vt
function vt = fv_t(t)
    vt = fz2_dt(t) + fz3_dt(t);
end
% wt
function wt = fw_t(t)
    wt = fz1_dt(t);
end
% z2 = x*cos(θ) + y*sin(θ)
function z2 = fz2(x,y,theta)
    z2 = x*cos(theta) + y*sin(theta);
end % end fz2
% z3 = x*sin(θ) + y*cos(θ)
function z3 = fz3(x,y,theta)
    z3 = x*sin(theta) + y*cos(theta);
end % end fz2
%ftheta_s
function thetas = ftheta_s(st)
    thetas = fz1_s(st);
end
%ftheta_t
function thetat = ftheta_t(t)
    st = fs_t(t);
    thetat = fz1_s(st);
end
% %fy_s
% function ys = fy_s(s)
%     ys = fz2_s(s)*sin(ftheta_s(s)) - fz3_s(s)*cos(ftheta_s(s));
% end
% %fx_s
% function xs = fx_s(s)
%     xs = fz2_s(s)*cos(ftheta_s(s)) - fz3_s(s)*sin(ftheta_s(s));
% end
% fs_t
function st=fs_t(t)
global k b;
    st = k*t+b;
end %end st
% fs_dt
function s_dt=fs_dt(~)
global k;
%     s_dt = 2*k^2*t;
s_dt = k;
end %end st
