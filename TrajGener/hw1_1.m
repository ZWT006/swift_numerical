clc;clear;close all;
% path = ginput() * 100.0;
path = [50,	                50;
        145.659906173533,	186.312819235911;
        191.836813342538,	272.271495258284;
        267.755189847367,	328.985619044511;
        392.144396900521,	385.147806029386;
        514.147782625023,	469.962638842209;
        613.408785998573,	537.502971052556;
        702.631862930347,	605.144253906707;
        794.630236687311,	668.449239620839;
        852.120251855336,	664.522777123541;
        940,	            550];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% segment polynomial optimization coefficients
% n > i > j
% polynomial order  : n
% input order   : i (i > j means input control minimum value)
% cost order    : j (cost function order)
% segmant num   : k
% segmant dof   : k*(n+1)
% constraint num: j(start state)+j(final state)+(k-1)(mid points pose)
% unknow num    : k*(n+1)
% equal constraint: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_order       = 7;          % 多项式的阶数
n_costorder   = 2;          % 最小化的求导阶次 0=posi;1=vel;2=acc;3=jerk;4=snap;
n_inputorder  = 3;          % 输入的阶次 可以理解为segment 之间满足等式约束的阶次
n_seg         = size(path,1)-1; % 分段数
n_poly_perseg = (n_order+1);    % 每段的系数个数

ts = zeros(n_seg, 1);

%根据距离长度计算每两个点之间的时间分配
dist     = zeros(n_seg, 1);
dist_sum = 0;
T        = 60;
t_sum    = 0;

for i = 1:n_seg
    dist(i) = sqrt((path(i+1, 1)-path(i, 1))^2 + (path(i+1, 2) - path(i, 2))^2);
    dist_sum = dist_sum+dist(i);
end
for i = 1:n_seg-1
    ts(i) = dist(i)/dist_sum*T;
    t_sum = t_sum+ts(i);
end
ts(n_seg) = T - t_sum;

% % 简单地设置每一段的时间分配为1
% for i = 1:n_seg
%     ts(i) = 1.0;
% end

% 分别对x和y方向求取对应的多项式系数
% poly_coef_x = MinimumSnapQPSolver(path(:, 1), ts, n_seg, n_order, n_costorder, n_inputorder);
% poly_coef_y = MinimumSnapQPSolver(path(:, 2), ts, n_seg, n_order, n_costorder, n_inputorder);
poly_coef_x = MinimumQPSolver(path(:, 1), ts, n_seg, n_order, n_costorder, n_inputorder);
poly_coef_y = MinimumQPSolver(path(:, 2), ts, n_seg, n_order, n_costorder, n_inputorder);


% 用于显示轨迹
X_n = [];
Y_n = [];
k = 1;
tstep = 0.01;
for i=0:n_seg-1
    %#####################################################
    % STEP 3: get the coefficients of i-th segment of both x-axis
    % and y-axis
    Pxi = poly_coef_x((n_order+1)*(i)+1:(n_order+1)*(i)+n_order+1); 
    Pyi = poly_coef_y((n_order+1)*(i)+1:(n_order+1)*(i)+n_order+1);
    
    for t = 0:tstep:ts(i+1)
        px=flip(Pxi);
        py=flip(Pyi);
        X_n(k)  = polyval(px, t);
        Y_n(k)  = polyval(py, t);
        pdx=polyder(px);
        pdy=polyder(py);
        X_dn(k)  = polyval(pdx, t);
        Y_dn(k)  = polyval(pdy, t);
        k = k + 1;
    end
end

figure (1)
plot(X_n, Y_n , 'g-');
hold on
scatter(path(1:size(path, 1), 1), path(1:size(path, 1), 2));
grid on
axis equal
xlim([0 1000]);
ylim([0 800]);
figure (2)
hold on
tv=0:tstep:(k-2)*tstep;
plot(tv,X_dn, 'r-');
plot(tv,Y_dn, 'b-');
legend('x vel','y vel');
grid on
xlim([0 T]);
% ylim([0 800]);
% Minisnap求解器
function poly_coef = MinimumSnapQPSolver(waypoints, ts, n_seg, n_order, n_costorder, n_inputorder)
    % 起点约束
    start_cond = [waypoints(1), 0, 0, 0];
    % 终点约束
    end_cond   = [waypoints(end), 0, 0, 0];
    %#####################################################
    % STEP 1: 计算Q矩阵
    Q = getQ(n_seg, n_order, n_costorder,ts);
    
    %#####################################################
    % STEP 2: 计算对应的约束矩阵A_beq
    [Aeq, beq] = getAbeq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond);
    
    f = zeros(size(Q,1),1);
    % 求解多项式系数
    poly_coef = quadprog(Q,f,[],[],Aeq, beq);
end

% Minisnap求解器
function poly_coef = MinimumQPSolver(waypoints, ts, n_seg, n_order,n_costorder,n_inputorder)
    % 起点约束
    start_cond = zeros(1,n_inputorder);
    start_cond(1) = waypoints(1);
    % 终点约束
    end_cond = zeros(1,n_inputorder);
    end_cond(1) = waypoints(end);
    %#####################################################
    % STEP 1: 计算Q矩阵
    Q = getQ(n_seg, n_order, n_costorder,ts);
    
    %#####################################################
    % STEP 2: 计算对应的约束矩阵A_beq
    [Aeq, beq] = getAbeq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond);
    
    f = zeros(size(Q,1),1);
    % 求解多项式系数
    poly_coef = quadprog(Q,f,[],[],Aeq, beq);
end