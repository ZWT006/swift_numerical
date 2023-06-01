%***************************************
%Author: Wentao Zhang
%Date: 2023-5
%E-mail: zwt190315@163.com
% Problem: 
% 1): seg 较多时无法求解,并且结果比较抽象
% 2): init values 问题,应该得是设置了初始点更快优化呀
% 3): low/up bound 问题,NLopt怕是不提供这功能呀
% 4): SupplyGradient 没有不太行呀优化非常慢
% 5): 
% TODO：需要完善的地方
% 1):[x] 倒数 Obstacle Cost 的 grad 计算
% 2):[x] 利用求解等式约束使优化问题降维
% 3):[x] 添加 Quadratic Programming 用于求解 Init Optimal Value
% 8): Quadratic init value 的超限处理,根据max上限进行缩放
% 9): 利用Equality Constrain使问题降维并去除该约束
% 4): 编写Constrain/Cost验证函数
% 5): 编写 NLopt 优化求解相关函数文件
% 6): 与前面 LazyKinoPRM 的工作拼接 实现完整的工作流程
% 7): OvalConstrain 的效果验证

%***************************************
% Segments Polynomial Trajectory Parameters
addpath("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\LazyPRM\")
addpath("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\TrajGener\")
close all;
clear;clc;

% map 加载地图与 LazyKinoPRM 中加载地图保持一致
map = imread('F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\map\map0.png');
load("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\path.mat");
RATION = 100;
path(:,1)=path(:,1)/RATION;
path(:,2)=path(:,2)/RATION;
[n_seg,~]=size(path);
n_seg = n_seg - 1;
% reference parameters
Vel_factor = 2; % reference Linear Velocity  2m/s
W_factor   = 2; % reference Angular Velocity rad/s
pv_max = Vel_factor*1.5;
pa_max = Vel_factor*6;
wv_max = W_factor*1.5;
wa_max = W_factor*6;
% Vel_factor = 3; % reference Linear Velocity  2m/s
% W_factor   = 2.5; % reference Angular Velocity rad/s

dist= zeros(n_seg, 1);
ts  = ones(n_seg, 1)*0.8;
% 计算参考线速度时间
for i = 1:n_seg
    dist(i) = sqrt((path(i+1, 1)-path(i, 1))^2 + (path(i+1, 2) - path(i, 2))^2);
    ts(i) = dist(i)/Vel_factor;
end
% 计算参考角速度时间/并选择较大值
for i = 1:n_seg
    t_temp = AngleDelta(path(i+1, 3),path(i,3))/W_factor;
    if (t_temp > ts(i))
        ts(i) = t_temp;
    end
end
clear t_temp;
%%%%%%% attention 这两个起始的时间非常影响优化的结果可以根据初末速度来设置一下
xvi=0;yvi=0;qvi=0;
ti = ts(1) * max([(Vel_factor-xvi)/Vel_factor,(Vel_factor-yvi)/Vel_factor,(W_factor-qvi)/W_factor]);
ts(1)   = ts(1) + ti;
xvf=0;yvf=0;qvf=0;
tf = ts(end) * max([(Vel_factor-xvf)/Vel_factor,(Vel_factor-yvf)/Vel_factor,(W_factor-qvf)/W_factor]);
ts(end) = ts(end) + tf;
T=sum(ts);

%##########################################################################
%Step1: 使用QP求解器求解多项式系数 获得初始点
%%%%%%%%%%%%%%%%
% 处理角度变化问题
path_q = path(:, 3);
path_deg=rad2deg(path_q);
path_m = zeros(n_seg+1,1);
path_m(1)=path_q(1);
for idx=2:n_seg+1
%     path_m(idx)=AngleDelta(path_q(idx-1),path_q(idx));
    path_m(idx)=path_m(idx-1)+AngleDelta(path_q(idx-1),path_q(idx));
end
path_deg_m=rad2deg(path_m);

path(:,3)=path_m;
%%%%%%%%
n_order       = 7;          % 多项式的阶数 自由度为 n_order+1
n_costorder   = 4;          % 最小化的求导阶次 0=posi;1=vel;2=acc;3=jerk;4=snap;
n_inputorder  = 4;          % 输入的阶次 可以理解为segment 之间满足等式约束的阶次

poly_coef_x = MinimumPolySolver(path(:, 1), ts, n_seg, n_order, n_costorder, n_inputorder);
poly_coef_y = MinimumPolySolver(path(:, 2), ts, n_seg, n_order, n_costorder, n_inputorder);
poly_coef_q = MinimumPolySolver(path(:, 3), ts, n_seg, n_order, n_costorder, n_inputorder);

%##########################################################################
pathstates = zeros(3,n_inputorder,n_seg+1);
pathstates(:,:,1) = [path(1,1),0,0,0;path(1,2),0,0,0;path(1,3),0,0,0];
pathstates(:,:,end) = [path(end,1),0,0,0;path(end,2),0,0,0;path(end,3),0,0,0];
for idx=0:n_seg-1
    %#####################################################
    % STEP 3: get the coefficients of i-th segment of both x-axis
    % and y-axis
    Pxi = poly_coef_x((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
    Pyi = poly_coef_y((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
    Pqi = poly_coef_q((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
    px=flip(Pxi);
    py=flip(Pyi);
    pq=flip(Pqi);
    t = ts(idx+1);
    pathstates(1,1,idx+2)  = polyval(px, t);
    pathstates(2,1,idx+2)  = polyval(py, t);
    pathstates(3,1,idx+2)  = polyval(pq, t);
    % velocity
    pdx=polyder(px);
    pdy=polyder(py);
    pdq=polyder(pq);
    pathstates(1,2,idx+2)  = polyval(pdx, t);
    pathstates(2,2,idx+2)  = polyval(pdy, t);
    pathstates(3,2,idx+2)  = polyval(pdq, t);
    % accelaration
    pddx=polyder(pdx);
    pddy=polyder(pdy);
    pddq=polyder(pdq);
    pathstates(1,3,idx+2)  = polyval(pddx, t);
    pathstates(2,3,idx+2)  = polyval(pddy, t);
    pathstates(3,3,idx+2)  = polyval(pddq, t);
    % jerk
    pdddx=polyder(pddx);
    pdddy=polyder(pddy);
    pdddq=polyder(pddq);
    pathstates(1,4,idx+2)  = polyval(pdddx, t);
    pathstates(2,4,idx+2)  = polyval(pdddy, t);
    pathstates(3,4,idx+2)  = polyval(pdddq, t);
end

QPX_n = [];
QPY_n = [];
QPQ_n = [];
QPX_dn = [];
QPY_dn = [];
QPQ_dn = [];
QPX_ddn = [];
QPY_ddn = [];
QPQ_ddn = [];
QP_k=1;
tstep = 0.01;
for idx=0:n_seg-1
    QPPxi = poly_coef_x((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
    QPPyi = poly_coef_y((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
    QPPqi = poly_coef_q((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
    for t = 0:tstep:ts(idx+1)
        % Quadratic Optimization ##########################################
        QPpx=flip(QPPxi);
        QPpy=flip(QPPyi);
        QPpq=flip(QPPqi);
        QPX_n(QP_k)  = polyval(QPpx, t)*RATION;
        QPY_n(QP_k)  = polyval(QPpy, t)*RATION;
        QPQ_n(QP_k)  = polyval(QPpq, t);
        % velocity
        QPpdx=polyder(QPpx);
        QPpdy=polyder(QPpy);
        QPpdq=polyder(QPpq);
        QPX_dn(QP_k)  = polyval(QPpdx, t);
        QPY_dn(QP_k)  = polyval(QPpdy, t);
        QPQ_dn(QP_k)  = polyval(QPpdq, t);
        % accelaration
        QPpddx=polyder(QPpdx);
        QPpddy=polyder(QPpdy);
        QPpddq=polyder(QPpdq);
        QPX_ddn(QP_k)  = polyval(QPpddx, t);
        QPY_ddn(QP_k)  = polyval(QPpddy, t);
        QPQ_ddn(QP_k)  = polyval(QPpddq, t);
        QP_k=QP_k+1;
    end
end
QP_k=QP_k-1;

%Step2: 根据选择分段重置path/ts/start_cond/goal_cond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seq_start = 1;
seq_end   = n_seg+1; % n_seg+1
n_seg = seq_end - seq_start;
path = path(seq_start:seq_end,:);
Tv0 = sum(ts(1:seq_start))-ts(seq_start);
QP_ts = ts;
ts   = ts(seq_start:seq_end-1);
tst  = ts; % 保存未优化的 ts
T=sum(ts);
dist = dist(seq_start:seq_end-1);
start_cond = squeeze(pathstates(:,:,seq_start))';
goal_cond  = squeeze(pathstates(:,:,seq_end))';
poly_coef_x = poly_coef_x((seq_start-1)*(n_order+1)+1:(seq_end-1)*(n_order+1));
poly_coef_y = poly_coef_y((seq_start-1)*(n_order+1)+1:(seq_end-1)*(n_order+1));
poly_coef_q = poly_coef_q((seq_start-1)*(n_order+1)+1:(seq_end-1)*(n_order+1));


segpoly.pv_max = pv_max;
segpoly.pa_max = pa_max;
segpoly.wv_max = wv_max;
segpoly.wa_max = wa_max;
segpoly.dyna_rate = 0.8;

segpoly.oval_rate = 0.8;
segpoly.ORIEN_VEL = 2;
segpoly.VERDIT_VEL = 1;
%##########################################################################
segpoly.DEBUG_PRINT = true;
segpoly.DEBUG_PLOT  = true;
segpoly.TIME_PRINT  = false;
segpoly.CHECK_PLOT  = false;
%%%%%%%%%%%%% GLOBAL DEFINE
MATLAB_SOLVER = 1;
NLOPT_SOLVER  = 2;
OPT_SOLVER = MATLAB_SOLVER;
PLOT_DEBUG = false;
QP_PLOT    = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 优化参数
n = 0;  % 打印figure开始序号
n_order = 8; % 7阶多项式
n_cost  = 4;
n_input = 4;
n_dim   = 3;
% 二次优化初始点
subOptValuesInit = true;
% 时间最优选项
TimeOptimal = true;
% 使用 CostFunction提供的梯度
ObjectiveGradient = true;
% 使用 constant low/up bounds
ConstantBounds    = true;
% 梯度检查选项 正常求解关闭
CheckGradients    = false;
segpoly.ds = 0.05; % 距离分辨率
segpoly.dt = 0.05; % 时间分辨率
segpoly.d_th = 0.5; % 距离代价的阈值
% 几种cost的权重
segpoly.lambda_smooth = 0.1;     % default 1
segpoly.lambda_obstacle =10;%0.01; % default 0.01
segpoly.lambda_dynamic = 1;%500;   % default 500
segpoly.lambda_time = 8000;%3000;     % default 2000 
segpoly.lambda_oval = 0;%10;       % default 10
% oval cost 和 oval constrain 选择一个起作用即可
segpoly.switch_ovalcon = true;
segpoly.switch_equacon = true;

%##########################################################################
%可行参数
% lambda_smooth 1   lambda_obstacle 0.01    lambda_dynamic  600 
% lambda_time   2000    obstacle_cost 无 gradt dynamic_cost 无gradc   

%##########################################################################

segpoly.norder  = n_order;
segpoly.ncost   = n_cost;
segpoly.seg     = n_seg; %pieceNum
segpoly.Dim     = n_dim; % 优化变量的维度
segpoly.coeffl  = n_seg * n_order * n_dim;
segpoly.ninput  = n_input; % 输入的阶次 默认与ncost一致 zeros(1,3)
segpoly.start_cond  = start_cond;
segpoly.goal_cond   = goal_cond;
% segpoly.start_cond  = [path(1,:)    ;0;0;0];
% segpoly.goal_cond   = [path(end,:)  ;0;0;0];
segpoly.waypoints   = path;
segpoly.dists   = dist;
segpoly.T   = ts; %Time Vector
segpoly.R   = 0.8; % theta的权重系数
[Aeq, beq] = getAbeqMatrix([],segpoly);
segpoly.Aeq = Aeq; % Ac=b
segpoly.beq = beq; % 
segpoly.c   = []; % 优化变量 
segpoly.JgdC = [];
segpoly.JgdT = [];
segpoly.Tpow = [];
segpoly.Map  = [];

% 时间最优选项
segpoly.TimeOptimal = TimeOptimal;

%#########################################################################%
% sdf 地图信息

n=n+1;
sdfmap = sdfMap(map);
fp=figure(n);
sdfmap.showSDFMap(fp);
polytraj = PolyTraj(segpoly);
polytraj = polytraj.setTarray(ts);

%##################### 功能函数配置
segpoly.sdf = sdfmap;
segpoly.traj = polytraj;
n=n+1;
optfp = figure(n);
segpoly.optfp   = optfp;
if (segpoly.DEBUG_PLOT)
global iter costArray;
iter = 0;
costArray = [];
end

x0 = ones(segpoly.coeffl,1)*2; % 起始迭代点
if (subOptValuesInit)
    for i = 0:n_seg-1
        x0(1+n_order*(i*n_dim):n_order*(i*n_dim+1))   =  poly_coef_x(1+n_order*i:n_order*(i+1));
        x0(1+n_order*(i*n_dim+1):n_order*(i*n_dim+2)) =  poly_coef_y(1+n_order*i:n_order*(i+1));
        x0(1+n_order*(i*n_dim+2):n_order*(i*n_dim+3)) =  poly_coef_q(1+n_order*i:n_order*(i+1));
    end
end
% 将 quadprog 的coeffs放入segpoly
segpoly.coeffs = x0;

lowb = ones(segpoly.coeffl,1)*-500;
upb = ones(segpoly.coeffl,1)*500;

if (TimeOptimal)
    [Aeq, beq] = getAbeqMatrix([],segpoly);
    bas = Aeq*x0 - beq;
    basum = sum(bas.^2);
    
    fprintf("Aeq bais sum = %d\n",basum);
    
    if (subOptValuesInit)
        t0 = log(ts);
    else
        t0 = zeros(n_seg,1)-0.5;
    end
    x0 = [x0;t0];
    disp("Init ts datas:");
    disp(ts');
    disp("Init t0 datas:");
    disp(t0');
    lowb = [lowb;ones(n_seg,1)*-3];
    upb  = [upb ;ones(n_seg,1)*1.6];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch OPT_SOLVER

    case MATLAB_SOLVER % MATLAB 非线性优化
% 优化选项
% Algorithm:
% trust-region-reflective %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 不接受非线性约束,可以使用梯度加速求解 'SpecifyObjectiveGradient',true,
% interior-point
% 内点法,具有lbfgs的选项 其实lbfgs就是一种求Hession
% sqp sqp-legacy
% reference：https://ww2.mathworks.cn/help/optim/ug/constrained-nonlinear-optimization-algorithms.html#brnpd5f
options = optimoptions('fmincon','display', 'iter');
options.Algorithm   = "interior-point";
options.MaxFunctionEvaluations = 256000; % Default=3000
options.SpecifyObjectiveGradient = ObjectiveGradient; % 不用梯度算得慢死 甚至算不出来结果
options.CheckGradients           = CheckGradients;
options.EnableFeasibilityMode = true; % 内点法找不到满足约束的解,可笑把初始值搞成非0就行了
options.FiniteDifferenceType  = "central"; % 计算一阶微分 forward | central
options.SubproblemAlgorithm = 'factorization'; % factorization 直接尝试牛顿步 | cp 允许共轭梯度
options.OptimalityTolerance = 3e-6;
options.ConstraintTolerance = 1e-5;
problem.options = options;
problem.solver = 'fmincon';
problem.objective = @(x)CostFunc(x,segpoly); %匿名函数可以使用额外参数
problem.x0 = x0; % 起始迭代点
if (ConstantBounds)
    problem.lb = lowb;
    problem.ub = upb;
end
if (TimeOptimal)
    problem.nonlcon = @(x)nonConstrain(x,segpoly);
    problem.Aeq = [];
    problem.beq = [];
else
    if (segpoly.switch_equacon)
        problem.nonlcon = @(x)nonConstrain(x,segpoly);
        problem.Aeq = [];
        problem.beq = [];
    else
        [Aeq, beq] = getAbeqMatrix([],segpoly);
        problem.Aeq = Aeq;
        problem.beq = beq;
    end
end
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
[coeffs, fval, exitflag, output]=fmincon(problem);

    case NLOPT_SOLVER % NLopt 非线性优化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low-stage BFGS:   NLOPT_LD_TNEWTON_PRECOND_RESTART;NLOPT_LD_TNEWTON_PRECOND
%                   NLOPT_LD_TNEWTON_RESTART        ;NLOPT_LD_TNEWTON 
opt.algorithm = NLOPT_LD_MMA;
opt.min_objective = @(x)smoothCost(x,segpoly); %匿名函数可以使用额外参数
nintlength = ones(segpoly.coeffl,1);
% reference：https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#nonlinear-constraints

% fc 是不等式约束
% opt.fc = { (@(x) myconstraint(x,2,0)), (@(x) myconstraint(x,-1,1)) };
% opt.fc_tol = 1e-8;

% h 是等式约束
opt.h = {@(x) equialConstrain(x,segpoly)}; %匿名函数可以使用额外参数
opt.h_tol = 1e-8;


opt.xtol_rel = (1e-4);
[xopt, fmin, retcode] = nlopt_optimize(opt, x0);

end

[obsCost,obsgrad]=obstacleCost(coeffs,segpoly);
[smoCost,smograd]=smoothCost(coeffs,segpoly);
[dynCost,dyngrad]=dynamicCost(coeffs,segpoly);
[timCost,timgrad]=timeCost(coeffs,segpoly);
[ovaCost,ovagrad]=ovalCost(coeffs,segpoly);

lobsCost = segpoly.lambda_obstacle * obsCost;
lsmoCost = segpoly.lambda_smooth   * smoCost;
ldynCost = segpoly.lambda_dynamic  * dynCost;
ltimCost = segpoly.lambda_time     * timCost;
lovaCost = segpoly.lambda_oval     * ovaCost;

fprintf("smoCost = %8.4f; obsCost = %8.4f; dynCost = %8.4f; ovaCost = %8.4f; timCost = %8.6f \n",lsmoCost,lobsCost,ldynCost,lovaCost,ltimCost);

% polytraj.showSegState();

if (TimeOptimal)
    disp(ts');
    ts = coeffs(end-n_seg+1:end);
    ts = exp(ts);
    disp(ts');
    coeffs(end-n_seg+1:end) = [];
    segpoly.T = ts;
    fprintf("Init T =%3.4f;",T);
    T = sum(ts);
    fprintf("Opt T =%3.4f\n",T);
    [Aeq, beq] = getAbeqMatrix([],segpoly);
    bas = Aeq*coeffs - beq;
    basum = sum(bas.^2);
    fprintf("Aeq bais sum = %d\n",basum);
end

% polynomial trajectory 多项式轨迹
% polytraj = PolyTraj(segpoly);
polytraj = polytraj.setCoeffs(coeffs);
polytraj = polytraj.setTarray(ts);
[pos,vel,acc]=polytraj.getStates();

if (segpoly.DEBUG_PLOT)
figure(optfp)
hold on
xlabel("iter");
ylabel("cost");
plot(costArray(:,1),'r-');
plot(costArray(:,2),'b-');
plot(costArray(:,3),'g-');
plot(costArray(:,4),'y-');
plot(costArray(:,5),'k-');
legend('smoCost','obsCost','dynCost','timCost','ovaCost');
grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 用于显示轨迹
RATION  = 100;
X_n = [];
Y_n = [];
Q_n = [];
X_dn = [];
Y_dn = [];
Q_dn = [];
X_ddn = [];
Y_ddn = [];
Q_ddn = [];
k = 1;
tstep = 0.01;
for idx=0:n_seg-1
    %#####################################################
    % STEP 3: get the coefficients of i-th segment of both x-axis and y-axis
    Pxi = coeffs(n_order*(idx*n_dim)+1            :n_order*(idx*n_dim)+n_order);
    Pyi = coeffs(n_order*(idx*n_dim)+n_order+1    :n_order*(idx*n_dim)+n_order*2);
    Pqi = coeffs(n_order*(idx*n_dim)+n_order*2+1  :n_order*(idx*n_dim)+n_order*3);
    for t = 0:tstep:ts(idx+1)
        % Nonlinear Optimization ##########################################
        px=flip(Pxi);
        py=flip(Pyi);
        pq=flip(Pqi);
        X_n(k)  = polyval(px, t)*RATION;
        Y_n(k)  = polyval(py, t)*RATION;
        Q_n(k)  = polyval(pq, t);
        % velocity
        pdx=polyder(px);
        pdy=polyder(py);
        pdq=polyder(pq);
        X_dn(k)  = polyval(pdx, t);
        Y_dn(k)  = polyval(pdy, t);
        Q_dn(k)  = polyval(pdq, t);
        % accelaration
        pddx=polyder(pdx);
        pddy=polyder(pdy);
        pddq=polyder(pdq);
        X_ddn(k)  = polyval(pddx, t);
        Y_ddn(k)  = polyval(pddy, t);
        Q_ddn(k)  = polyval(pddq, t);
        k = k + 1;
    end
%     fprintf('t = %2.6f;ts(idx) =%2.6f \n',t,ts(idx+1))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xy map & palstance
figure (fp)
% subplot(2,3,iter);
hold on
datename  = " nonlinear";
pic_title = strcat("xy map & palstance",datename);
title(pic_title);
AngUNIT=50;
[path_length,~] = size(path);
for idx=1:path_length
%     plot([path(idx,1),path(idx,1)+cos(path(idx,3))*AngUNIT],[path(idx,2),path(idx,2)+sin(path(idx,3))*AngUNIT],'-','Color','b','LineWidth',1);
    plot([path(idx,1)*RATION,path(idx,1)*RATION+cos(path(idx,3))*AngUNIT],[path(idx,2)*RATION,path(idx,2)*RATION+sin(path(idx,3))*AngUNIT],'-','Color','b','LineWidth',1);
%     plot([path(idx,1)*RATION,path(idx,1)*RATION+cos(path_m(idx))*AngUNIT],[path(idx,2)*RATION,path(idx,2)*RATION+sin(path_m(idx))*AngUNIT],'c-','LineWidth',2);
end
for idx=1:path_length-1
    plot([path(idx,1)*RATION,path(idx+1,1)*RATION],[path(idx,2)*RATION,path(idx+1,2)*RATION],'g-','LineWidth',1);
end
AngUNIT=5;
k = k - 1;
for idx=1:k
    detl = Q_dn(idx)*AngUNIT;
    detl = abs(detl);
    theta = Q_n(idx);
    plot([X_n(idx),X_n(idx)+cos(theta)*detl],[Y_n(idx),Y_n(idx)+sin(theta)*detl],'b-','LineWidth',1);
end
plot(X_n, Y_n , 'r-');
% scatter(path(1:size(path, 1), 1), path(1:size(path, 1), 2),'ok');
scatter(path(1:size(path, 1), 1)*RATION, path(1:size(path, 1), 2)*RATION,'ob');
if (QP_PLOT)
    plot(QPX_n,QPY_n,'m--');
end
scatter(pos(:,1)*RATION,pos(:,2)*RATION,25,'.k');

grid on
hold off
axis equal
xlim([0 sdfmap.rows]);
ylim([0 sdfmap.cols]);
set(gcf,'Position', [100, 100, 100+sdfmap.rows, 100+sdfmap.cols]);

T = sum(QP_ts);
if(PLOT_DEBUG) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT_DEBUG
% xy velocity
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("xy velocity",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
tv = tv + Tv0;
plot(tv,X_dn, 'r-');
plot(tv,Y_dn, 'b-');
if (QP_PLOT)
    tv=0:tstep:(QP_k-1)*tstep;
    plot(tv,QPX_dn, 'r--');
    plot(tv,QPY_dn, 'b--');
end
t_temp=0;
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp+Tv0,X_dn(qdn_idx),'*r');
    scatter(t_temp+Tv0,Y_dn(qdn_idx),'*b');
end
%##############################################
if (QP_PLOT)
    t_temp=0;
    for idx=1:length(QP_ts)
        t_temp = t_temp + QP_ts(idx);
        qdn_idx = ceil(t_temp/tstep);
        scatter(t_temp,QPX_dn(qdn_idx),'xr');
        scatter(t_temp,QPY_dn(qdn_idx),'xb');
    end
end
% legend('x vel','y vel');
grid on
xlim([0 ceil(T)]);
% ylim([-150 250]);
ylim([floor(min(Y_dn))-2 ceil(max(Y_dn))+2]);
line([0,ceil(T)],[pv_max,pv_max],'linestyle','--','color','c','LineWidth',2);
line([0,ceil(T)],[-pv_max,-pv_max],'linestyle','--','color','c','LineWidth',2);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% angular velocity
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("angular velocity",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
tv = tv + Tv0;
plot(tv,Q_dn, 'b-');
if(QP_PLOT)
    tv=0:tstep:(QP_k-1)*tstep;
    plot(tv,QPQ_dn, 'b--');
end
t_temp = 0;
scatter(t_temp+Tv0,Q_dn(1),'*r');
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp+Tv0,Q_dn(qdn_idx),'*r');
end
%##############################################
if(QP_PLOT)
    t_temp=0;
    scatter(t_temp,QPQ_dn(1),'*r');
    for idx=1:length(QP_ts)
        t_temp = t_temp + QP_ts(idx);
        qdn_idx = ceil(t_temp/tstep);
        scatter(t_temp,QPQ_dn(qdn_idx),'xr');
    end
end
grid on
xlim([0 ceil(T)]);
ylim([floor(min(Q_dn)) ceil(max(Q_dn))]);
line([0,ceil(T)],[wv_max,wv_max],'linestyle','--','color','c','LineWidth',2);
line([0,ceil(T)],[-wv_max,-wv_max],'linestyle','--','color','c','LineWidth',2);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% angular position
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("angular position",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
tv = tv + Tv0;
plot(tv,rad2deg(Q_n), 'b-');
%##############################################
if(QP_PLOT)
    tv=0:tstep:(QP_k-1)*tstep;
    plot(tv,rad2deg(QPQ_n), 'b--');
end
t_temp = 0;
path_m = path(:,3);
scatter(t_temp+Tv0,rad2deg(path_m(1)),'*r');
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    scatter(t_temp+Tv0,rad2deg(path_m(idx+1)),'*r');
end
grid on
xlim([0 ceil(T)]);
ylim([-200 200]);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% xy acceleration
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("xy acceleration",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
tv = tv + Tv0;
plot(tv,X_ddn, 'r-');
plot(tv,Y_ddn, 'b-');
if (QP_PLOT)
    tv=0:tstep:(QP_k-1)*tstep;
    plot(tv,QPX_ddn, 'r--');
    plot(tv,QPY_ddn, 'b--');
end
t_temp=0;
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp+Tv0,X_ddn(qdn_idx),'*r');
    scatter(t_temp+Tv0,Y_ddn(qdn_idx),'*b');
end
%##############################################
if (QP_PLOT)
    t_temp=0;
    for idx=1:length(QP_ts)
        t_temp = t_temp + QP_ts(idx);
        qdn_idx = ceil(t_temp/tstep);
        scatter(t_temp,QPX_ddn(qdn_idx),'xr');
        scatter(t_temp,QPY_ddn(qdn_idx),'xb');
    end
end
% legend('x acc','y acc');
grid on
xlim([0 ceil(T)]);
ylim([floor(min(Y_ddn))-2 ceil(max(Y_ddn))+2]);
line([0,ceil(T)],[pa_max,pa_max],'linestyle','--','color','c','LineWidth',2);
line([0,ceil(T)],[-pa_max,-pa_max],'linestyle','--','color','c','LineWidth',2);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% angular acceleration
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("angular acceleration",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
tv = tv + Tv0;
plot(tv,Q_ddn, 'b-');
if (QP_PLOT)
    tv=0:tstep:(QP_k-1)*tstep;
    plot(tv,QPQ_ddn, 'b--');
end
t_temp = 0;
scatter(t_temp+Tv0,Q_ddn(1),'*r');
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp+Tv0,Q_ddn(qdn_idx),'*r');
end
%##############################################
if (QP_PLOT)
    t_temp=0;
    scatter(t_temp,Q_ddn(1),'*r');
    for idx=1:length(QP_ts)
        t_temp = t_temp + QP_ts(idx);
        qdn_idx = ceil(t_temp/tstep);
        scatter(t_temp,QPQ_ddn(qdn_idx),'xr');
    end
end
grid on
xlim([0 ceil(T)]);
ylim([floor(min(Q_ddn)) ceil(max(Q_ddn))]);
line([0,ceil(T)],[wa_max,wa_max],'linestyle','--','color','c','LineWidth',2);
line([0,ceil(T)],[-wa_max,-wa_max],'linestyle','--','color','c','LineWidth',2);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT_DEBUG
