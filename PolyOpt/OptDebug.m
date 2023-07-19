%%%% 加载 waypoints
close all; clear; clc;
load("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\pathnode.mat");
% for AngleDelta() function
addpath("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\LazyPRM\");
% for QP solver
addpath("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\TrajGener\");

SAVE_CSV = true;

%%%% 计算每段的时间
RATION = 100;
path(:,1)=path(:,1)/RATION;
path(:,2)=path(:,2)/RATION;

[n_seg,~]=size(path);
n_seg = n_seg - 1;
% reference parameters
Vel_factor = 1.4; % reference Linear Velocity  2m/s
W_factor   = 1.4; % reference Angular Velocity rad/s
pv_max = Vel_factor*1.5;
pa_max = Vel_factor*1.5;
wv_max = W_factor*1.5;
wa_max = W_factor*1.5;


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
%%%% 初末段时间的缩放
xvi=0;yvi=0;qvi=0;
ti = ts(1) * max([(Vel_factor-xvi)/Vel_factor,(Vel_factor-yvi)/Vel_factor,(W_factor-qvi)/W_factor]);
ts(1)   = ts(1) + ti;
clear xvi yvi qvi ti
xvf=0;yvf=0;qvf=0;
tf = ts(end) * max([(Vel_factor-xvf)/Vel_factor,(Vel_factor-yvf)/Vel_factor,(W_factor-qvf)/W_factor]);
ts(end) = ts(end) + tf;
clear xvf yvf qvf tf


%%%% 处理角度变化问题
path_q = path(:, 3);
path_m = zeros(n_seg+1,1);
path_m(1)=path_q(1);
for idx=2:n_seg+1
%     path_m(idx)=AngleDelta(path_q(idx-1),path_q(idx));
    path_m(idx)=path_m(idx-1)+AngleDelta(path_q(idx-1),path_q(idx));
end
path(:,3)=path_m;

clear path_q path_m

%%%%%%%% QP Solver 
n_order       = 7;          % 多项式的阶数 自由度为 n_order+1
n_costorder   = 4;          % 最小化的求导阶次 0=posi;1=vel;2=acc;3=jerk;4=snap;
n_inputorder  = 4;          % 输入的阶次 可以理解为segment 之间满足等式约束的阶次

% 二次优化 dynamic limit
QPdynamiclimit = true;
% Quadratic Debiug Structure
OP_structure.QP_inequality = QPdynamiclimit;

%%%% QP 
OP_structure.v_max = pv_max;
OP_structure.a_max = pa_max;
poly_coef_x = MinimumPolySolver(path(:, 1), ts, n_seg, n_order, n_costorder, n_inputorder,OP_structure);
OP_structure.v_max = pv_max;
OP_structure.a_max = pa_max;
poly_coef_y = MinimumPolySolver(path(:, 2), ts, n_seg, n_order, n_costorder, n_inputorder,OP_structure);
OP_structure.v_max = wv_max;
OP_structure.a_max = wa_max;
poly_coef_q = MinimumPolySolver(path(:, 3), ts, n_seg, n_order, n_costorder, n_inputorder,OP_structure);

n_order = 8; % 7阶多项式
n_cost  = 4;
n_input = 4;
n_dim   = 3;
opt_dof = n_dim*n_seg*n_order - 2*n_input*n_dim - n_input*(n_seg-1)*n_dim - (n_seg-1)*n_dim;

QPcoeffs = ones(24*n_seg,1)*2; % 起始迭代点
for i = 0:n_seg-1
    QPcoeffs(1+n_order*(i*n_dim):n_order*(i*n_dim+1))   =  poly_coef_x(1+n_order*i:n_order*(i+1));
    QPcoeffs(1+n_order*(i*n_dim+1):n_order*(i*n_dim+2)) =  poly_coef_y(1+n_order*i:n_order*(i+1));
    QPcoeffs(1+n_order*(i*n_dim+2):n_order*(i*n_dim+3)) =  poly_coef_q(1+n_order*i:n_order*(i+1));
end

%%%% EDF map and polynomial trajectory
map = imread('F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\map\map10.png');
sdfmap = sdfMap(map);
%%%% segpoly
bound_rate = 0.8;
oval_rate  = 0.8;

segpoly.pv_max = pv_max;
segpoly.pa_max = pa_max;
segpoly.wv_max = wv_max;
segpoly.wa_max = wa_max;
segpoly.dyna_rate = bound_rate;
segpoly.oval_rate = oval_rate;
segpoly.ORIEN_VEL = 2;
segpoly.VERDIT_VEL = 1;
segpoly.ds = 0.05; % 距离分辨率
segpoly.dt = 0.05; % 时间分辨率
segpoly.d_th = 0.5; % 距离代价的阈值
% 几种cost的权重
segpoly.lambda_smooth = 0.1;        % default 1     0.1
segpoly.lambda_obstacle =1;%0.01;   % default 0.01  1       1
segpoly.lambda_dynamic = 10;%500;   % default 500   1       10
segpoly.lambda_time = 8000;%3000;   % default 2000  8000    
segpoly.lambda_oval = 0;%10;       % default 10
% oval cost 和 oval constrain 选择一个起作用即可
segpoly.switch_ovalcon = true;
% Nonlinear equality Constrain
segpoly.switch_equacon = true;
% Using Equality Constrain Reduce Debiug DOF 
segpoly.ReduceOptimalValue  = false;

segpoly.norder  = n_order;
segpoly.ncost   = n_cost;
segpoly.seg     = n_seg; %pieceNum
segpoly.Dim     = n_dim; % 优化变量的维度
segpoly.coeffl  = n_seg * n_order * n_dim;
segpoly.dof     = opt_dof;
segpoly.ninput  = n_input; % 输入的阶次 默认与ncost一致 zeros(1,3)
segpoly.start_cond  = [path(1,:)    ;0 0 0;0 0 0;0 0 0];
segpoly.goal_cond   = [path(end,:)  ;0 0 0;0 0 0;0 0 0];
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
segpoly.TimeOptimal = true;

polytraj = PolyTraj(segpoly);
polytraj = polytraj.setTarray(ts);
polytraj = polytraj.setCoeffs(QPcoeffs);

TimeMat = [];
StatesTraj = [];

for idi = 1:n_seg
    [pos,vel,acc] = polytraj.getStates(idi);
    jer = polytraj.getJerk(idi);
    StatesTraj = [StatesTraj;pos,vel,acc,jer];
    timemat = getCoeffCons(ts(idi),7,4);
    TimeMat = [TimeMat;timemat];
end

Discrete = [ts,polytraj.dists,polytraj.bars,polytraj.bardt];
reduceCoeffs = polytraj.getReduceOptVelue(QPcoeffs);

[Re_Aeq,Re_beq] = polytraj.SolveAeqbeq(reduceCoeffs,segpoly);
Re_x = pinv(Re_Aeq)*Re_beq;

segpoly.sdf = sdfmap;
segpoly.traj = polytraj;

% OptVelue = QPcoeffs;
% re_OptVelue = zeros(opt_dof,1);
% fprintf("optIndex : ");
% dof_num = 0;
% for id_seg=0:segpoly.seg-1 
%     if (dof_num >= segpoly.dof)
%         break;
%     end % 5~8 阶次系数
%     for id_dim = 0:segpoly.Dim-1
%         if (dof_num >= segpoly.dof)
%             break;
%         end
%         for id_order=segpoly.ninput+2:segpoly.norder
%             % 当前优化变量的行序号
%             row_num = id_seg*segpoly.Dim*segpoly.norder + id_dim*segpoly.norder + id_order;
%             dof_num = dof_num + 1;
%             re_OptVelue(dof_num) = OptVelue(row_num);
%             fprintf(" %d ",row_num);
%             if (dof_num >= segpoly.dof)
%                 break;
%             end
%         end
%     end
% end
% fprintf("\n");
% clear OptVelue id_order id_seg id_dim

%%%% 扩充 coeffs 为 优化变量
QPcoeffs = [QPcoeffs;log(ts)];

TimePow = [ts';power(ts,2)';power(ts,3)';power(ts,4)';power(ts,5)';power(ts,6)';power(ts,7)'];


[smoCost,smograd]=smoothCost(QPcoeffs,segpoly);
[obsCost,obsgrad]=obstacleCost(QPcoeffs,segpoly);
[dynCost,dyngrad]=dynamicCost(QPcoeffs,segpoly);
[timCost,timgrad]=timeCost(QPcoeffs,segpoly);
[ovaCost,ovagrad]=ovalCost(QPcoeffs,segpoly);




if (SAVE_CSV)
%%%% 保存离散参数
filename = "E:\datas\Swift\Debug\MATLABDiscrete.csv";
TRAJ_DATA = Discrete;
TRAJ_DATA = round(TRAJ_DATA,4);
writematrix(TRAJ_DATA,filename);
%%%% 保存离散状态
filename = "E:\datas\Swift\Debug\MATLABStatesTraj.csv";
TRAJ_DATA = StatesTraj;
TRAJ_DATA = round(TRAJ_DATA,4);
writematrix(TRAJ_DATA,filename);
%%%% 保存时间多次幂
filename = "E:\datas\Swift\Debug\MATLABTimeMat.csv";
TRAJ_DATA = TimeMat;
TRAJ_DATA = round(TRAJ_DATA,4);
writematrix(TRAJ_DATA,filename);
%%%% 保存cost和grad
filename = "E:\datas\Swift\Debug\MATLABCostGrad.csv";
TRAJ_DATA = [[smoCost;0;smograd], [obsCost;0;obsgrad], [dynCost;0;dyngrad], [timCost;0;timgrad], [ovaCost;0;ovagrad]];
TRAJ_DATA = round(TRAJ_DATA,4);
writematrix(TRAJ_DATA,filename);
%%%% 保存Ax=b求解
filename = "E:\datas\Swift\Debug\MATLABMatVec.csv";
TRAJ_DATA = [Re_Aeq,Re_x,Re_beq];
TRAJ_DATA = round(TRAJ_DATA,4);
writematrix(TRAJ_DATA,filename);
%%%% 保存原始矩阵和coeffs
filename = "E:\datas\Swift\Debug\MATLABAbeq.csv";
TRAJ_DATA = [Aeq,beq];
TRAJ_DATA = round(TRAJ_DATA,4);
writematrix(TRAJ_DATA,filename);

% QPcoeffs(end-2:end)=[];
% filename = "E:\datas\Swift\Debug\MATLABQPCoeffs.csv";
% TRAJ_DATA = [QPcoeffs,[poly_coef_x;poly_coef_y;poly_coef_q],coeffs];
% TRAJ_DATA = round(TRAJ_DATA,4);
% writematrix(TRAJ_DATA,filename);
end
