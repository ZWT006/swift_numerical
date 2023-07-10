[n_seg,~]=size(path);
n_seg = n_seg - 1;
%Step2: 根据选择分段重置path/ts/start_cond/goal_cond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_order = 7;
seq_start = 1;
seq_end   = n_seg+1; % n_seg+1
n_seg = seq_end - seq_start;
path_seg = path(seq_start:seq_end,:);
Tv0 = sum(ts(1:seq_start))-ts(seq_start);
QP_ts = ts;
ts_seg   = ts(seq_start:seq_end-1);
T=sum(ts_seg);
dist_seg = dist(seq_start:seq_end-1);
start_cond = squeeze(pathstates(:,:,seq_start))';
goal_cond  = squeeze(pathstates(:,:,seq_end))';
poly_coef_x_seg = poly_coef_x((seq_start-1)*(n_order+1)+1:(seq_end-1)*(n_order+1));
poly_coef_y_seg = poly_coef_y((seq_start-1)*(n_order+1)+1:(seq_end-1)*(n_order+1));
poly_coef_q_seg = poly_coef_q((seq_start-1)*(n_order+1)+1:(seq_end-1)*(n_order+1));

n_order = 8; % 7阶多项式
n_cost  = 4;
n_input = 4;
n_dim   = 3;
% opt 问题的自由度 整个自由度 减去 等式约束 (waypoints连续性约束+waypoints固定点约束)
% opt value = seg * dimension * order - n_dim*(n_seg + 1)*n_input
opt_dof = n_dim*n_seg*n_order - 2*n_input*n_dim - n_input*(n_seg-1)*n_dim - (n_seg-1)*n_dim;

segpoly.norder  = n_order;
segpoly.ncost   = n_cost;
segpoly.seg     = n_seg; %pieceNum
segpoly.Dim     = 3; % 优化变量的维度
segpoly.coeffl  = n_seg * n_order * n_dim;
segpoly.dof     = opt_dof;
segpoly.ninput  = n_input; % 输入的阶次 默认与ncost一致 zeros(1,3)
segpoly.start_cond  = start_cond;
segpoly.goal_cond   = goal_cond;
% segpoly.start_cond  = [path(1,:)    ;0;0;0];
% segpoly.goal_cond   = [path(end,:)  ;0;0;0];
segpoly.waypoints   = path_seg;
segpoly.dists   = dist_seg;
segpoly.T   = ts_seg; %Time Vector
segpoly.R   = 0.8; % theta的权重系数
[Aeq, beq] = getAbeqMatrix([],segpoly);
segpoly.Aeq = Aeq; % Ac=b
segpoly.beq = beq; % 
segpoly.c   = []; % 优化变量 
segpoly.JgdC = [];
segpoly.JgdT = [];
segpoly.Tpow = [];
segpoly.Map  = [];

segpoly.ds = 0.05; % 距离分辨率
segpoly.dt = 0.05; % 时间分辨率
segpoly.d_th = 0.5; % 距离代价的阈值


x00 = ones(segpoly.coeffl,1)*2; % 起始迭代点
for i = 0:n_seg-1
    x00(1+n_order*(i*n_dim):n_order*(i*n_dim+1))   =  poly_coef_x_seg(1+n_order*i:n_order*(i+1));
    x00(1+n_order*(i*n_dim+1):n_order*(i*n_dim+2)) =  poly_coef_y_seg(1+n_order*i:n_order*(i+1));
    x00(1+n_order*(i*n_dim+2):n_order*(i*n_dim+3)) =  poly_coef_q_seg(1+n_order*i:n_order*(i+1));
end

segpoly.coeffs = x00;
lowb = ones(segpoly.coeffl,1)*-500;
upb = ones(segpoly.coeffl,1)*500;
polytraj = PolyTraj(segpoly);
polytraj = polytraj.setTarray(ts_seg);
x00   = polytraj.getReduceOptVelue(x00);
lowb = polytraj.getReduceOptVelue(lowb);
upb  = polytraj.getReduceOptVelue(upb);

% 关于 Equality Constrain 简化优化问题的一些测试
[Re_Aeq,Re_beq] = polytraj.SolveAeqbeq(x00,segpoly);
[Aeq, beq] = getAbeqMatrix([],segpoly);

fprintf("Optimal Problem DOF = %4d \n",segpoly.dof);

rank_Aeq = rank(Aeq);
fprintf("rank of Aeq = %4d \n",rank_Aeq);

rank_Aeq_beq = rank([Aeq,beq]);
fprintf("rank of beq = %4d \n",rank_Aeq_beq);

rank_Re_Aeq = rank(Re_Aeq);
fprintf("rank of Reduce Aeq = %4d \n",rank_Re_Aeq);

rank_Re_Aeq_beq = rank([Re_Aeq,Re_beq]);
fprintf("rank of Reduce beq = %4d \n",rank_Re_Aeq_beq);

% 不同求解方法
fprintf("linsolve \n");
test_coeffs = linsolve(Re_Aeq,Re_beq);
fprintf("inv solve \n");
inv_coeffs = pinv(Re_Aeq)*Re_beq;
fprintf("div solve \n");
div_coeffs = Re_Aeq\Re_beq;
fprintf("LU solve \n");
[Aeq_L,Aeq_U] = lu(Re_Aeq);
figure(1)
subplot(1,3,1)
spy(Aeq_L)
title('L factor')
subplot(1,3,2)
spy(Aeq_U)
title('U factor')
subplot(1,3,3)
spy(Re_Aeq)
title('A factor')

figure(2)
spy(Re_Aeq)
title('Re_Aeq')

% LU_coeffs = 

% lu_solver = luSolver(Re_Aeq,Re_beq,50,50);
% lu_solver = lu_solver.factorizeLU();
% lu_solver = lu_solver.Solve();
% LU_coeffs = lu_solver.x;

% 求解 coeffs 的 A*x=b 结果偏差
fprintf("Equality Constrain bias:\n");
test_beq = Re_Aeq*segpoly.coeffs;
test_bias = test_beq-Re_beq;
% 经过测试这样处理是能够满足 Ax=b 的线性方程的,但是最的Re_Aeq仿真非奇异搞不定呀
test_bias_sqr = sum(test_bias.^2);
fprintf("Quadratic bias: %6.4f\n",test_bias_sqr);

test_beq = Re_Aeq*test_coeffs;
test_bias = test_beq-Re_beq;
test_bias_sqr = sum(test_bias.^2);
fprintf("linsolve bias: %6.4f\n",test_bias_sqr);

test_beq = Re_Aeq*inv_coeffs;
test_bias = test_beq-Re_beq;
test_bias_sqr = sum(test_bias.^2);
fprintf("inv bias: %6.4f\n",test_bias_sqr);

test_beq = Re_Aeq*div_coeffs;
test_bias = test_beq-Re_beq;
test_bias_sqr = sum(test_bias.^2);
fprintf("div bias: %6.4f\n",test_bias_sqr);

% test_beq = Re_Aeq*LU_coeffs;
% test_bias = test_beq-Re_beq;
% test_bias_sqr = sum(test_bias.^2);
% fprintf("LU bias: %6.4f\n",test_bias_sqr);

% 求解的coeffs 偏差 segpoly coeffs
fprintf("segpoly coeffs bias:\n");
test_bias = segpoly.coeffs-test_coeffs;
test_bias_sqr = sum(test_bias.^2);
fprintf("linsolve coeffs bias: %6.4f\n",test_bias_sqr);

test_bias = segpoly.coeffs-inv_coeffs;
test_bias_sqr = sum(test_bias.^2);
fprintf("INV coeffs bias: %6.4f\n",test_bias_sqr);

test_bias = segpoly.coeffs-div_coeffs;
test_bias_sqr = sum(test_bias.^2);
fprintf("DIV coeffs bias: %6.4f\n",test_bias_sqr);

% test_bias = segpoly.coeffs-LU_coeffs;
% test_bias_sqr = sum(test_bias.^2);
% fprintf("LU coeffs bias: %6.4f\n",test_bias_sqr);

% 求解的coeffs 偏差 linsolve
fprintf("linsolve coeffs bias:\n");
test_bias = test_coeffs-inv_coeffs;
test_bias_sqr = sum(test_bias.^2);
fprintf("INV coeffs bias: %6.4f\n",test_bias_sqr);

test_bias = test_coeffs-div_coeffs;
test_bias_sqr = sum(test_bias.^2);
fprintf("DIV coeffs bias: %6.4f\n",test_bias_sqr);

% test_bias = test_coeffs-LU_coeffs;
% test_bias_sqr = sum(test_bias.^2);
% fprintf("LU coeffs bias: %6.4f\n",test_bias_sqr);
