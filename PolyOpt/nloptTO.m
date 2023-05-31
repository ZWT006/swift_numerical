%***************************************
%Author: Wentao Zhang
%Date: 2023-3-11
%E-mail: zwt190315@163.com
%Reference: 
% https://nlopt.readthedocs.io/en/latest/
% https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
% https://zhaohailong.blog.csdn.net/article/details/119928345
% https://blog.csdn.net/qq_32761549/article/details/119911096
%Problem:
% TODO：感觉可以提升或者可以尝试的地方 
%***************************************
% 使用NLopt进行非线性优化,NLopt的基础使用

%%% algorithm的选择
% 命名方式：NLOPT_{G,L}{N,D}_xxxx
% G/L denotes global/local optimization; N/D denotes derivative-free/gradient-based algorithms
%   enum algorithm {
%      GN_DIRECT = 0,
%      GN_DIRECT_L,
%      GN_DIRECT_L_RAND,
%      GN_DIRECT_NOSCAL,
%      GN_DIRECT_L_NOSCAL,
%      GN_DIRECT_L_RAND_NOSCAL,
%      GN_ORIG_DIRECT,
%      GN_ORIG_DIRECT_L,
%      GD_STOGO,
%      GD_STOGO_RAND,
%      LD_LBFGS_NOCEDAL,
%      LD_LBFGS,
%      LN_PRAXIS,
%      LD_VAR1,
%      LD_VAR2,
%      LD_TNEWTON,
%      LD_TNEWTON_RESTART,
%      LD_TNEWTON_PRECOND,
%      LD_TNEWTON_PRECOND_RESTART,
%      GN_CRS2_LM,
%      GN_MLSL,
%      GD_MLSL,
%      GN_MLSL_LDS,
%      GD_MLSL_LDS,
%      LD_MMA,
%      LN_COBYLA,
%      LN_NEWUOA,
%      LN_NEWUOA_BOUND,
%      LN_NELDERMEAD,
%      LN_SBPLX,
%      LN_AUGLAG,
%      LD_AUGLAG,
%      LN_AUGLAG_EQ,
%      LD_AUGLAG_EQ,
%      LN_BOBYQA,
%      GN_ISRES,
%      AUGLAG,
%      AUGLAG_EQ,
%      G_MLSL,
%      G_MLSL_LDS,
%      LD_SLSQP,
%      LD_CCSAQ,
%      GN_ESCH,
%      NUM_ALGORITHMS /*不是一种算法 只是算法的数量*/
%   };

opt.algorithm = NLOPT_LD_MMA;
opt.lower_bounds = [-inf, 0];
opt.min_objective = @myfunc;
% reference：https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#nonlinear-constraints
% fc 是不等式约束
opt.fc = { (@(x) myconstraint(x,2,0)), (@(x) myconstraint(x,-1,1)) };
opt.fc_tol = [1e-8, 1e-8];
opt.xtol_rel = 1e-4;
[xopt, fmin, retcode] = nlopt_optimize(opt, [3.356 5.678]);