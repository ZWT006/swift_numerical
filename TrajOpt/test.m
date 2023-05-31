%% oval constrain calculate
% close all
% clear
% clc
% vx=1;
% vy=3;
% phi=2;
% a=1;
% b=2;
% Rq = [cos(phi),sin(phi);-sin(phi),cos(phi)];
% v = [vx;vy];
% vB = Rq*v;
% vBx = vB(1);
% vBy = vB(2);
% vB_L2 = vBx^2/a^2+vBy^2/b^2;
% RM=[sin(phi)^2/b^2+cos(phi)^2/a^2,cos(phi)*sin(phi)/a^2-cos(phi)*sin(phi)/b^2;
%     cos(phi)*sin(phi)/a^2-cos(phi)*sin(phi)/b^2,sin(phi)^2/a^2+cos(phi)^2/b^2];
% vT=[vx,vy];
% vB_L1 = vT*RM*v;
% fprintf('vB norm 2= %6.4f; Eq1 \n',vB_L1);
% fprintf('vB norm 2= %6.4f; Eq2 \n',vB_L2);
opt.algorithm = NLOPT_LD_MMA;
opt.lower_bounds = [-inf, 0];
opt.min_objective = @myfunc;
opt.fc = { (@(x) myconstraint(x,2,0)), (@(x) myconstraint(x,-1,1)) };
opt.fc_tol = [1e-8, 1e-8];
opt.xtol_rel = 1e-4;
[xopt, fmin, retcode] = nlopt_optimize(opt, [1.234 5.678]);
