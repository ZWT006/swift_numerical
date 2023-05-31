function poly_coef = MinimumPolySolver(waypoints, ts, n_seg, n_order,n_costorder,n_inputorder)
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
