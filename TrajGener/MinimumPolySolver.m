function poly_coef = MinimumPolySolver(waypoints, ts, n_seg, n_order,n_costorder,n_inputorder,OP_structure)
    % 参数结构体
    QP_inequality = OP_structure.QP_inequality;
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
    [Aieq, bieq] = getAbieq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond,OP_structure);
    
    f = zeros(size(Q,1),1);
    options = optimoptions('quadprog','MaxIterations',6000);
    % 求解多项式系数
    if (QP_inequality)
        poly_coef = quadprog(Q,f,Aieq,bieq,Aeq, beq,[],[],[],options);
    else
        poly_coef = quadprog(Q,f,[],[],Aeq, beq,[],[],[],options);
    end
end
