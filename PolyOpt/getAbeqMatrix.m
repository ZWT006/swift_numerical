function [Aeq, beq] = getAbeqMatrix(coeffs,segpoly)
% Minimun Snap Trajectory Generation P35 L5.pdf
% 注意 Aeq_start;Aeq_end;Aeq_wp;Aeq_con 的列数相同为总优化变量个数n_seg*(n_order+1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 参数设置
    n_seg = segpoly.seg;
    n_order = segpoly.norder - 1;
    n_inputorder = segpoly.ninput;
    n_dim   = segpoly.Dim;
    ts = segpoly.T;
    start_cond = segpoly.start_cond;
    end_cond   = segpoly.goal_cond;
    waypoints  = segpoly.waypoints;
    n_all_poly = n_seg*(n_order+1)*n_dim;
    coeff_n = n_order+1;
    %#####################################################
    % p,v,a,j 的起点约束, 约束数量与n_inputorder有关 [x,y,q]
    Aeq_start = zeros(n_inputorder*n_dim, n_all_poly);
    beq_start = zeros(n_inputorder*n_dim,1);
    % STEP 2.1: write expression of Aeq_start and beq_start
    for i=1:n_dim
        Aeq_start(1+(i-1)*n_inputorder:i*n_inputorder,1+(i-1)*coeff_n:i*coeff_n) = getCoeffCons(0,n_order,n_inputorder);
        vector =  start_cond(:,i);
        beq_start(1+(i-1)*n_inputorder:i*n_inputorder) =  vector;
    end
    
    %#####################################################
    % p,v,a,j 的终端约束
    Aeq_end = zeros(n_inputorder*n_dim, n_all_poly);
    beq_end = zeros(n_inputorder*n_dim,1);
    t = ts(end);
    for i=1:n_dim
        Aeq_end(1+(i-1)*n_inputorder:i*n_inputorder,(end-(n_dim-i+1)*coeff_n+1):(end-(n_dim-i)*coeff_n)) = getCoeffCons(t,n_order,n_inputorder);
        vector =  end_cond(:,i);
        beq_end(1+(i-1)*n_inputorder:i*n_inputorder) =  vector;
    end
    %#####################################################
    % 中点的位置约束
    Aeq_wp = zeros((n_seg-1)*n_dim, n_all_poly);
    beq_wp = zeros((n_seg-1)*n_dim, 1);
    
    for k = 0:1:n_seg-2
        coeff = getCoeffCons(ts(k+1),n_order,1); % 中间点只固定 pose
        for i = 1:n_dim
            beq_wp(k*n_dim+i, 1) = waypoints(k+2,i);
            % 1:t:t^2:t^3:t^4…… * [p0 p1 p2 p3……]T
            Aeq_wp(k*n_dim+i, 1+(k*n_dim+i-1)*coeff_n:(k*n_dim + i)*coeff_n) = coeff(1, :);  
        end
    end
    
    %#####################################################
    % 连续性约束
    row = n_inputorder*(n_seg-1)*n_dim;
    col = n_all_poly;
    Aeq_con = zeros(row,col);
    beq_con = zeros(row,1);
    
    for k = 0:1:n_seg-2 
        Mpositive = getCoeffCons(ts(k+1),n_order,n_inputorder);
        Mnegative = -getCoeffCons(0,n_order,n_inputorder);
        for i = 0:n_dim-1
            Aeq_con(1+n_inputorder*(k*n_dim+i):n_inputorder*(k*n_dim+i+1),1+coeff_n*(k*n_dim+i):coeff_n*(k*n_dim+i+1)) = Mpositive;
            Aeq_con(1+n_inputorder*(k*n_dim+i):n_inputorder*(k*n_dim+i+1),1+coeff_n*((k+1)*n_dim+i):coeff_n*((k+1)*n_dim+i+1)) = Mnegative;
        end
    end
    
    %#####################################################
    % 构造约束矩阵
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end

