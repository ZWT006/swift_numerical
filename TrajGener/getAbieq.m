function [Aeq,beq]= getAbieq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond,OP_structure)
% Minimun Snap Trajectory Generation P35 L5.pdf
% 注意 Aeq_start;Aeq_end;Aeq_wp;Aeq_con 的列数相同为总优化变量个数n_seg*(n_order+1)
    n_all_poly = n_seg*(n_order+1);
    coeff_n = n_order+1;
    if (nargin == 8)
        v_max = OP_structure.v_max;
        a_max = OP_structure.a_max;
    end
    %#####################################################
    % 连续性约束
    row = n_inputorder*(n_seg-1);
    col = n_all_poly;
    Aeq_con = zeros(row*2,col);
    beq_con = zeros(row*2,1);
    n_inputorder = 8;
    for k = 0:1:n_seg-2 
        Mpositive = getCoeffCons(ts(k+1),n_order,n_inputorder);
        Mnegative = -getCoeffCons(0,n_order,n_inputorder);
        velaccMpositive = Mpositive(2:3,:); % 正速度和加速度
        velaccMnegative = Mnegative(2:3,:); % 负速度和加速度
        velaccMatrix = [velaccMpositive;velaccMnegative];
        Aeq_con(1+n_inputorder*k:4+n_inputorder*k,1+coeff_n*k:coeff_n*(k+1)) = velaccMatrix;
        Aeq_con(5+n_inputorder*k:n_inputorder*(k+1),1+coeff_n*(k+1):coeff_n*(k+2)) = velaccMatrix;
        beq_con(1+n_inputorder*k:n_inputorder*(k+1)) = [v_max;a_max;v_max;a_max;v_max;a_max;v_max;a_max];
    end
    
    %#####################################################
    % 构造约束矩阵
    Aeq = Aeq_con;
    beq = beq_con;
end