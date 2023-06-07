function [Aeq, beq] = getAbieqMatrix(coeffs,segpoly)
% Minimun Snap Trajectory Generation P35 L5.pdf
% 注意 Aeq_start;Aeq_end;Aeq_wp;Aeq_con 的列数相同为总优化变量个数n_seg*(n_order+1)
% 调整稀疏矩阵 Aeq 的非零元素位置 提高 Ax=b的求解速率
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 参数设置
    n_seg = segpoly.seg;
    n_order = segpoly.norder - 1;
    n_inputorder = segpoly.ninput;
    n_dim   = segpoly.Dim;
    ts = segpoly.T;
    n_all_poly = n_seg*(n_order+1)*n_dim;
    coeff_n = n_order+1;

    pv_max = segpoly.pv_max;
    pa_max = segpoly.pa_max;
    wv_max = segpoly.wv_max;
    wa_max = segpoly.wa_max;

    xy_max = [pv_max;pa_max;pv_max;pa_max;pv_max;pa_max;pv_max;pa_max];
    q_max  = [wv_max;wa_max;wv_max;wa_max;wv_max;wa_max;wv_max;wa_max];

    %#####################################################
    % waypoint 的 dynamic limit 约束
    row = 2*(n_seg-1)*n_dim;
    col = n_all_poly;
    % -max <  每段起点和终点 < max
    Aeq_con = zeros(row*4,col);
    beq_con = zeros(row*4,1);
    
    for k = 0:1:n_seg-2 
        Mpositive = getCoeffCons(ts(k+1),n_order,n_inputorder);
        Mnegative = -getCoeffCons(0,n_order,n_inputorder);
        velaccMpositive = Mpositive(2:3,:); % 正速度和加速度
        velaccMnegative = Mnegative(2:3,:); % 负速度和加速度
        % x ##############################################################
        Aeq_con(1+24*k:2+24*k,1+coeff_n*(k*n_dim+0):coeff_n*(k*n_dim+1)) = velaccMpositive;
        Aeq_con(3+24*k:4+24*k,1+coeff_n*(k*n_dim+0):coeff_n*(k*n_dim+1)) = -velaccMpositive;
        Aeq_con(5+24*k:6+24*k,1+coeff_n*((k+1)*n_dim+0):coeff_n*((k+1)*n_dim+1)) = velaccMnegative;
        Aeq_con(7+24*k:8+24*k,1+coeff_n*((k+1)*n_dim+0):coeff_n*((k+1)*n_dim+1)) = -velaccMnegative;
        % y ##############################################################
        Aeq_con(9+24*k:10+24*k,1+coeff_n*(k*n_dim+1):coeff_n*(k*n_dim+2)) = velaccMpositive;
        Aeq_con(11+24*k:12+24*k,1+coeff_n*(k*n_dim+1):coeff_n*(k*n_dim+2)) = -velaccMpositive;
        Aeq_con(13+24*k:14+24*k,1+coeff_n*((k+1)*n_dim+1):coeff_n*((k+1)*n_dim+2)) = velaccMnegative;
        Aeq_con(15+24*k:16+24*k,1+coeff_n*((k+1)*n_dim+1):coeff_n*((k+1)*n_dim+2)) = -velaccMnegative;
        % q ##############################################################
        Aeq_con(17+24*k:18+24*k,1+coeff_n*(k*n_dim+2):coeff_n*(k*n_dim+3)) = velaccMpositive;
        Aeq_con(19+24*k:20+24*k,1+coeff_n*(k*n_dim+2):coeff_n*(k*n_dim+3)) = -velaccMpositive;
        Aeq_con(21+24*k:22+24*k,1+coeff_n*((k+1)*n_dim+2):coeff_n*((k+1)*n_dim+3)) = velaccMnegative;
        Aeq_con(23+24*k:24+24*k,1+coeff_n*((k+1)*n_dim+2):coeff_n*((k+1)*n_dim+3)) = -velaccMnegative;
        beq_con(1+24*k:24+24*k) = [xy_max;xy_max;q_max];
    end
    %#####################################################
    % 构造约束矩阵
    Aeq = Aeq_con;
    beq = beq_con;
end

