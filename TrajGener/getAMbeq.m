function [Aeq,beq]= getAMbeq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond)
% Minimun Snap Trajectory Generation P35 L5.pdf ;Bezier Polynomial
% 注意 Aeq_start;Aeq_end;Aeq_wp;Aeq_con 的列数相同为总优化变量个数n_seg*(n_order+1)
    n_all_poly = n_seg*(n_order+1);
    coeff_n = n_order+1;
    %#####################################################
    % p,v,a,j 的起点约束, 约束数量与n_inputorder有关
    Aeq_start = zeros(n_inputorder, n_all_poly);
    % STEP 2.1: write expression of Aeq_start and beq_start
    polyCoeff =  getCoeffCons(0,n_order,n_inputorder);
    bezierCoeff = polyCoeff*BerneteinCoeff(n_order);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%处理贝塞尔曲线时间归一化的问题
    t = ts(1);
    power=0:1:n_order;
    scal_t=t.^power;
    Aeq_start(1:n_inputorder,1:coeff_n) = bezierCoeff./scal_t;
    beq_start =  start_cond';
    
    %#####################################################
    % p,v,a,j 的终端约束
    t = ts(end);
    Aeq_end = zeros(n_inputorder, n_all_poly);
    polyCoeff =  getCoeffCons(t,n_order,n_inputorder);
    bezierCoeff = polyCoeff*BerneteinCoeff(n_order);
    power=0:1:n_order;
    scal_t=t.^power;
    Aeq_end(1:n_inputorder,(end-coeff_n+1):end) = bezierCoeff./scal_t;
    beq_end =  end_cond';
    
    %#####################################################
    % 中点的位置约束
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    
    for k = 0:1:n_seg-2
        beq_wp(k+1, 1) = waypoints(k+2);
        t=ts(k+1);
%         coeff = getCoeffCons(ts(k+1),n_order,1); % 中间点只固定 pose
        % 注意bezier曲线的时间为[0 1]
        polyCoeff =  getCoeffCons(t,n_order,n_inputorder);
        bezierCoeff = polyCoeff*BerneteinCoeff(n_order);
        power=0:1:n_order;
        scal_t=t.^power;
        % 1:t:t^2:t^3:t^4…… * [p0 p1 p2 p3……]T
        Aeq_wp(k+1, 1+k*coeff_n:(1+k)*coeff_n) = bezierCoeff(1, :)./scal_t;
    end
    beq_wp=waypoints(2:end-1);
    %#####################################################
    % 连续性约束
%     Aeq_con_p = zeros(n_seg-1, n_all_poly);
%     beq_con_p = zeros(n_seg-1, 1);
%     Aeq_con_v = zeros(n_seg-1, n_all_poly);
%     beq_con_v = zeros(n_seg-1, 1);
%     Aeq_con_a = zeros(n_seg-1, n_all_poly);
%     beq_con_a = zeros(n_seg-1, 1);
%     Aeq_con_j = zeros(n_seg-1, n_all_poly);
%     beq_con_j = zeros(n_seg-1, 1);
%     Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
%     beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    row = n_inputorder*(n_seg-1);
    col = n_all_poly;
    Aeq_con = zeros(row,col);
    beq_con = zeros(row,1);
    power=0:1:n_order;
    
    for k = 0:1:n_seg-2 
        t=ts(k+1);
        polyCoeff =  getCoeffCons(t,n_order,n_inputorder);
        bezierCoeff = polyCoeff*BerneteinCoeff(n_order);
%         bezierCoeff(1,:)=[];
        scal_t=t.^power;
        Aeq_con(1+n_inputorder*k:n_inputorder*(k+1),1+coeff_n*k:coeff_n*(k+1)) =bezierCoeff./scal_t;
%         Aeq_con(1+n_inputorder*k:n_inputorder*(k+1),1+coeff_n*k:coeff_n*(k+1)) = getCoeffCons(ts(k+1),n_order,n_inputorder);
        t=ts(k+2);
        polyCoeff =  getCoeffCons(0,n_order,n_inputorder);
        bezierCoeff = polyCoeff*BerneteinCoeff(n_order);
%         bezierCoeff(1,:)=[];
        scal_t=t.^power;
        Aeq_con(1+n_inputorder*k:n_inputorder*(k+1),1+coeff_n*(k+1):coeff_n*(k+2)) = -bezierCoeff./scal_t;
%         Aeq_con(1+n_inputorder*k:n_inputorder*(k+1),1+coeff_n*(k+1):coeff_n*(k+2)) = -getCoeffCons(0,n_order,n_inputorder);            
    end
    
    %#####################################################
    % 构造约束矩阵
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end