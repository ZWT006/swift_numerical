function [Aieq, bieq] = getAMbieq(n_seg, n_order, ts, v_max, a_max)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % STEP 3.2.1 位置约束 暂时不需要

%     v = ones(n_all_poly,1);
%     for i = 0:n_seg - 1
%     % 为啥这里空间位置要乘以速度矢量因子呢？贝塞尔曲线的t的变化范围不是[0 1]吗？
%         v(i *(n_order + 1)+1:i *(n_order + 1)+n_order + 1) = v(i *(n_order + 1)+1:i *(n_order + 1)+n_order + 1) * ts(i+1);
%     end
%     
%     
%     Aieq_p = diag(v);
%     bieq_p = zeros(n_all_poly,1);
%     for i = 0:n_seg - 1
%         bieq_p(i *(n_order + 1)+1:i *(n_order + 1)+n_order + 1) = corridor_range(i+1,2);
%     end
%     
%     Aieq_p = [Aieq_p;-Aieq_p];
%     bieq_p = [bieq_p;bieq_p];
%     
%     for i = 0:n_seg - 1
%         bieq_p(n_all_poly + i *(n_order + 1)+1:n_all_poly + i *(n_order + 1)+n_order + 1)...
%             = (-1) * corridor_range(i + 1,1);
%     end
%     
% 
    %#####################################################
    % STEP 3.2.2 v constraint   
    Aieq_v = zeros(n_order * n_seg, n_all_poly);
    j = 1;
    for i = 1:n_order * n_seg               %这里速度应该除以时间缩放因子的1次方
        Aieq_v(i,j:j+1) = n_order * [-1, 1]/ts(floor(j/(n_order + 1)) + 1);
        if mod(j + 1, n_order + 1) == 0
            % 新段,挪动两位从c0开始
            j = j + 2;
        else
            % 单段区间内遍历c0~cn
            j = j + 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -1   1   0   0   0   0   0   0   0   0   0   0
    %  0  -1   1   0   0   0   0   0   0   0   0   0
    %  0   0  -1   1   0   0   0   0   0   0   0   0
    %  0   0   0  -1   1   0   0   0   0   0   0   0   0
    %  0   0   0   0  -1   1   0   0   0   0   0   0   0   0
    %  0   0   0   0   0  -1   1   0   0   0   0   0   0   0   0
    %  0   0   0   0   0   0  -1   1   0   0   0   0   0   0   0   0
    %  0   0   0   0   0   0   0   0  -1   1   0   0   0   0   0   0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Aieq_v = [Aieq_v;-Aieq_v];
    bieq_v = ones(2 * n_order * n_seg,1)* v_max;

    %#####################################################
    % STEP 3.2.3 a constraint   
    Aieq_a = zeros((n_order - 1) * n_seg, n_all_poly);
    j = 1;
    for i = 1:(n_order - 1) * n_seg                         %感觉这里加速度应该除以速度缩放因子的平方呀
        Aieq_a(i,j:j+2) = n_order * (n_order - 1)*[1, -2, 1]/ts(floor(j/(n_order + 1)) + 1)^2;
        if mod(j + 2, n_order + 1) == 0
            j = j + 3;
        else
            j = j + 1;
        end
    end
    
    Aieq_a = [Aieq_a;-Aieq_a];
    bieq_a = ones((n_order - 1) * n_seg * 2,1)*a_max;
    %#####################################################
    % combine all components to form Aieq and bieq   
%     Aieq = [Aieq_p; Aieq_v; Aieq_a];
%     bieq = [bieq_p; bieq_v; bieq_a];
    Aieq = [Aieq_v; Aieq_a];
    bieq = [bieq_v; bieq_a];
end