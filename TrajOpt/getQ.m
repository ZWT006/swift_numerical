function Q = getQ(n_seg, n_order, n_costorder, ts)
% Minimun n_costorder Polynimal Trajectory Generation P33 L5.pdf
% 
    Q = [];
    % k = segments
    for k = 1:n_seg
        Q_k = zeros(n_order + 1, n_order + 1);
        %#####################################################
        % STEP 1.1: 计算第k段的矩阵Q_k 
        for i = n_costorder:n_order
            for j = n_costorder:n_order
                % 根据离散的多项式cost计算minisnap应该是固定的(i+j-7) 因为是原多项式求4阶导
                % 另外真实的积分应该是有T(j)^(i+j-7)-T(j-1)^(i+j-7),这里的时间使用的是相对的所以后项应该是0
                % 如果 cost 的阶次由 n_costorder 设定,则Q矩阵大小可变 row = n_costorder * k;
                % row = (n_order + 1 - n_costorder) * k
                Q_k(i+1,j+1) = factorial(i)/factorial(i-n_costorder)*factorial(j)/factorial(j-n_costorder)/(i+j-n_costorder*2+1)*ts(k)^(i+j-n_costorder*2+1);
            end
        end
        % 构建对角矩阵
        Q = blkdiag(Q, Q_k);
    end
end