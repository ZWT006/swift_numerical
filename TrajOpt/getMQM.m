function MQM = getMQM(n_seg, n_order, n_costorder, ts)
% Minimun n_costorder Polynimal Trajectory Generation P33 L5.pdf
% 这里的计算应该是没有问题的,根据原本多项式优化的p'Qp
% 添加这里bezier到幂基的矩阵M, p=m*c 转化为这里的 c'm'Qmc
% 如果从多项式的角度转化,考虑n_costorder 和时间缩放因子的ts
    MQM = [];
    % k = segments
    for k = 1:n_seg
        Q_k = zeros(n_order + 1, n_order + 1);
        M_k = BerneteinCoeff(n_order);
        %#####################################################
        % STEP 1.1: 计算第k段的矩阵Q_k 
        t_s=ts(k);
        for i = n_costorder:n_order
            for j = n_costorder:n_order
                % 根据离散的多项式cost计算minisnap应该是固定的(i+j-7) 因为是原多项式求4阶导
                % 根据时间归一化t_s^(2*n_costorder-1);这里可能有问题 这里的cost是对应的uTRu 再积分
                % 另外真实的积分应该是有T(j)^(i+j-7)-T(j-1)^(i+j-7),这里的时间使用的是相对的所以后项应该是0
                % 如果 cost 的阶次由 n_costorder 设定,则Q矩阵大小可变 row = n_costorder * k;
                % row = (n_order + 1 - n_costorder) * k
                %%%%%%%% 这里的计算有两种形式,1)将控制点转化为多项式系数,带入时间缩放因子,再根据多项式优化的公式进行计算
%                 t_k=t_s;
%                 Q_k(i+1,j+1) = factorial(i)/factorial(i-n_costorder)*factorial(j)/factorial(j-n_costorder)/...
%                     (i+j-n_costorder*2+1)*t_k^(i+j-n_costorder*2+1)/t_s^(i+j-1); %根据时间因子转化bernstein basis (1/s)^i i只与最早的n阶多项式的系数所在阶次有关
                %%%%%%%% 这里的计算有两种形式,2)直接根据Bezier的求导特性,将时间缩放因子在控制点处就考虑
                t_k=1;
                Q_k(i+1,j+1) = factorial(i)/factorial(i-n_costorder)*factorial(j)/factorial(j-n_costorder)/...
                    (i+j-n_costorder*2+1)*t_k^(i+j-n_costorder*2+1)/t_s^(2*n_costorder); %根据时间因子转化bernstein basis (1/s)^i i只与最早的n阶多项式的系数所在阶次有关
            end
        end
        % 构建对角矩阵
        MQM = blkdiag(MQM, M_k'*Q_k*M_k);
    end
end

%##########################################################################
% p0  p1  p2      p3  ......  p(n-1)      p(n)
% 1   1/s 1/s^2   1/s^3       1/s^(n-1)   1/s^(n)         pt
% 0   1/s 1/s^2   1/s^3       1/s^(n-1)   1/s^(n)         p't
% 0   0   1/s^2   1/s^3       1/s^(n-1)   1/s^(n)         p''t
% 0   0   0       1/s^3       1/s^(n-1)   1/s^(n)         p'''t
%##########################################################################