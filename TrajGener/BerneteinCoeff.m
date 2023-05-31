function M_c = BerneteinCoeff(n)
%BERNETEINCOEFF 根据控制点的数量生成Bezier curve coefficients
% 将Bernstein Basis转化为 Power Basis(幂基)
% @reference: https://www.qiujiawei.com/bezier-1/
M_c=zeros(n,n);%对数组进行初始化

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n+1
    for k=i:n+1
        M_c(k,i)=factorial(n)/factorial(i-1)/factorial(n-k+1)/factorial(k-i)*(-1)^(k-i);
    end
end
end

