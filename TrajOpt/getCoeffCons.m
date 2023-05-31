function coeff = getCoeffCons(t, n_order, n_inputorder)
% 返回多项式的系数矩阵
coeff=zeros(n_inputorder,n_order+1);
for row =1:n_inputorder
    for coll = row:n_order+1
        % factorial(coll-1)/factorial(coll-row) = n_order的k阶导的系数
        % t^(coll-row) 该项阶次
        coeff(row,coll)=factorial(coll-1)/factorial(coll-row)*t^(coll-row);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coffeicents matrixs 的形式为A
% q0 q1 q2 q3 q4 q5 q6 q7 q8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A   coeff = [1,  1*t,  1*t^2,  1*t^3,  1*t^4,  1*t^5,  1*t^6,  1*t^7;
%              0,  1,    2*t,    3*t^2,  4*t^3,  5*t^4,  6*t^5,  7*t^6;
%              0,  0,    2,      6*t,    12*t^2, 20*t^3, 30*t^4, 42*t^5;
%              0,  0,    0,      6,      24*t,   60*t^2, 120*t^3,210*t^4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B   =  左右翻转A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
