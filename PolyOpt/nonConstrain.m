function [c,ceq] = nonConstrain(coeffs,segpoly)
%NONCONSTRAIN 非线性约束
%  只要涉及到与时间有关的多项式优化问题,不可避免就变成了非线性的了
% 1. 等式约束 waypoints的状态连续性约束
% 2. 等式约束 waypoints的pos固定约束
n_seg = segpoly.seg;
dim = segpoly.Dim;
n_order = segpoly.norder;

%##########################################################################
% Aeq*x=beq nonlinear equality constraints 
TimeOptimal = segpoly.TimeOptimal;
TIME_PRINT = segpoly.TIME_PRINT;

% 等式约束降维选项
ReduceOptimalValue = segpoly.ReduceOptimalValue;


if (TimeOptimal)
    ts = coeffs(end-n_seg + 1:end); % 最后的n_seg个变量是tau
    coeffs(end-n_seg + 1:end)=[];
    polytraj = segpoly.traj.settauarray(ts);
    ts = exp(ts);
    segpoly.T = ts;
    T = segpoly.T;
    if (TIME_PRINT)
        fprintf("Init ts datas:");
        disp(ts');
        fprintf("Init T  datas:");
        disp(T');
        Tmin = min(T);
        if (Tmin <= 0.05)
            fprintf("<=======>ts min = %2.4f \n",Tmin);
        end
    end
end
switch_equacon = segpoly.switch_equacon;
ts = segpoly.T;
% 使用 ReduceOptimalValue 后默认关闭非线性等式约束
if (ReduceOptimalValue)
    segpoly.T = ts;
    coeffs = segpoly.traj.SolveCoeffs(coeffs,segpoly);
    switch_equacon = false;
end

polytraj = segpoly.traj.setCoeffs(coeffs);

grad = zeros(segpoly.coeffl,1);
gradt = zeros(n_seg,1);

if (switch_equacon)
    [Aeq, beq] = getAbeqMatrix([],segpoly);
    ceq = Aeq*coeffs - beq;
else
    ceq = [];
end

% bias = coeffs - segpoly.coeffs;
% biassum = sum(bias.^2);
% fprintf("coeffes bias sum = %d\n",biassum);
% cepbasum = sum(ceq.^2);
% fprintf("ceq bias sum = %d\n",cepbasum);

%##########################################################################
% oval velocity constraints nonlinear inequality constraints

switch_ovalcon = segpoly.switch_ovalcon;
ORIEN_VEL = segpoly.ORIEN_VEL;
VERDIT_VEL = segpoly.VERDIT_VEL;
if (switch_ovalcon)
    c = zeros(sum(segpoly.traj.bars)+n_seg,1);
    k = 1;
    for idi = 1:n_seg
        [pos,vel,acc] = polytraj.getStates(idi);
        Tdt = polytraj.bardt(n_seg);
        xgrad = zeros(n_order,1);
        ygrad = zeros(n_order,1);
        qgrad = zeros(n_order,1);
        for idj = 0:polytraj.bars(idi)
            Mt = getCoeffCons(idi*Tdt,n_order,4);
            vect1 = Mt(1,:);
            vect2 = Mt(2,:);
            vel_I = vel(idj+1,1:2); % 惯性系下的速度
            acc_I = acc(idj+1,1:2); % 惯性系下的加速度
            yaw = pos(idj+1,3);
            Rmatrix = [cos(yaw),sin(yaw);-sin(yaw),cos(yaw)];
            vel_B = Rmatrix * vel_I'; % 载体系下的速度
            acc_B = Rmatrix * acc_I'; % 载体系下的加速度
            vel_B2 = vel_B(1)^2/ORIEN_VEL^2+vel_B(2)^2/VERDIT_VEL^2;
            c(k,1) = vel_B2 - 1;
            k=k+1;
            gradt(idi) = gradt(idi) + ...
                (2*vel_B(2)*(acc_B(2) - cos(yaw)*vel(1)*vel(3) - sin(yaw)*vel(2)*vel(3))/VERDIT_VEL^2 + ...
                2*vel_B(1)*(acc_B(1) - sin(yaw)*vel(1)*vel(3) + cos(yaw)*vel(2)*vel(3))/ORIEN_VEL^2) * idi * Tdt;
            xgrad = xgrad + ...
                (2*cos(yaw)*vel_B(1)/ORIEN_VEL^2 - 2*sin(yaw)*vel_B(2)/VERDIT_VEL^2) * vect2';
            ygrad = ygrad + ...
                (2*sin(yaw)*vel_B(1)/ORIEN_VEL^2 + 2*cos(yaw)*vel_B(2)/VERDIT_VEL^2) * vect2';
            qgrad = qgrad + ...
                (2*vel_B(1)*vel_B(2)/ORIEN_VEL^2 + 2*(-vel_B(1))*vel_B(2)/VERDIT_VEL^2) * vect1';
        end
        grad((idi-1)*dim*n_order + 1:((idi-1)*dim+1)*n_order)     = xgrad;
        grad(((idi-1)*dim+1)*n_order + 1:((idi-1)*dim+2)*n_order) = ygrad;
        grad(((idi-1)*dim+2)*n_order + 1:((idi-1)*dim+3)*n_order) = qgrad;
    end

else
    c=[];
end
if(TimeOptimal)
    grad = [grad;gradt];
end
end

%% 多项式系数矩阵：由时间t和多项式阶次,输入阶次计算多项式系数矩阵
function coeff = getCoeffCons(t, n_order, n_inputorder)
    % 返回多项式的系数矩阵
    coeff=zeros(n_inputorder,n_order);
    for row =1:n_inputorder
        for coll = row:n_order
            % factorial(coll-1)/factorial(coll-row) = n_order的k阶导的系数
            % t^(coll-row) 该项阶次
            coeff(row,coll)=factorial(coll-1)/factorial(coll-row)*t^(coll-row);
        end
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % coffeicents matrixs 的形式为A
    % q0 q1 q2 q3 q4 q5 q6 q7 q8
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     coeff = [1,  1*t,  1*t^2,  1*t^3,  1*t^4,  1*t^5,  1*t^6,  1*t^7;
    %              0,  1,    2*t,    3*t^2,  4*t^3,  5*t^4,  6*t^5,  7*t^6;
    %              0,  0,    2,      6*t,    12*t^2, 20*t^3, 30*t^4, 42*t^5;
    %              0,  0,    0,      6,      24*t,   60*t^2, 120*t^3,210*t^4];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
