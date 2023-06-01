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
if (TimeOptimal)
    ts = coeffs(end-n_seg + 1:end); % 最后的n_seg个变量是tau
    coeffs(end-n_seg + 1:end)=[];
    polytraj = segpoly.traj.setCoeffs(coeffs);
    polytraj = polytraj.settauarray(ts);
    segpoly.T = exp(ts);
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
else
    polytraj = segpoly.traj.setCoeffs(coeffs);
end

grad = zeros(segpoly.coeffl,1);
gradt = zeros(n_seg,1);

switch_equacon = segpoly.switch_equacon;

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
        xgrad = zeros(n_order,1);
        ygrad = zeros(n_order,1);
        qgrad = zeros(n_order,1);
        for idj = 0:polytraj.bars(idi)
            vel_I = vel(idj+1,1:2); % 惯性系下的速度
            acc_I = acc(idj+1,1:2); % 惯性系下的加速度
            yaw = pos(idj+1,3);
            Rmatrix = [cos(yaw),sin(yaw);-sin(yaw),cos(yaw)];
            vel_B = Rmatrix * vel_I'; % 载体系下的速度
            acc_B = Rmatrix * acc_I'; % 载体系下的加速度
            vel_B2 = vel_B(1)^2/ORIEN_VEL^2+vel_B(2)^2/VERDIT_VEL^2;
            c(k,1) = vel_B2 - 1;
            k=k+1;
        end
    end
else
    c=[];
end

end

