function [cost,grad] = ovalCost(coeffs,segpoly)
%OVALCOST 此处显示有关此函数的摘要
%   此处显示详细说明
n_seg = segpoly.seg;
dim = segpoly.Dim;
n_order = segpoly.norder;
ORIEN_VEL = segpoly.ORIEN_VEL;
VERDIT_VEL = segpoly.VERDIT_VEL;
oval_rate = segpoly.oval_rate;

%##########################################################################
TimeOptimal = segpoly.TimeOptimal;
if (TimeOptimal)
    ts = coeffs(end-n_seg + 1:end); % 最后的n_seg个变量是tau
    coeffs(end-n_seg + 1:end)=[];
    polytraj = segpoly.traj.setCoeffs(coeffs);
    polytraj = polytraj.settauarray(ts);
else
    polytraj = segpoly.traj.setCoeffs(coeffs);
end

cost = 0;
grad = zeros(segpoly.coeffl,1);
gradt = zeros(n_seg,1);

for idi = 1:n_seg
    [pos,vel,acc] = polytraj.getStates(idi);
    Tdt = polytraj.bardt(n_seg);
    xgrad = zeros(n_order,1);
    ygrad = zeros(n_order,1);
    qgrad = zeros(n_order,1);
    for idj = 0:polytraj.bars(idi)
        vel_I = vel(1:2); % 惯性系下的速度
        acc_I = acc(1:2); % 惯性系下的加速度
        Mt = getCoeffCons(idi*Tdt,n_order,4);
        vect1 = Mt(1,:);
        vect2 = Mt(2,:);
        yaw = pos(3);
        Rmatrix = [cos(yaw),sin(yaw);-sin(yaw),cos(yaw)];
        vel_B = Rmatrix * vel_I'; % 载体系下的速度
        acc_B = Rmatrix * acc_I'; % 载体系下的加速度
        vel_B2 = vel_B(1)^2/ORIEN_VEL^2+vel_B(2)^2/VERDIT_VEL^2;
        if (vel_B2 > oval_rate)
            cost = cost + vel_B2 - oval_rate;
% \left\{ \begin{array}{c}
% 	\frac{2(\mathrm{vx(}t)\cos\mathrm{(}q(t))+\mathrm{vy(}t)\sin\mathrm{(}q(t)))\left( \mathrm{vx(}t)q\prime(t)(-\sin\mathrm{(}q(t)))+\mathrm{vy(}t)q\prime(t)\cos\mathrm{(}q(t))+\cos\mathrm{(}q(t))\mathrm{vx}\prime(t)+\sin\mathrm{(}q(t))\mathrm{vy}\prime(t) \right)}{\mathrm{vela}}+\\
% 	\frac{2(\mathrm{vy(}t)\cos\mathrm{(}q(t))-\mathrm{vx(}t)\sin\mathrm{(}q(t)))\left( \mathrm{vx(}t)q\prime(t)(-\cos\mathrm{(}q(t)))-\mathrm{vy(}t)q\prime(t)\sin\mathrm{(}q(t))-\sin\mathrm{(}q(t))\mathrm{vx}\prime(t)+\cos\mathrm{(}q(t))\mathrm{vy}\prime(t) \right)}{\mathrm{velb}}\\
% \end{array} \right\} 
            gradt(idi) = gradt(idi) + ...
                (2*vel_B(2)*(acc_B(2) - cos(yaw)*vel(1)*vel(3) - sin(yaw)*vel(2)*vel(3))/VERDIT_VEL^2 + ...
                2*vel_B(1)*(acc_B(1) - sin(yaw)*vel(1)*vel(3) + cos(yaw)*vel(2)*vel(3))/ORIEN_VEL^2) * idi * Tdt;
% \left\{ \begin{array}{c}
% 	\frac{2\cos\mathrm{(}q(t))(\mathrm{vx(}t)\cos\mathrm{(}q(t))+\mathrm{vy(}t)\sin\mathrm{(}q(t)))}{\mathrm{vela}}-\frac{2\sin\mathrm{(}q(t))(\mathrm{vy(}t)\cos\mathrm{(}q(t))-\mathrm{vx(}t)\sin\mathrm{(}q(t)))}{\mathrm{velb}}\\
% 	\frac{2\sin\mathrm{(}q(t))(\mathrm{vx(}t)\cos\mathrm{(}q(t))+\mathrm{vy(}t)\sin\mathrm{(}q(t)))}{\mathrm{vela}}+\frac{2\cos\mathrm{(}q(t))(\mathrm{vy(}t)\cos\mathrm{(}q(t))-\mathrm{vx(}t)\sin\mathrm{(}q(t)))}{\mathrm{velb}}\\
% 	\frac{2(\mathrm{vy(}t)\cos\mathrm{(}q(t))-\mathrm{vx(}t)\sin\mathrm{(}q(t)))(\mathrm{vx(}t)\cos\mathrm{(}q(t))+\mathrm{vy(}t)\sin\mathrm{(}q(t)))}{\mathrm{vela}}+\\
% 	\frac{2(\mathrm{vy(}t)\cos\mathrm{(}q(t))-\mathrm{vx(}t)\sin\mathrm{(}q(t)))(\mathrm{vx(}t)(-\cos\mathrm{(}q(t)))-\mathrm{vy(}t)\sin\mathrm{(}q(t)))}{\mathrm{velb}}\\
% \end{array} \right\} 
            xgrad = xgrad + ...
                (2*cos(yaw)*vel_B(1)/ORIEN_VEL^2 - 2*sin(yaw)*vel_B(2)/VERDIT_VEL^2) * vect2';
            ygrad = ygrad + ...
                (2*sin(yaw)*vel_B(1)/ORIEN_VEL^2 + 2*cos(yaw)*vel_B(2)/VERDIT_VEL^2) * vect2';
            qgrad = qgrad + ...
                (2*vel_B(1)*vel_B(2)/ORIEN_VEL^2 + 2*(-vel_B(1))*vel_B(2)/VERDIT_VEL^2) * vect1';
        else
%             cost = cost + 0;
        end
    end
    grad((idi-1)*dim*n_order + 1:((idi-1)*dim+1)*n_order)     = xgrad;
    grad(((idi-1)*dim+1)*n_order + 1:((idi-1)*dim+2)*n_order) = ygrad;
    grad(((idi-1)*dim+2)*n_order + 1:((idi-1)*dim+3)*n_order) = qgrad;
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