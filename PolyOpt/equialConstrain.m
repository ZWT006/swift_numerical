function [ceq,grad] = equialConstrain(coeffs,segpoly)
%EQUIALCONSTRAIN Trajectory Optimization Equail Constrain
%   此处显示详细说明
n_seg = segpoly.seg;
dim = segpoly.Dim;
n_order = segpoly.norder;
%##########################################################################
% Aeq*x=beq nonlinear equality constraints 
TimeOptimal = segpoly.TimeOptimal;

if (TimeOptimal)
    ts = coeffs(end-n_seg + 1:end); % 最后的n_seg个变量是tau
    coeffs(end-n_seg + 1:end)=[];
    polytraj = segpoly.polytraj.settauarray(ts);
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


grad = zeros(segpoly.coeffl,1);
gradt = zeros(n_seg,1);

[Aeq, beq] = getAbeqMatrix([],segpoly);
ceq = Aeq*coeffs - beq;


end

