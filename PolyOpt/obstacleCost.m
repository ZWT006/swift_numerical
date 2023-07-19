function [cost,grad] = obstacleCost(coeffs,segpoly)
%OBSTACLECOST 此处显示有关此函数的摘要
%   计算ESDF的cost
% 根据ds离散计算obstacle cost
% ObstacleCost = (dist-d_th)^2 * s
sdfmap = segpoly.sdf;
d_th = segpoly.d_th;
n_seg = segpoly.seg;
dim = segpoly.Dim;
n_order = segpoly.norder;

%##########################################################################
TimeOptimal = segpoly.TimeOptimal;
% 等式约束降维选项
ReduceOptimalValue = segpoly.ReduceOptimalValue;
%%%%%%%%%%%%%%%%%%% //Debug: 保存数据
if isfield(segpoly,'SAVE_DATA')
    SAVE_DATA = segpoly.SAVE_DATA;
else
    SAVE_DATA = false;
end
gdobs = [];

if (TimeOptimal)
    ts = coeffs(end-n_seg + 1:end); % 最后的n_seg个变量是tau
    coeffs(end-n_seg + 1:end)=[];
    ts = exp(ts);
else
    ts = segpoly.T;
end

if (ReduceOptimalValue)
    segpoly.T = ts;
    coeffs = segpoly.traj.SolveCoeffs(coeffs,segpoly);
end

polytraj = segpoly.traj.setCoeffs(coeffs);
% 此时的ts已经是变换后的时间了
polytraj = polytraj.setTarray(ts);

cost = 0;
grad = zeros(segpoly.coeffl,1);
gradt = zeros(n_seg,1); % 时间相关的梯度
PLOT_OBS = false;

for idi = 1:n_seg
    [pos,vel,acc] = polytraj.getStates(idi);
    Tdt = polytraj.bardt(idi);
    xgrad = zeros(segpoly.norder,1);
    ygrad = zeros(segpoly.norder,1);
    xgradt = 0; ygradt = 0;
    veccost = 0;
    for idj = 0:polytraj.bars(idi) % 比bars的长度+1
        dotpos = pos(idj+1,:);
        dotvel = vel(idj+1,:);
        dotacc = acc(idj+1,:);
        [dotdist,dotxgrad,dotygrad] = sdfmap.getDistAndGrad(dotpos(1),dotpos(2));
%         fprintf("dotxgrad=%2.4f; dotygrad=%2.4f\n",dotxgrad,dotygrad);
        if (dotdist < d_th ) % 超过阈值才计算cost %|| dotdist > 99
            % 使用倒数计算使曲线光滑
            veccost = veccost + 1/(2*(dotdist))^2 * (dotvel(1)^2+dotvel(2)^2) * Tdt;
            tvec = polytraj.getTvec(idj*Tdt);
            xgrad = xgrad + dotxgrad * tvec' * (-2*(dotdist)^3);
            ygrad = ygrad + dotygrad * tvec' * (-2*(dotdist)^3);
            % 感觉sdf对时间的偏导数是0呀
%             if (TimeOptimal)
%                 xgradt = xgradt + dotxgrad * vel(1) * idj*Tdt * (-2*(dotdist+0.01)^3);
%                 ygradt = ygradt + dotygrad * vel(2) * idj*Tdt * (-2*(dotdist+0.01)^3);
%             end
            % 这种计算二次距离差的方法造成cost曲线不够光滑
%             veccost = veccost + (dotdist - d_th)^2 * (dotvel(1)^2+dotvel(2)^2) * Tdt;
%             xgrad = xgrad + dotxgrad * tvec';
%             ygrad = ygrad + dotygrad * tvec';
%             if (TimeOptimal)
%                 xgradt = xgradt + dotxgrad * vel(1) * idj*Tdt;
%                 ygradt = ygradt + dotygrad * vel(2) * idj*Tdt;
%             end
        end
        if PLOT_OBS
            if (dotdist > 99)
                fprintf("obstacleCost in obstacle\n");
            end
        end
        %%%% //[xpos ypos dist xgrad ygrad Tdt cost]
        gdobs = [gdobs;dotpos(1),dotpos(2),dotdist,dotxgrad,dotygrad,idj*Tdt, 1/(2*(dotdist+0.01))^2 * (dotvel(1)^2+dotvel(2)^2) * Tdt];
    end
    cost = cost + veccost;
    grad((idi-1)*dim*n_order + 1:((idi-1)*dim+1)*n_order) = xgrad;
    grad(((idi-1)*dim+1)*n_order + 1:((idi-1)*dim+2)*n_order) = ygrad;
    % gradt 直接置零
    gradt(idi) = xgradt + ygradt;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% save debug data
if (SAVE_DATA)
    filename = "E:\datas\Swift\Debug\MATLABgdobs.csv";
    TRAJ_DATA = gdobs;
    TRAJ_DATA = round(TRAJ_DATA,6);
    writematrix(TRAJ_DATA,filename);
    fprintf("save gdobs debug datas, 上天保佑 \n");
end


if(ReduceOptimalValue)
    grad = segpoly.traj.getReduceOptVelue(grad);
end
if(TimeOptimal)
    grad = [grad;gradt];
end
cost = double(cost);
end


function grad = getGradatT(cp,c,t)
    cx = c(1:8);cy=c(9:16);
    cx(1)=[]; cy(1)=[];
grad = cp*(cx(1)^2 + cy(1)^2) + 2*cp*(2*cx(1)*cx(2) + 2*cy(1)*cy(2))*t + ...
   3*cp*((2/3)*(2*cx(2)^2 + 3*cx(1)*cx(3)) + (2/3)*(2*cy(2)^2 + 3*cy(1)*cy(3)))*t^2 + ...
   4*cp*(3*cx(2)*cx(3) + 2*cx(1)*cx(4) + 3*cy(2)*cy(3) + 2*cy(1)*cy(4))*t^3 + ...
   5*cp*((1/5)*(9*cx(3)^2 + 16*cx(2)*cx(4) + 10*cx(1)*cx(5)) + (1/5)*(9*cy(3)^2 + 16*cy(2)*cy(4) + 10*cy(1)*cy(5)))*t^4 + ...
   6*cp*((2/3)*(6*cx(3)*cx(4) + 5*cx(2)*cx(5) + 3*cx(1)*cx(6)) + (2/3)*(6*cy(3)*cy(4) + 5*cy(2)*cy(5) + 3*cy(1)*cy(6)))*t^5 + ...
   7*cp*((2/7)*(8*cx(4)^2 + 15*cx(3)*cx(5) + 12*cx(2)*cx(6) + 7*cx(1)*cx(7)) + (2/7)*(8*cy(4)^2 + 15*cy(3)*cy(5) + 12*cy(2)*cy(6) + 7*cy(1)*cy(7)))*t^6 + ...
   8*cp*((1/2)*(10*cx(4)*cx(5) + 9*cx(3)*cx(6) + 7*cx(2)*cx(7)) + (1/2)*(10*cy(4)*cy(5) + 9*cy(3)*cy(6) + 7*cy(2)*cy(7)))*t^7 + ...
   9*cp*((1/9)*(25*cx(5)^2 + 48*cx(4)*cx(6) + 42*cx(3)*cx(7)) + (1/9)*(25*cy(5)^2 + 48*cy(4)*cy(6) + 42*cy(3)*cy(7)))*t^8 + ...
   10*cp*((2/5)*(15*cx(5)*cx(6) + 14*cx(4)*cx(7)) + (2/5)*(15*cy(5)*cy(6) + 14*cy(4)*cy(7)))*t^9 + ...
   11*cp*((2/11)*(18*cx(6)^2 + 35*cx(5)*cx(7)) + (2/11)*(18*cy(6)^2 + 35*cy(5)*cy(7)))*t^10 + ...
   12*cp*(7*cx(6)*cx(7) + 7*cy(6)*cy(7))*t^11 + 13*cp*((49*cx(7)^2)/13 + (49*cy(7)^2)/13)*t^12;
end
