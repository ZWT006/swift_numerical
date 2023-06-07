function [cost,grad] = dynamicCost(coeffs,segpoly)
%DYNAMICCOST 此处显示有关此函数的摘要
%   此处显示详细说明

n_seg = segpoly.seg;
dim = segpoly.Dim;
n_order = segpoly.norder;

%##########################################################################
% 时间优化选项
TimeOptimal = segpoly.TimeOptimal;
% 等式约束降维选项
ReduceOptimalValue = segpoly.ReduceOptimalValue;

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
gradt = zeros(n_seg,1);

dynamic_rate = segpoly.dyna_rate;
pv_max = dynamic_rate * segpoly.pv_max;
pa_max = dynamic_rate * segpoly.pa_max;
wv_max = dynamic_rate * segpoly.wv_max;
wa_max = dynamic_rate * segpoly.wa_max;

velth = [pv_max,pv_max,wv_max];
accth = [pa_max,pa_max,wa_max];

for idi = 1:n_seg
    [pos,vel,acc] = polytraj.getStates(idi);
    jer = polytraj.getJerk(idi);
    Tdt = polytraj.bardt(n_seg);
    gradzeros = zeros(segpoly.norder,1);
    xvgrad = gradzeros;
    yvgrad = gradzeros;
    qvgrad = gradzeros;
    xagrad = gradzeros;
    yagrad = gradzeros;
    qagrad = gradzeros;
    temxvgrad = gradzeros;
    temyvgrad = gradzeros;
    temqvgrad = gradzeros;
    temxagrad = gradzeros;
    temyagrad = gradzeros;
    temqagrad = gradzeros;
    velcost = 0;acccost = 0;
    vgradt  = 0;agradt  = 0;
    for idj = 0:polytraj.bars(idi) % 比bars的长度+1
        % 计算梯度时不能取绝对值
        dotvel = vel(idj+1,:);
        dotacc = acc(idj+1,:);
        dotjer = jer(idj+1,:);
        delvel = [0,0,0];
        delacc = [0,0,0];
        Mt = getCoeffCons(idi*Tdt,n_order,4);
        vect2 = Mt(2,:);
        vect3 = Mt(3,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_{\mathrm{vel}}=\left( v-v_{\mathrm{th}} \right) ^2
% \\
% \nabla _{\mathrm{c}}f_{\mathrm{vel}}=2\cdot \left( v-v_{\mathrm{th}} \
% \right) \cdot \nabla _{\mathrm{c}}v=2\cdot \left( v-v_{\mathrm{th}} \
% \right) \cdot \frac{\partial v}{\partial c}
% \\
% \nabla _{\mathrm{t}}f_{\mathrm{vel}}=2\cdot \left( v-v_{\mathrm{th}} \
% \right) \cdot \n abla _{\mathrm{t}}v=2\cdot \left( v-v_{\mathrm{th}} \
% \right) \cdot a
% \\
% f_{\mathrm{acc}}=\left( a-a_{\mathrm{th}} \right) ^2
% \\
% \nabla _{\mathrm{c}}f_{\mathrm{acc}}=2\cdot \left( a-a_{\mathrm{th}} \
% \right) \cdot \n abla _{\mathrm{c}}a=2\cdot \left( a-a_{\mathrm{th}} \
% \right) \cdot \frac{\partial a}{\partial c}
% \\
% \nabla _{\mathrm{t}}f_{\mathrm{acc}}=2\cdot \left( a-a_{\mathrm{th}} \
% \right) \cdot \n abla _{\mathrm{t}}a=2\cdot \left( a-a_{\mathrm{th}} \
% \right) \cdot j,j \rightarrow jerk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 计算cost和gradient
        if (dotvel(1) > pv_max) 
            temxvgrad = 2 * (dotvel(1) - pv_max)* vect2';
            delvel(1) = dotvel(1) - pv_max;
        elseif (dotvel(1) < -pv_max) 
            temxvgrad = 2 * (dotvel(1) + pv_max)* vect2';
            delvel(1) = dotvel(1) + pv_max;
        else
            temxvgrad = gradzeros;
            delvel(1) = 0;
        end
        if (dotvel(2) > pv_max) 
            temyvgrad = 2 * (dotvel(2) - pv_max)* vect2';
            delvel(2) = dotvel(2) - pv_max;
        elseif (dotvel(2) < -pv_max)
            temyvgrad = 2 * (dotvel(2) + pv_max)* vect2';
            delvel(2) = dotvel(2) + pv_max;
        else
            temyvgrad = gradzeros;
            delvel(2) = 0;
        end
        if (dotvel(3) > wv_max) 
            temqvgrad = 2 * (dotvel(3) - wv_max)* vect2';
            delvel(3) = dotvel(3) - wv_max;
        elseif (dotvel(3) < -wv_max)
            temqvgrad = 2 * (dotvel(3) + wv_max)* vect2';
            delvel(3) = dotvel(3) + wv_max;
        else
            temqvgrad = gradzeros;
            delvel(3) = 0;
        end
        if (dotacc(1) > pa_max) 
            temxagrad = 2 * (dotacc(1) - pa_max)* vect3';
            delacc(1) = dotacc(1) - pa_max;
        elseif (dotacc(1) < -pa_max)
            temxagrad = 2 * (dotacc(1) + pa_max)* vect3';
            delacc(1) = dotacc(1) + pa_max;
        else
            temxagrad = gradzeros;
            delacc(1) = 0;
        end
        if (dotacc(2) > pa_max) 
            temyagrad = 2 * (dotacc(2) - pa_max)* vect3';
            delacc(2) = dotacc(2) - pa_max;
        elseif (dotacc(2) < -pa_max)
            temyagrad = 2 * (dotacc(2) + pa_max)* vect3';
            delacc(2) = dotacc(2) + pa_max;
        else
            temyagrad = gradzeros;
            delacc(2) = 0;
        end
        if (dotacc(3) > wa_max) 
            temqagrad = 2 * (dotacc(3) - wa_max)* vect3';
            delacc(3) = dotacc(3) - wa_max;
        elseif (dotacc(3) < -wa_max)
            temqagrad = 2 * (dotacc(3) + wa_max)* vect3';
            delacc(3) = dotacc(3) + wa_max;
        else
            temqagrad = gradzeros;
            delacc(3) = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % temxvgrad = 2 * dotvel(1) * vect2';
        % temyvgrad = 2 * dotvel(2) * vect2';
        % temqvgrad = 2 * dotvel(3) * vect2';
        % temxagrad = 2 * dotacc(1) * vect3';
        % temyagrad = 2 * dotacc(2) * vect3';
        % temqagrad = 2 * dotacc(3) * vect3';
        % dotvel = abs(vel(idj+1,:));
        % dotacc = abs(acc(idj+1,:));
        % detvel = dotvel-velth;
        % detacc = dotacc-accth;
        % if(dotvel(1) < pv_max) detvel(1)=0; temxvgrad=gradzeros; end
        % if(dotvel(2) < pv_max) detvel(2)=0; temxvgrad=gradzeros; end
        % if(dotvel(3) < wv_max) detvel(3)=0; temxvgrad=gradzeros; end
        % if(dotacc(1) < pa_max) detacc(1)=0; temxvgrad=gradzeros; end
        % if(dotacc(2) < pa_max) detacc(2)=0; temxvgrad=gradzeros; end
        % if(dotacc(3) < wa_max) detacc(3)=0; temxvgrad=gradzeros; end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 计算cost时取绝对值
        velcost = velcost + sum(delvel.^2);
        acccost = acccost + sum(delacc.^2);
        xvgrad = xvgrad + temxvgrad;
        yvgrad = yvgrad + temyvgrad;
        qvgrad = qvgrad + temqvgrad;
        xagrad = xagrad + temxagrad;
        yagrad = yagrad + temyagrad;
        qagrad = qagrad + temqagrad;
        % 计算时间梯度
        if(TimeOptimal)
            vgradt = vgradt + 2 * delvel * dotacc' * idi * Tdt;
            agradt = agradt + 2 * delacc * dotjer' * idi * Tdt; % 自然数对数的导数为自身
        end
    end % grad的计算
    cost = cost + velcost;% + acccost;
    % 这个梯度有问题呀,添加后就报错：Converged to an infeasible point
    grad((idi-1)*dim*n_order + 1:((idi-1)*dim+1)*n_order)     = xvgrad;% + xagrad;
    grad(((idi-1)*dim+1)*n_order + 1:((idi-1)*dim+2)*n_order) = yvgrad;% + yagrad;
    grad(((idi-1)*dim+2)*n_order + 1:((idi-1)*dim+3)*n_order) = qvgrad;% + qagrad;
    % gradt 直接置零
    if (TimeOptimal)
        gradt(idi) = vgradt + agradt;
    end
end
if(ReduceOptimalValue)
    grad = segpoly.traj.getReduceOptVelue(grad);
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
