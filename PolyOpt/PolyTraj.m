classdef PolyTraj
    %POLYTRAJ 多项式轨迹优化的类
    %   此处显示详细说明
    properties
        n_seg; %piecen_segum
        Dim; % 优化变量的维度
        nodr; % 多项式的阶次 固定为8？
        ncost; % mini 的阶次 固定为4？
        ninput; % 输入的阶次 默认与ncost一致
        start_cond;
        goal_cond;
        T; %Time Vector
        tau; % T = exp(tau)
        R; % theta的权重系数
        Aeq; % Ac=b
        beq; % 
        coeffs; % 
        JgdC;
%         gdP;
        JgdT;
        Tpow;
        ds; % 距离分辨率
        dt; % 时间分辨率(用于优化问题计算) 而不是多项式的计算精度
        dists; % [x,y]的欧氏距离
        d_th; % 距离代价的阈值
        bars;   % 每段seg根据ds分bars
        dotl; % 离散点的总个数
        bardt; % 每段seg的bars的dt
        coeffl; % 系数变量的个数
        dof; % 优化问题系数自由度
    end
    methods
        %% 功能函数区 init set 
        function obj = PolyTraj(segpoly)
            %POLYTRAJ 
            %   矩阵初始化
            obj.n_seg = segpoly.seg;
            obj.Dim = segpoly.Dim;
            obj.nodr = segpoly.norder;
            obj.ncost = segpoly.ncost;
            obj.ninput = segpoly.ninput;
            obj.T = segpoly.T;
            obj.R = segpoly.R;
            obj.ds = segpoly.ds;
            obj.dt = segpoly.dt;
            obj.dists = segpoly.dists;
            obj.bars = round(obj.dists./obj.ds);
            obj.dotl = sum(obj.bars) + obj.n_seg;
            obj.coeffl = segpoly.coeffl;
            obj.d_th = segpoly.d_th;
            obj.dof  = segpoly.dof;
        end
        function obj = setCoeffs(obj,coeffs)
            % 设置[x,y,q]的系数
            obj.coeffs = coeffs(:);
        end
        function obj = setTarray(obj,T)
            % 设置每段seg的时间 T = t
            obj.T = T;
            obj.bardt = T./obj.bars;
        end
        function obj = settauarray(obj,tau)
            % 设置每段seg的时间 T = et
            obj.tau = tau;
            obj.T = exp(tau);
            obj.bardt = obj.T./obj.bars;
        end
        %% 辅助函数区
        function showSegState(obj)
            for seg = 1:obj.n_seg
                fprintf("seg = %2d, T = %2.4f, bart = %2.4f, dist = %2.4f bars = %2d\n",seg,obj.T(seg),obj.bardt(seg),obj.dists(seg),obj.bars(seg));
            end
        end
        function Tvec = tau2T(tau)
            Tvec = exp(tau);
        end
        function tauvec = T2tau(T)
            tauvec = log(T);
        end
        function tvec = getTvec(~,t)
            tvec = getCoeffCons(t,8,1);
        end
        function tvec = gettauvec(~,tau)
            t = exp(tau);
            tvec = getCoeffCons(t,8,1);
        end
        %% 功能函数区 get
        function [cost,grad] = getsmoothcost(obj,coeffes)
            % 计算smoothcost 以及对应的grad
%             cost = 0;
%             for j = 1: obj.n_seg
%                 % 取出 nodr长度的系数
%                 xc = c((j-1)*obj.nodr*obj.Dim + 1:(j-1)*obj.nodr*obj.Dim + obj.nodr);
%                 xcost = xc'*getSegPolyCostMatrix(obj.nodr,obj.ncost,obj.T(j))*xc;
%                 yc = c((j-1)*obj.nodr*obj.Dim + obj.nodr + 1:(j-1)*obj.nodr*obj.Dim + obj.nodr*2);
%                 ycost = yc'*getSegPolyCostMatrix(obj.nodr,obj.ncost,obj.T(j))*yc;
%                 qc = c((j-1)*obj.nodr*obj.Dim + obj.nodr*2 + 1:(j-1)*obj.nodr*obj.Dim + obj.nodr*3);
%                 qcost = qc'*getSegPolyCostMatrix(obj.nodr,obj.ncost,obj.T(j))*qc;
%                 cost = xcost + ycost + obj.R * qcost;
% %                 fprintf("Matrix smooth xcost :%f \n",xcost);
% %                 fprintf("Matrix smooth ycost :%f \n",ycost);
% %                 fprintf("Matrix smooth qcost :%f \n",qcost);
% %                 fprintf("<============================>\n");
%             end
%             fprintf("Matrix smooth cost :%f \n",cost);
%             fprintf("<============================>\n");
            cost = 0;
            for j = 1: obj.n_seg
                % 取出 nodr长度的系数
                xc = coeffes((j-1)*obj.nodr*obj.Dim + 1:(j-1)*obj.nodr*obj.Dim + obj.nodr);
                xcost = getSegPolyCost(xc,obj.T(j));
                yc = coeffes((j-1)*obj.nodr*obj.Dim + obj.nodr + 1:(j-1)*obj.nodr*obj.Dim + obj.nodr*2);
                ycost = getSegPolyCost(yc,obj.T(j));
                qc = coeffes((j-1)*obj.nodr*obj.Dim + obj.nodr*2 + 1:(j-1)*obj.nodr*obj.Dim + obj.nodr*3);
                qcost = getSegPolyCost(qc,obj.T(j));
                cost = cost + xcost + ycost + obj.R * qcost;
%                 fprintf("Equation smooth xcost :%f \n",xcost);
%                 fprintf("Equation smooth ycost :%f \n",ycost);
%                 fprintf("Equation smooth qcost :%f \n",qcost);
%                 fprintf("<============================>\n");
            end
%             fprintf("Equation smooth cost :%f \n",cost);
%             fprintf("<============================>\n");
            grad = zeros(1,obj.nodr*obj.Dim);
            gradt = zeros(1,obj.n_seg);
            for j = 1: obj.n_seg
                % 取出 nodr长度的系数
                xc = coeffes((j-1)*obj.nodr*obj.Dim + 1:(j-1)*obj.nodr*obj.Dim + obj.nodr);
                xgrad = getSegPolyGardatc(xc,obj.T(j));
                yc = coeffes((j-1)*obj.nodr*obj.Dim + obj.nodr + 1:(j-1)*obj.nodr*obj.Dim + obj.nodr*2);
                ygrad = getSegPolyGardatc(yc,obj.T(j));
                qc = coeffes((j-1)*obj.nodr*obj.Dim + obj.nodr*2 + 1:(j-1)*obj.nodr*obj.Dim + obj.nodr*3);
                qgrad = getSegPolyGardatc(qc,obj.T(j));
                grad = [grad ,xgrad,ygrad,obj.R * qgrad];
                xgradt = getSegPolyGardatT(xc,obj.T(j));
                ygradt = getSegPolyGardatT(yc,obj.T(j));
                qgradt = obj.R * getSegPolyGardatT(qc,obj.T(j));
                gradt(j) = (xgradt + ygradt + qgradt)*obj.T(j); % 对数的导数还是自身
            end
            grad=grad';
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [pos,vel,acc]= getStates(obj,seg)
            % 返回[x,y,q]的pos矩阵,离散程度用 dist/ds
            if (nargin == 1) % 没有选择seg 就返回所有seg的pos
                pos = zeros(obj.dotl,3);
                vel = pos; acc = pos; 
                for seg = 1:obj.n_seg
                    Tts = linspace(0,obj.T(seg),obj.bars(seg)+1);
%                     fprintf("seg = %2d, T = %2.4f, dist = %2.4f bars = %2d\n",seg,obj.T(seg),obj.dists(seg),obj.bars(seg));
                    xc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr);
                    yc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + obj.nodr + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr*2);
                    qc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + obj.nodr*2 + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr*3);
                    states = getstatexyq(xc,yc,qc,Tts);
                    idg = sum(obj.bars(1:seg)) + seg;
                    ids = idg - obj.bars(seg);
                    pos(ids:idg,:) = squeeze(states(1,:,:))';
                    vel(ids:idg,:) = squeeze(states(2,:,:))';
                    acc(ids:idg,:) = squeeze(states(3,:,:))';
                end
            elseif(nargin == 2) % 选择seg 只返回对应seg的pos
                if (seg > obj.n_seg)
                    fprintf("getPos seg = %d, over range%d",seg,obj.n_seg);
                    return;
                end
                Tts = linspace(0,obj.T(seg),obj.bars(seg)+1);
                xc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr);
                yc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + obj.nodr + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr*2);
                qc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + obj.nodr*2 + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr*3);
                states = getstatexyq(xc,yc,qc,Tts);
                pos = squeeze(states(1,:,:))';
                vel = squeeze(states(2,:,:))';
                acc = squeeze(states(3,:,:))';
            end
%             fprintf("getPos inputs = %d\n",nargin);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function jer = getJerk(obj,seg)
            % 返回[x,y,q]的jer,离散程度用 dist/ds
            if (nargin == 1) % 没有选择seg 就返回所有seg的pos
                jer = zeros(obj.dotl,3);
                for seg = 1:obj.n_seg
                    Tts = linspace(0,obj.T(seg),obj.bars(seg)+1);
%                     fprintf("seg = %2d, T = %2.4f, dist = %2.4f bars = %2d\n",seg,obj.T(seg),obj.dists(seg),obj.bars(seg));
                    xc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr);
                    yc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + obj.nodr + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr*2);
                    qc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + obj.nodr*2 + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr*3);
                    states = getJerk(xc,yc,qc,Tts);
                    idg = sum(obj.bars(1:seg)) + seg;
                    ids = idg - obj.bars(seg);
                    jer(ids:idg,:) = states;
                end
            elseif(nargin == 2) % 选择seg 只返回对应seg的pos
                if (seg > obj.n_seg)
                    fprintf("getPos seg = %d, over range%d",seg,obj.n_seg);
                    return;
                end
                Tts = linspace(0,obj.T(seg),obj.bars(seg)+1);
                xc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr);
                yc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + obj.nodr + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr*2);
                qc = obj.coeffs((seg-1)*obj.nodr*obj.Dim + obj.nodr*2 + 1:(seg-1)*obj.nodr*obj.Dim + obj.nodr*3);
                jer = getJerk(xc,yc,qc,Tts);
            end
%             fprintf("getPos inputs = %d\n",nargin);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function coeffs = SolveCoeffs(obj,x,segpoly)
            [xl,~] = size(x);
            [Aeq, beq] = getAbeqMatrix([],segpoly);
            dof_num = 0;
            for id_seg=0:segpoly.seg-1
                if (dof_num >= segpoly.dof)
                    break;
                end
                for id_dim = 0:segpoly.Dim-1
                    if (dof_num >= segpoly.dof)
                        break;
                    end
                    for id_order=segpoly.ninput+2:segpoly.norder
                        row_I = zeros(1,obj.coeffl); % row 向量
                        % 当前优化变量的行序号
                        row_num = id_seg*segpoly.Dim*segpoly.norder + id_dim*segpoly.norder + id_order;
                        row_I(row_num) = 1;
                        dof_num = dof_num + 1;
                        Aeq = insertRow(Aeq,row_I,row_num);
                        beq = insertRow(beq,x(dof_num),row_num);
                        if (dof_num >= segpoly.dof)
                            break;
                        end
                    end
                end
            end
%             coeffs = linsolve(Aeq,beq);
            coeffs = pinv(Aeq)*beq;
        end
        function [Aeq,beq] = SolveAeqbeq(obj,x,segpoly)
            [Aeq, beq] = getAbeqMatrix([],segpoly);
            dof_num = 0;
            for id_seg=0:segpoly.seg-1 
                if (dof_num >= segpoly.dof)
                    break;
                end % 5~8 阶次系数
                for id_dim = 0:segpoly.Dim-1
                    if (dof_num >= segpoly.dof)
                        break;
                    end
                    for id_order=segpoly.ninput+2:segpoly.norder
                        row_I = zeros(1,obj.coeffl); % row 向量
                        % 当前优化变量的行序号
                        row_num = id_seg*segpoly.Dim*segpoly.norder + id_dim*segpoly.norder + id_order;
                        row_I(row_num) = 1;
                        dof_num = dof_num + 1;
                        Aeq = insertRow(Aeq,row_I,row_num);
                        beq = insertRow(beq,x(dof_num),row_num);
                        if (dof_num >= segpoly.dof)
                            break;
                        end
                    end
                end
            end
        end
        function re_OptVelue = getReduceOptVelue(obj,OptVelue)
%             grad(obj.dof+1:end) = [];
            dof_num = 0;
            [row,col]=size(OptVelue);
            if (row > 1)
                row = obj.dof;
            end
            if (col > 1)
                col = obj.dof;
            end
            re_OptVelue = zeros(row,col);
            for id_seg=0:obj.n_seg-1 
                if (dof_num >= obj.dof)
                    break;
                end % 5~8 阶次系数
                for id_dim = 0:obj.Dim-1
                    if (dof_num >= obj.dof)
                        break;
                    end
                    for id_order=obj.ninput+2:obj.nodr
                        % 当前优化变量的行序号
                        row_num = id_seg*obj.Dim*obj.nodr + id_dim*obj.nodr + id_order;
                        dof_num = dof_num + 1;
                        re_OptVelue(dof_num) = OptVelue(row_num);
                        if (dof_num >= obj.dof)
                            break;
                        end
                    end
                end
            end
        end
    end
end

%% 计算等式约束Ac=b的 A矩阵 和 b矩阵
function [Aeq,beq]= getAbeq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond)
% Minimun Snap Trajectory Generation P35 L5.pdf
% 注意 Aeq_start;Aeq_end;Aeq_wp;Aeq_con 的列数相同为总优化变量个数n_seg*(n_order+1)
    n_all_poly = n_seg*(n_order+1);
    coeff_n = n_order+1;
    %#####################################################
    % p,v,a,j 的起点约束, 约束数量与n_inputorder有关
    Aeq_start = zeros(n_inputorder, n_all_poly);
    % STEP 2.1: write expression of Aeq_start and beq_start
    Aeq_start(1:n_inputorder,1:coeff_n) = getCoeffCons(0,n_order,n_inputorder);
    beq_start =  start_cond';
    
    %#####################################################
    % p,v,a,j 的终端约束
    Aeq_end = zeros(n_inputorder, n_all_poly);
    t = ts(end);
    Aeq_end(1:n_inputorder,(end-coeff_n+1):end) = getCoeffCons(t,n_order,n_inputorder);
    beq_end =  end_cond';
    
    %#####################################################
    % 中点的位置约束
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    
    for k = 0:1:n_seg-2
        beq_wp(k+1, 1) = waypoints(k+2);
        coeff = getCoeffCons(ts(k+1),n_order,1); % 中间点只固定 pose
        % 1:t:t^2:t^3:t^4…… * [p0 p1 p2 p3……]T
        Aeq_wp(k+1, 1+k*coeff_n:(1+k)*coeff_n) = coeff(1, :);  
    end
    
    %#####################################################
    % 连续性约束
    row = n_inputorder*(n_seg-1);
    col = n_all_poly;
    Aeq_con = zeros(row,col);
    beq_con = zeros(row,1);
    
    for k = 0:1:n_seg-2 
        Aeq_con(1+n_inputorder*k:n_inputorder*(k+1),1+coeff_n*k:coeff_n*(k+1)) = getCoeffCons(ts(k+1),n_order,n_inputorder);
        Aeq_con(1+n_inputorder*k:n_inputorder*(k+1),1+coeff_n*(k+1):coeff_n*(k+2)) = -getCoeffCons(0,n_order,n_inputorder);            
    end
    
    %#####################################################
    % 构造约束矩阵
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
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
% A   coeff = [1,  1*t,  1*t^2,  1*t^3,  1*t^4,  1*t^5,  1*t^6,  1*t^7;
%              0,  1,    2*t,    3*t^2,  4*t^3,  5*t^4,  6*t^5,  7*t^6;
%              0,  0,    2,      6*t,    12*t^2, 20*t^3, 30*t^4, 42*t^5;
%              0,  0,    0,      6,      24*t,   60*t^2, 120*t^3,210*t^4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B   =  左右翻转A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 计算SmoothCostd的
function Q = getSegPolyCostMatrix(n_order,n_costorder,ts)
Q = zeros(n_order, n_order );
for i = n_costorder:n_order-1
    for j = n_costorder:n_order-1
        % 根据离散的多项式cost计算minisnap应该是固定的(i+j-7) 因为是原多项式求4阶导
        % 另外真实的积分应该是有T(j)^(i+j-7)-T(j-1)^(i+j-7),这里的时间使用的是相对的所以后项应该是0
        % 如果 cost 的阶次由 n_costorder 设定,则Q矩阵大小可变 row = n_costorder * k;
        % row = (n_order + 1 - n_costorder) * k
        Q(i+1,j+1) = factorial(i)/factorial(i-n_costorder)*factorial(j)/factorial(j-n_costorder)/(i+j-n_costorder*2+1)*ts^(i+j-n_costorder*2+1);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = getSegPolyCost(c,t)
c(1)=[]; %MATLAB从1开始索引,所以c(7)应该是符号公式中的第8个系数
cost = 576 *c(4)^2*t + 2880*c(4)*c(5)*t^2 + 4800*c(5)^2*t^3 + 5760*c(4)*c(6)*t^3 + ... 
       21600*c(5)*c(6)*t^4 + 10080*c(4)*c(7)*t^4 + 25920*c(6)^2*t^5 + ... 
       40320*c(5)*c(7)*t^5 + 100800*c(6)*c(7)*t^6 + 100800*c(7)^2*t^7 ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gard = getSegPolyGardatc(c,t)
c(1)=[];
gard =  [0, 0, 0, 0, 1152*c(4)*t + 2880*c(5)*t^2 + 5760*c(6)*t^3 + 10080*c(7)*t^4,... 
         2880*c(4)*t^2 + 9600*c(5)*t^3 + 21600*c(6)*t^4 + 40320*c(7)*t^5, ...
         5760*c(4)*t^3 + 21600*c(5)*t^4 + 51840*c(6)*t^5 + 100800*c(7)*t^6, ...
         10080*c(4)*t^4 + 40320*c(5)*t^5 + 100800*c(6)*t^6 + 201600*c(7)*t^7];
end
function gard = getSegPolyGardatT(c,t)
c(1)=[];
gard =   576*c(4)^2 + 5760*c(4)*c(5)*t + 14400*c(5)^2*t^2 + 17280*c(4)*c(6)*t^2 + ...
         86400*c(5)*c(6)*t^3 + 40320*c(4)*c(7)*t^3 + 129600*c(6)^2*t^4 + ...
         201600*c(5)*c(7)*t^4 + 604800*c(6)*c(7)*t^5 + 705600*c(7)^2*t^6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 根据多项式系数和时间计算[pos,vel,acc]
function state = getstatexyq(xc,yc,qc,t)
% TODO: 这里是固定了阶次,之后可以优化代码
tl = length(t);
state(3,3,tl) = 0;
    for i=1:tl
        M = getCoeffCons(t(i),8,3);
        xstate = M*xc;
        ystate = M*yc;
        qstate = M*qc;
        state(:,:,i)=[xstate,ystate,qstate];
    end
end
function state = getJerk(xc,yc,qc,t)
tl = length(t);
state(tl,3) = 0;
    for i=1:tl
        M = getCoeffCons(t(i),8,4);
        xstate = M*xc;
        ystate = M*yc;
        qstate = M*qc;
        state(i,:)=[xstate(4),ystate(4),qstate(4)];
    end
end
%% 针对矩阵进行操作 插入行
function Matrix = insertRow(Matrix,row,idx)
    % insert row = idx
    Matrix =  [Matrix(1:(idx-1),:);row;Matrix(idx:end,:)];
end