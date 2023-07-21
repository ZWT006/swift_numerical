%***************************************
%@Author: Wentao Zhang
%@Date: 2023-3-1
%@E-mail: zwt190315@163.com
%@Reference: ##########
%@Problems: 
%@Description:
%@TODO：##########
%***************************************
%%%%%%%%
% Step1: 设置优化问题初始参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
% 加载waypoints
load("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\path.mat");
addpath("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\LocalOBVP\")
RATION = 100;
path(:,1)=path(:,1)/RATION;
path(:,2)=path(:,2)/RATION;
BEZIER_POLY = 1;
NORMAL_POLY = 2;
POLY_TYPE = NORMAL_POLY;

OP_structure.QP_inequality = false;

% for i=1:10
%     path(end,:)=[];
% end

% path(2,:)=[];

% path(2,:)=[];
% path(1,:)=[];
% path(end,:)=[];
% path = ginput() * 100.0;
% path = [50,	                50;
%         145.659906173533,	186.312819235911;
%         191.836813342538,	272.271495258284;
%         267.755189847367,	328.985619044511;
%         392.144396900521,	385.147806029386;
%         514.147782625023,	469.962638842209;
%         613.408785998573,	537.502971052556;
%         702.631862930347,	605.144253906707;
%         794.630236687311,	668.449239620839;
%         852.120251855336,	664.522777123541;
%         940,	            550];
% path(:,3)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% segment polynomial optimization coefficients
% n > i > j
% polynomial order  : n
% input order   : i (i > j means input control minimum value)
% cost order    : j (cost function order)
% segmant num   : k
% segmant dof   : k*(n+1)
% constraint num: j(start state)+j(final state)+(k-1)(mid points pose)
% unknow num    : k*(n+1)
% equal constraint: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_order       = 7;          % 多项式的阶数
n_costorder   = 4;          % 最小化的求导阶次 0=posi;1=vel;2=acc;3=jerk;4=snap;
n_inputorder  = 4;          % 输入的阶次 可以理解为segment 之间满足等式约束的阶次
n_seg         = size(path,1)-1; % 分段数
n_poly_perseg = (n_order+1);    % 每段的系数个数

%%%%%%%%
%Step2: 根据距离长度计算每两个点之间的时间分配
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vel_factor = 1.4; % reference Linear Velocity  2m/s
W_factor   = 1.4; % reference Angular Velocity rad/s
pv_max = Vel_factor*1.5;
pa_max = Vel_factor*1.5;
wv_max = W_factor*1.5;
wa_max = W_factor*1.5;
dist= zeros(n_seg, 1);
ts  = ones(n_seg, 1)*0.8;

% 计算参考线速度时间
for i = 1:n_seg
    dist(i) = sqrt((path(i+1, 1)-path(i, 1))^2 + (path(i+1, 2) - path(i, 2))^2);
    ts(i) = dist(i)/Vel_factor;
end
% 计算参考角速度时间/并选择较大值
for i = 1:n_seg
    t_temp = AngleDelta(path(i+1, 3),path(i,3))/W_factor;
    if (t_temp > ts(i))
        ts(i) = t_temp;
    end
end


ts(1)   = ts(1)*2;
ts(end) = ts(end)*2;
T = sum(ts);

%%%%%%%%%%%%%%%%
% 处理角度变化问题
path_q = path(:, 3);
path_deg=rad2deg(path_q);
path_m = zeros(n_seg+1,1);
path_m(1)=path_q(1);
for idx=2:n_seg+1
%     path_m(idx)=AngleDelta(path_q(idx-1),path_q(idx));
    path_m(idx)=path_m(idx-1)+AngleDelta(path_q(idx-1),path_q(idx));
end
path_deg_m=rad2deg(path_m);

path(:,3)=path_m;
%%%%%%%%
%Step3: 使用QP求解器求解多项式系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Clock #############################
tic;
switch POLY_TYPE
    case BEZIER_POLY

        bezier_coef_x = MinimumBesierSolver(path(:, 1), ts, n_seg, n_order, n_costorder, n_inputorder, pv_max, pa_max);
        bezier_coef_y = MinimumBesierSolver(path(:, 2), ts, n_seg, n_order, n_costorder, n_inputorder, pv_max, pa_max);
        bezier_coef_q = MinimumBesierSolver(path(:, 3), ts, n_seg, n_order, n_costorder, n_inputorder, wv_max, wa_max);
        M_c = BerneteinCoeff(n_order);
        power=0:1:n_order;
        for idx=0:n_seg-1
            t = ts(idx+1);
            scal_t=t.^power;
            scal_t=scal_t';
            xseg_coeff = bezier_coef_x((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
            poly_coef_x((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1)=M_c*xseg_coeff./scal_t;
            yseg_coeff = bezier_coef_y((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
            poly_coef_y((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1)=M_c*yseg_coeff./scal_t;
            qseg_coeff = bezier_coef_q((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
            poly_coef_q((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1)=M_c*qseg_coeff./scal_t;
        end

    case NORMAL_POLY
        OP_structure.v_max = pv_max;
        OP_structure.a_max = pa_max;
        poly_coef_x = MinimumPolySolver(path(:, 1), ts, n_seg, n_order, n_costorder, n_inputorder,OP_structure);
%         poly_coef_x = MinimumPolySolver(path(:, 1), ts, n_seg, n_order, n_costorder, n_inputorder);
        OP_structure.v_max = pv_max;
        OP_structure.a_max = pa_max;
        poly_coef_y = MinimumPolySolver(path(:, 2), ts, n_seg, n_order, n_costorder, n_inputorder,OP_structure);
%         poly_coef_y = MinimumPolySolver(path(:, 2), ts, n_seg, n_order, n_costorder, n_inputorder);
        OP_structure.v_max = wv_max;
        OP_structure.a_max = wa_max;
        poly_coef_q = MinimumPolySolver(path(:, 3), ts, n_seg, n_order, n_costorder, n_inputorder,OP_structure);
%         poly_coef_q = MinimumPolySolver(path(:, 3), ts, n_seg, n_order, n_costorder, n_inputorder);
end
datename = strcat(num2str(n_costorder),"-cost-",num2str(n_inputorder),"-input");
% filenema = strcat("opt_coeff_",num2str(n_costorder),"_cost_",num2str(n_inputorder),"_input");
% save(filenema);

%%
% 用于显示轨迹
X_n = [];
Y_n = [];
Q_n = [];
X_dn = [];
Y_dn = [];
Q_dn = [];
X_ddn = [];
Y_ddn = [];
Q_ddn = [];
bezier_coef_xd=[];
bezier_coef_yd=[];
bezier_coef_qd=[];
bezier_coef_xdd=[];
bezier_coef_ydd=[];
bezier_coef_qdd=[];

k = 1;
tstep = 0.01;
switch POLY_TYPE
    case NORMAL_POLY
        for idx=0:n_seg-1
            %#####################################################
            % STEP 3: get the coefficients of i-th segment of both x-axis
            % and y-axis
            Pxi = poly_coef_x((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
            Pyi = poly_coef_y((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
            Pqi = poly_coef_q((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
            for t = 0:tstep:ts(idx+1)
                px=flip(Pxi);
                py=flip(Pyi);
                pq=flip(Pqi);
                X_n(k)  = polyval(px, t)*RATION;
                Y_n(k)  = polyval(py, t)*RATION;
                Q_n(k)  = polyval(pq, t);
                % velocity
                pdx=polyder(px);
                pdy=polyder(py);
                pdq=polyder(pq);
                X_dn(k)  = polyval(pdx, t);
                Y_dn(k)  = polyval(pdy, t);
                Q_dn(k)  = polyval(pdq, t);
                % accelaration
                pddx=polyder(pdx);
                pddy=polyder(pdy);
                pddq=polyder(pdq);
                X_ddn(k)  = polyval(pddx, t);
                Y_ddn(k)  = polyval(pddy, t);
                Q_ddn(k)  = polyval(pddq, t);
                k = k + 1;
            end
            fprintf('t = %2.6f;ts(idx) =%2.6f \n',t,ts(idx+1))
        end

    case BEZIER_POLY
        for idx=0:n_seg-1
            %##############################################################
            %#用bezier曲线控制点计算多阶导数 p=位置
            M_c = BerneteinCoeff(n_order);
            power=0:1:n_order;
            t = ts(idx+1);
            scal_t=t.^power;
            scal_t=scal_t';
            xseg_coeff = bezier_coef_x((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
            px=M_c*xseg_coeff./scal_t;
            yseg_coeff = bezier_coef_y((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
            py=M_c*yseg_coeff./scal_t;
            qseg_coeff = bezier_coef_q((n_order+1)*(idx)+1:(n_order+1)*(idx)+n_order+1);
            pq=M_c*qseg_coeff./scal_t;
            %##############################################################
            %#用bezier曲线控制点计算多阶导数 v=速度
            M_c = BerneteinCoeff(n_order-1);
            power=0:1:n_order-1;
            t = ts(idx+1);
            scal_t=t.^power;
            scal_t=scal_t';
            xdseg_coeff = n_order*(xseg_coeff(2:end)-xseg_coeff(1:end-1))/t;
            bezier_coef_xd = [bezier_coef_xd ;xdseg_coeff];
            pdx=M_c*xdseg_coeff./scal_t;
            ydseg_coeff = n_order*(yseg_coeff(2:end)-yseg_coeff(1:end-1))/t;
            bezier_coef_yd = [bezier_coef_yd ;ydseg_coeff];
            pdy=M_c*ydseg_coeff./scal_t;
            qdseg_coeff = n_order*(qseg_coeff(2:end)-qseg_coeff(1:end-1))/t;
            bezier_coef_qd = [bezier_coef_qd ;qdseg_coeff];
            pdq=M_c*qdseg_coeff./scal_t;
            %##############################################################
            %#用bezier曲线控制点计算多阶导数 a=加速度
            M_c = BerneteinCoeff(n_order-2);
            power=0:1:n_order-2;
            t = ts(idx+1);
            scal_t=t.^power;
            scal_t=scal_t';
            xddseg_coeff = (n_order-1)*(xdseg_coeff(2:end)-xdseg_coeff(1:end-1))/t;
            bezier_coef_xdd = [bezier_coef_xdd ;xddseg_coeff];
            pddx=M_c*xddseg_coeff./scal_t;
            yddseg_coeff = (n_order-1)*(ydseg_coeff(2:end)-ydseg_coeff(1:end-1))/t;
            bezier_coef_ydd = [bezier_coef_ydd ;yddseg_coeff];
            pddy=M_c*yddseg_coeff./scal_t;
            qddseg_coeff = (n_order-1)*(qdseg_coeff(2:end)-qdseg_coeff(1:end-1))/t;
            bezier_coef_qdd = [bezier_coef_qdd ;qddseg_coeff];
            pddq=M_c*qddseg_coeff./scal_t;
            px=flip(px);
            py=flip(py);
            pq=flip(pq);

            pdx=flip(pdx);
            pdy=flip(pdy);
            pdq=flip(pdq);

            pddx=flip(pddx);
            pddy=flip(pddy);
            pddq=flip(pddq);
            for t = 0:tstep:ts(idx+1)
%             for t = 0:tstep:1
                
                X_n(k)  = polyval(px, t)*RATION;
                Y_n(k)  = polyval(py, t)*RATION;
                Q_n(k)  = polyval(pq, t);
                % velocity
                
                X_dn(k)  = polyval(pdx, t);
                Y_dn(k)  = polyval(pdy, t);
                Q_dn(k)  = polyval(pdq, t);
                % accelaration
                
                X_ddn(k)  = polyval(pddx, t);
                Y_ddn(k)  = polyval(pddy, t);
                Q_ddn(k)  = polyval(pddq, t);
                k = k + 1;
            end
        end
end

k=k-1;

%%%%%%%%
%Step4: 绘制结果图例
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%<------------------------------------------------------------------->%%%
n=0;
% xy map & palstance
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on

pic_title = strcat("xy map & palstance",datename);
title(pic_title);
AngUNIT=50;
[path_length,~] = size(path);
for idx=1:path_length
%     plot([path(idx,1),path(idx,1)+cos(path(idx,3))*AngUNIT],[path(idx,2),path(idx,2)+sin(path(idx,3))*AngUNIT],'-','Color','b','LineWidth',1);
    plot([path(idx,1)*RATION,path(idx,1)*RATION+cos(path(idx,3))*AngUNIT],[path(idx,2)*RATION,path(idx,2)*RATION+sin(path(idx,3))*AngUNIT],'-','Color','b','LineWidth',1);
%     plot([path(idx,1)*RATION,path(idx,1)*RATION+cos(path_m(idx))*AngUNIT],[path(idx,2)*RATION,path(idx,2)*RATION+sin(path_m(idx))*AngUNIT],'c-','LineWidth',2);
end
for idx=1:path_length-1
    plot([path(idx,1)*RATION,path(idx+1,1)*RATION],[path(idx,2)*RATION,path(idx+1,2)*RATION],'g-','LineWidth',1);
end
AngUNIT=1;
for idx=1:k
    detl = Q_dn(idx)*AngUNIT+2;
    theta = Q_n(idx);
    plot([X_n(idx),X_n(idx)+cos(theta)*detl],[Y_n(idx),Y_n(idx)+sin(theta)*detl],'c-','LineWidth',1);
end
if (POLY_TYPE == BEZIER_POLY)
    for idx=1:n_seg
        if(mod(idx,2))
            scatter(bezier_coef_x((idx-1)*(n_order+1)+1:idx*(n_order+1))*RATION,bezier_coef_y((idx-1)*(n_order+1)+1:idx*(n_order+1))*RATION,100,"g");
        else
            scatter(bezier_coef_x((idx-1)*(n_order+1)+1:idx*(n_order+1))*RATION,bezier_coef_y((idx-1)*(n_order+1)+1:idx*(n_order+1))*RATION,100,"c");
        end
    end
end
plot(X_n, Y_n , 'r-');
% scatter(path(1:size(path, 1), 1), path(1:size(path, 1), 2),'ok');
scatter(path(1:size(path, 1), 1)*RATION, path(1:size(path, 1), 2)*RATION,'ok');
grid on

hold off
axis equal
xlim([0 1200]);
ylim([0 800]);
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% xy velocity
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("xy velocity",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
plot(tv,X_dn, 'r-');
plot(tv,Y_dn, 'b-');
t_temp=0;
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp,X_dn(qdn_idx),'*r');
end
t_temp=0;
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp,Y_dn(qdn_idx),'*b');
end
% legend('x vel','y vel');
grid on
xlim([0 ceil(T)]);
% ylim([-150 250]);
ylim([floor(min(Y_dn))-2 ceil(max(Y_dn))+2]);
line([0,ceil(T)],[pv_max,pv_max],'linestyle','--','color','c','LineWidth',2);
line([0,ceil(T)],[-pv_max,-pv_max],'linestyle','--','color','c','LineWidth',2);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% angular velocity
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("angular velocity",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
plot(tv,Q_dn, 'b-');
t_temp = 0;
scatter(t_temp,Q_dn(1),'*r');
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp,Q_dn(qdn_idx),'*r');
end
grid on
xlim([0 ceil(T)]);
ylim([floor(min(Q_dn)) ceil(max(Q_dn))]);
line([0,ceil(T)],[wv_max,wv_max],'linestyle','--','color','c','LineWidth',2);
line([0,ceil(T)],[-wv_max,-wv_max],'linestyle','--','color','c','LineWidth',2);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% angular position
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("angular position",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
plot(tv,rad2deg(Q_n), 'b-');
t_temp = 0;
scatter(t_temp,rad2deg(path_m(1)),'*r');
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    scatter(t_temp,rad2deg(path_m(idx+1)),'*r');
end
if (POLY_TYPE == BEZIER_POLY)
    sum_t=0;
    for idx=1:n_seg
        if(mod(idx,2))
            for idy=1:(n_order+1)
                plot([sum_t sum_t+ts(idx)],rad2deg([bezier_coef_q((idx-1)*(n_order+1)+idy) bezier_coef_q((idx-1)*(n_order+1)+idy)]),'g');
            end
        else
            for idy=1:8
                plot([sum_t sum_t+ts(idx)],rad2deg([bezier_coef_q((idx-1)*(n_order+1)+idy) bezier_coef_q((idx-1)*(n_order+1)+idy)]),'c');
            end
        end
        sum_t = sum_t + ts(idx);
    end
end
grid on
xlim([0 ceil(T)]);
ylim([-200 200]);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% xy acceleration
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("xy acceleration",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
plot(tv,X_ddn, 'r-');
plot(tv,Y_ddn, 'b-');
t_temp=0;
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp,X_ddn(qdn_idx),'*r');
end
t_temp=0;
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp,Y_ddn(qdn_idx),'*b');
end
% legend('x acc','y acc');
grid on
xlim([0 ceil(T)]);
ylim([floor(min(Y_ddn))-2 ceil(max(Y_ddn))+2]);
line([0,ceil(T)],[pa_max,pa_max],'linestyle','--','color','c','LineWidth',2);
line([0,ceil(T)],[-pa_max,-pa_max],'linestyle','--','color','c','LineWidth',2);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% angular acceleration
n=n+1;
figure (n)
% subplot(2,3,iter);
hold on
pic_title = strcat("angular acceleration",datename);
title(pic_title);
tv=0:tstep:(k-1)*tstep;
plot(tv,Q_ddn, 'b-');
t_temp = 0;
scatter(t_temp,Q_ddn(1),'*r');
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    qdn_idx = ceil(t_temp/tstep);
    scatter(t_temp,Q_ddn(qdn_idx),'*r');
end
grid on
xlim([0 ceil(T)]);
ylim([floor(min(Q_ddn)) ceil(max(Q_ddn))]);
line([0,ceil(T)],[wa_max,wa_max],'linestyle','--','color','c','LineWidth',2);
line([0,ceil(T)],[-wa_max,-wa_max],'linestyle','--','color','c','LineWidth',2);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%%%%%%
%StepN: Minisnap求解器
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function poly_coef = MinimumPolySolver(waypoints, ts, n_seg, n_order,n_costorder,n_inputorder,OP_structure)
    if (nargin == 7)
        QP_inequality = OP_structure.QP_inequality;
    else
        QP_inequality = false;
    end
    
    % 起点约束
    start_cond = zeros(1,n_inputorder);
    start_cond(1) = waypoints(1);
    % 终点约束
    end_cond = zeros(1,n_inputorder);
    end_cond(1) = waypoints(end);
    %#####################################################
    % STEP 1: 计算Q矩阵
    Q = getQ(n_seg, n_order, n_costorder,ts);
    
    %#####################################################
    % STEP 2: 计算对应的约束矩阵A_beq
    [Aeq, beq] = getAbeq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond);
    
    if (nargin == 7)
        [Aieq, bieq] = getAbieq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond,OP_structure);
    end

    f = zeros(size(Q,1),1);
    options = optimoptions('quadprog','MaxIterations',6000,'Display','iter');
    % 求解多项式系数
    if (QP_inequality)
%         poly_coef = quadprog(Q,f,Aieq,bieq,Aeq, beq,[],[],[],options);
        poly_coef = quadprog(Q,f,Aieq,bieq,Aeq, beq,[],[],[]);
    else
%         poly_coef = quadprog(Q,f,[],[],Aeq, beq,[],[],[],options);
        poly_coef = quadprog(Q,f,[],[],Aeq, beq,[],[],[]);
    end
%     poly_coef = quadprog(Q,f,[],[],Aeq, beq);
end

function poly_coef = MinimumBesierSolver(waypoints, ts, n_seg, n_order,n_costorder,n_inputorder, v_max, a_max)
    % 起点约束
    start_cond = zeros(1,n_inputorder);
    start_cond(1) = waypoints(1);
    % 终点约束
    end_cond = zeros(1,n_inputorder);
    end_cond(1) = waypoints(end);
    %#####################################################
%     ts=ts./ts;
    % STEP 1: 计算Q矩阵
    Q = getMQM(n_seg, n_order, n_costorder,ts);
    Q_0 = nearestSPD(Q);
    %#####################################################
    % STEP 2: 计算对应的约束矩阵A_beq
%     [Aeq, beq] = getAMbeq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond);
    [Aeq, beq] = getBAbeq(n_seg, n_order, n_inputorder, waypoints, ts, start_cond, end_cond);
    %#####################################################
    % STEP 3: 获取飞行走廊中各个区域的范围以及相应的不等式约束
    
    % STEP 3.1: get corridor_range of x-axis or y-axis,
    % you can define corridor_range as [p1_min, p1_max;
    %                                   p2_min, p2_max;
    %                                   ...,
    %                                   pn_min, pn_max ];
%     corridor_range = [];
%     for k = 1:n_seg
%         corridor_range(k,:) = [corridor(axis,k) - corridor(axis+2,k), corridor(axis,k) + corridor(axis+2,k)];
%     end
    % STEP 3.2: 获取相应的不等式约束
%     [Aieq, bieq] = getAMbieq(n_seg, n_order, ts, v_max, a_max);
    [Aieq, bieq] = getBAbieq(n_seg, n_order, waypoints, ts, v_max, a_max);
    f = zeros(size(Q,1),1);
    for idx=1:n_seg
        f((idx-1)*(n_order+1)+1)=waypoints(idx);
        f((idx)*(n_order+1))=waypoints(idx+1);
    end
    % 求解多项式系数
%     poly_coef = quadprog(Q,f,Aieq, bieq, Aeq, beq);
%     poly_coef = quadprog(Q_0,f,Aieq, bieq, Aeq, beq);
%     poly_coef = quadprog(Q,f,[],[],Aeq, beq);
    poly_coef = quadprog(Q_0,f,[],[],Aeq, beq);
end