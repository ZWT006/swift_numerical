% n_costorder   = 4;          % 最小化的求导阶次 0=posi;1=vel;2=acc;3=jerk;4=snap;
% n_inputorder  = 4;          % 输入的阶次 可以理解为segment 之间满足等式约束的阶次

close all
clear
clc

iter = 0;
for n_costorder=2:4
    for n_inputorder = n_costorder:4
iter = iter +1;

filenema = strcat("opt_coeff_",num2str(n_costorder),"_cost_",num2str(n_inputorder),"_input");
datename = strcat(num2str(n_costorder),"-cost-",num2str(n_inputorder),"-input");
% save(filenema,'poly_coef_q', 'poly_coef_x', 'poly_coef_y');
load(filenema);

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
k = 1;
tstep = 0.01;
for i=0:n_seg-1
    %#####################################################
    % STEP 3: get the coefficients of i-th segment of both x-axis
    % and y-axis
    Pxi = poly_coef_x((n_order+1)*(i)+1:(n_order+1)*(i)+n_order+1); 
    Pyi = poly_coef_y((n_order+1)*(i)+1:(n_order+1)*(i)+n_order+1);
    Pqi = poly_coef_q((n_order+1)*(i)+1:(n_order+1)*(i)+n_order+1);
    
    for t = 0:tstep:ts(i+1)
        % position
        px=flip(Pxi);
        py=flip(Pyi);
        pq=flip(Pqi);
        X_n(k)  = polyval(px, t);
        Y_n(k)  = polyval(py, t);
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
end
k=k-1;

%%%%%%%%
%Step4: 绘制结果图例
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%<------------------------------------------------------------------->%%%
% xy map & palstance
figure (1)
subplot(2,3,iter);
hold on
pic_title = strcat("xy map & palstance",datename);
subtitle(pic_title);
AngUNIT=50;
for idx=1:k
    detl = Q_dn(idx)*AngUNIT;
    theta = Q_n(idx);
    plot([X_n(idx),X_n(idx)+cos(theta)*detl],[Y_n(idx),Y_n(idx)+sin(theta)*detl],'c-','LineWidth',1);
end
plot(X_n, Y_n , 'r-');
scatter(path(1:size(path, 1), 1), path(1:size(path, 1), 2),'ok');
grid on
[path_length,~] = size(path);
for idx=1:path_length
    plot([path(idx,1),path(idx,1)+cos(path(idx,3))*AngUNIT],[path(idx,2),path(idx,2)+sin(path(idx,3))*AngUNIT],'-','Color','b','LineWidth',1);
end
hold off
axis equal
xlim([0 1000]);
ylim([0 800]);
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% xy velocity
figure (2)
subplot(2,3,iter);
hold on
pic_title = strcat("xy velocity",datename);
subtitle(pic_title);
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
legend('x vel','y vel');
grid on
xlim([0 ceil(T)]);
ylim([-150 250]);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% angular velocity
figure(3)
subplot(2,3,iter);
hold on
pic_title = strcat("angular velocity",datename);
subtitle(pic_title);
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
ylim([-1.5 0.5]);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% angular position
figure(4)
subplot(2,3,iter);
hold on
pic_title = strcat("angular position",datename);
subtitle(pic_title);
tv=0:tstep:(k-1)*tstep;
plot(tv,rad2deg(Q_n), 'b-');
t_temp = 0;
scatter(t_temp,rad2deg(path_q(1)),'*r');
for idx=1:length(ts)
    t_temp = t_temp + ts(idx);
    scatter(t_temp,rad2deg(path_q(idx+1)),'*r');
end
grid on
xlim([0 ceil(T)]);
ylim([-60 65]);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% xy acceleration
figure (5)
subplot(2,3,iter);
hold on
pic_title = strcat("xy acceleration",datename);
subtitle(pic_title);
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
legend('x acc','y acc');
grid on
xlim([0 ceil(T)]);
ylim([-400 400]);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);

%%%<------------------------------------------------------------------->%%%
% angular acceleration
figure(6)
subplot(2,3,iter);
hold on
pic_title = strcat("angular acceleration",datename);
subtitle(pic_title);
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
ylim([-2 2]);
hold off
set(gcf,'Position', [100, 100, 1300, 700]);
    end % end n_inputorder
end % end n_costorder

%%%%%%%%%%%%%%%%%%%%%%%%%
% pic_name = strcat('TrajGener\','angular acceleration','.eps');
% saveas(gcf,pic_name,'epsc');
%保存图片
figure(1)
pic_name = strcat('TrajGener\','xy map & palstance','.eps');
saveas(gcf,pic_name,'epsc');
figure(2)
pic_name = strcat('TrajGener\','xy velocity','.eps');
saveas(gcf,pic_name,'epsc');
figure(3)
pic_name = strcat('TrajGener\','angular velocity','.eps');
saveas(gcf,pic_name,'epsc');
figure(4)
pic_name = strcat('TrajGener\','angular position','.eps');
saveas(gcf,pic_name,'epsc');
figure(5)
pic_name = strcat('TrajGener\','xy acceleration','.eps');
saveas(gcf,pic_name,'epsc');
figure(6)
pic_name = strcat('TrajGener\','angular acceleration','.eps');
saveas(gcf,pic_name,'epsc');

%%%%%%%%
%StepN: Minisnap求解器
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function poly_coef = MinimumQPSolver(waypoints, ts, n_seg, n_order,n_costorder,n_inputorder)
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
    
    f = zeros(size(Q,1),1);
    % 求解多项式系数
    poly_coef = quadprog(Q,f,[],[],Aeq, beq);
end