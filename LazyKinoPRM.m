%***************************************
%Author: Wentao Zhang
%Date: 2023-3-1
%E-mail: zwt190315@163.com
%Reference: https://www.cnblogs.com/tiandsp/p/12248779.html
%Problem: 正常的PRM要求所有节点的有向图，如果有N个节点的有向图，就得构建一个N*N的矩阵
% 来迭代进行Dijkstar算法来操作这个矩阵，耗时比较长并且随着采样点的增多，时间几何倍提升
% 还有一个随机采点的问题，造成结果的不确定性。从测试来看，貌似给固定的采样点影响也不大
% 但是在进行这种操作之前，最好是先把地图给腐蚀膨胀一下，拓展一下obstacle的外径，即使擦边也不会碰撞
% 从最终的结果来看，prm的最短路径效果是真的很不错的，是真的短。
% TODO：感觉可以提升或者可以尝试的地方 
% 1): 随机撒点和规律撒点，有啥区别吗？是不是也可以考虑随机和规律结合的那种？限制在栅格中撒点
% 2): 采样点Point类中的状态保存,修改的属性和类需要完善,目前存在问题,需要完善
% 3): 最主要的还是BVP求解轨迹添加进入搜索中
% 4): OBVP已经添加进入搜索,但是一些细节还没有处理:终止是两端状态都确定应该搞成一个BVP问题求解
% 5): 添加终点多个备选方案,优化最终结果,造成问题就是存在找不到解的情况,暂时想到的解决方法 Ⅰ)修改参数再搜索(耗时)
% Ⅱ)增加直连的标志,可以有不是那么好的解决方法
%***************************************

close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% search debug flag
SEARCH_DEBUG = false;

% 实时刷新图像
PLOT_DEBUG = true;

% 打印路径的状态信息,速度,加速度等
OBVP_DEBUG = false;

% 固定采样点,观察obvp效果
ESTEND_FLAG = false;

% 随机采样
SAMPLE_RANDOM = true;

% MAP erode
MAP_ERODE = true;
erode_range = 20; % pixel
% figure sequence
n = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('LazyPRM\');
addpath('LocalOBVP\');
addpath('map\');
addpath('TrajGener\');
addpath('PolyOpt\');
% addpath('TrajOpt\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local function difine
angle = @(x1,y1,x2,y2) atan2((y2-y1),(x2-x1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 加载地图
mapaddress = "map18.png"; % 有多个可选地图位于/map 文件夹
img = imread(mapaddress);   %空间地图
img = rgb2gray(img);   %空间地图转换为灰度图   
% img = transpose(img); 
initsdfmap = sdfMap(img);   % 导入计算 sdf
if (MAP_ERODE)
    se=strel('disk',erode_range);
    img=imerode(img, se);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 设置图层
fp = figure(n);
hold on
grid on
axis equal
axis on
set(gca,'Ydir','normal');
% title('xy position')
% axis on
grid on

if (OBVP_DEBUG)
    fv = figure(2);
    grid on
    hold on
    title('xy velocity')
    fa = figure(3);
    grid on
    hold on
    title('xy acceleration')
    fqp = figure(4);
    grid on
    hold on
    title('theta position')
    fqv = figure(5);
    grid on
    hold on
    title('theta velocity')
    fqa = figure(6);
    grid on
    hold on
    title('theta acceleration')
end

figure(fp)
% imshow(img);
initsdfmap.showSDFMap(fp);
hold on
grid on
axis equal
% axis on
set(gca,'Ydir','normal');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 设置起始点
% 行=y 列=x
[h,w,~]=size(img);
% p=ginput();                 %选取起始与结束位置
% p=[50,50;700,700];
if (mapaddress == "map4.png" || mapaddress == "map5.png" || ...
    mapaddress == "map9.png" || mapaddress == "map10.png")
    Point_i=[130,225,deg2rad(-45)];
%     Point_i=[200,280,deg2rad(90)];
%     Point_i=[540,340,deg2rad(45)];
%     Point_i=[200,450,deg2rad(90)];

%     Point_f=[935,580,deg2rad(-45)];
%     Point_f=[942,542,deg2rad(135)];
    Point_f=[130,580,deg2rad(45)];
%     Point_f=[200,200,deg2rad(-45)];
%     Point_f=[300,550,deg2rad(0)];

%     Point_i=[535,250,deg2rad(45)];
%     Point_f=[1000,100,deg2rad(-45)];
    xy_grid = 30;          %采样栅格
    xy_resolution = 20;    %随机采样的xy分辨率
    aix_resolution = 1;    %yaw角采样分辨率
    YAW_BIAS = pi;       % yaw角拓展的最大偏差
    DIR_GRID_M = 3;          % 节点拓展的分辨率
elseif (mapaddress == "map6.png")
    Point_i=[50,50,deg2rad(45)];
%     Point_f=[907,582,deg2rad(-45)];
    Point_f=[940,550,deg2rad(-45)];
elseif (mapaddress == "map1.png" || mapaddress == "map2.png" || ...
        mapaddress == "map3.png" || mapaddress == "map0.png")
    Point_i=[50,50,deg2rad(45)];
%     Point_f=[907,582,deg2rad(-45)];
    Point_f=[750,750,deg2rad(-45)];
    xy_grid = 30;          %采样栅格
    xy_resolution = 15;    %随机采样的xy分辨率
    aix_resolution = 1;    %yaw角采样分辨率
    YAW_BIAS = pi;       % yaw角拓展的最大偏差
    DIR_GRID_M = 3;          % 节点拓展的分辨率
elseif (mapaddress == "map7.png" || mapaddress == "map8.png")
    Point_i=[50,50,deg2rad(45)];
    Point_f=[2513,136,deg2rad(-45)];
    xy_grid = 40;          %采样栅格
    xy_resolution = 20;    %随机采样的xy分辨率
    aix_resolution = 1;    %yaw角采样分辨率
    YAW_BIAS = pi;       % yaw角拓展的最大偏差
    DIR_GRID_M = 3;          % 节点拓展的分辨率
elseif (mapaddress == "map11.png" || mapaddress == "map12.png")
    Point_i=[1000,1000,deg2rad(45)];
    Point_f=[300,1500,deg2rad(-45)];
    xy_grid = 30;          %采样栅格
    xy_resolution = 20;    %随机采样的xy分辨率
    aix_resolution = 1;    %yaw角采样分辨率
    YAW_BIAS = pi;       % yaw角拓展的最大偏差
    DIR_GRID_M = 3;          % 节点拓展的分辨率
elseif (mapaddress == "map13.png" || mapaddress == "map14.png" || ...
        mapaddress == "map15.png" || mapaddress == "map16.png" || ...
        mapaddress == "map17.png" || mapaddress == "map18.png" || ...
        mapaddress == "map20.png")

%     Point_i=[130,225,deg2rad(0)];
%     Point_f=[130,580,deg2rad(180)];

%     Point_i=[130,300,deg2rad(0)];
%     Point_f=[130,500,deg2rad(180)];

%     Point_i=[130,150,deg2rad(0)];
%     Point_f=[130,650,deg2rad(180)];

%     Point_i=[150,150,deg2rad(0)];
%     Point_f=[650,650,deg2rad(90)];

%     Point_i=[400,150,deg2rad(0)];
%     Point_f=[400,650,deg2rad(180)];

    Point_i=[150,650,deg2rad(0)];
    Point_f=[1050,150,deg2rad(0)];

    xy_grid = 30;          %采样栅格
    xy_resolution = 20;    %随机采样的xy分辨率
    aix_resolution = 1;    %yaw角采样分辨率
    YAW_BIAS = pi;       % yaw角拓展的最大偏差
    DIR_GRID_M = 3;          % 节点拓展的分辨率
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obvp的参数
Vel_factor = 1.4; % x,y 参考线速度
W_factor   = 2; % q 参考角速度
Racc       = 1; % cost = x + y + Racc*q (SampleOBVP计算的代价权重)
dt         = 0.01;
RATIO      = 100;   % 图片pixel 与真实距离 m 的比例 100 pixel = 1 m
% 计算node的sdf代价的阈值
dist_th    = 0.4;
heur_factor = 1.2;
% 是否使用sdf代价选择node(使node与obstacles保持一定距离)
SDF_COST   = false;
if (~SDF_COST)
heur_factor = 1; % 1.2  启发项的权重,影响搜索过程中对goal的趋近程度
end
sdf_factor  = 1.5; % 计算 trajectory cost 时 sdf 的影响程度,具体作用方式见 getSDFcost
c_angle = xy_grid/pi*0.8;   % 计算 trajectory cost 时 yaw 角的 cost 权重
GoalNum = 1;

plot(Point_i(1),Point_i(2),'o','Color','g','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
plot([Point_i(1),Point_i(1)+cos(Point_i(3))*12],[Point_i(2),Point_i(2)+sin(Point_i(3))*12],'-','Color','k','LineWidth',1);
plot(Point_f(1),Point_f(2),'o','Color','y','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','y');
plot([Point_f(1),Point_f(1)+cos(Point_f(3))*12],[Point_f(2),Point_f(2)+sin(Point_f(3))*12],'-','Color','k','LineWidth',1);
%实时刷新图像
if (PLOT_DEBUG)
    drawnow;
end
%% Step1: 空间中采样初始点
x_points = floor(w/xy_grid);
y_points = floor(h/xy_grid);
% point map的边界
MAX_MAP_Y = y_points;
MAX_MAP_X = x_points;
pose_map(y_points,x_points,aix_resolution) = Point();
pose_x=0;pose_y=0;pose_a=0;
% indexes = @(rows,cols,GIRD) (rows-1)*GIRD+cols;
% 有点特殊 但是MATLAB的坐标系也是正常的 y=行=rows; x=列=cols
%%%% sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_temp=1:y_points
    for j_temp=1:x_points
        for k_temp=1:aix_resolution
            if(SAMPLE_RANDOM)
                pose_y = (i_temp-1)*xy_grid + (rand()-0.5)*xy_resolution + xy_grid/2;
                pose_x = (j_temp-1)*xy_grid + (rand()-0.5)*xy_resolution + xy_grid/2;
                rand_pose = [pose_x,pose_y,0];
            else
                pose_y = (i_temp-1)*xy_grid + k_temp*xy_resolution;
                pose_x = (j_temp-1)*xy_grid + k_temp*xy_resolution;
                rand_pose = [pose_x,pose_y,0];
            end
            % sample yaw angle 没啥用,所以point的深度也进行xy采样
%             rand_pose(3) = rand()*2*pi;
            % 255 就是图像中的白色
            pose_map(i_temp,j_temp,k_temp)=pose_map(i_temp,j_temp,k_temp).setGrid(DIR_GRID_M,aix_resolution);
            [sdf,~,~] = initsdfmap.getDistAndGrad(rand_pose(1)/RATIO,rand_pose(2)/RATIO);
            pose_map(i_temp,j_temp,k_temp)=pose_map(i_temp,j_temp,k_temp).setSDF(sdf);
            if(img(floor(pose_y)+1,floor(pose_x)+1)==255)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 调整采样xy也随机
                % 无障碍的点为Type=M
                pose_map(i_temp,j_temp,k_temp)=pose_map(i_temp,j_temp,k_temp).setType('M');
                pose_map(i_temp,j_temp,k_temp)=pose_map(i_temp,j_temp,k_temp).setPose(rand_pose,[i_temp,j_temp,k_temp]);
            else
                % 有障碍的点为Type=W
                pose_map(i_temp,j_temp,k_temp)=pose_map(i_temp,j_temp,k_temp).setType('W');
                pose_map(i_temp,j_temp,k_temp)=pose_map(i_temp,j_temp,k_temp).setPose(rand_pose,[i_temp,j_temp,k_temp]);
            end
        end % end k_temp
    end % end j_temp
end % end i_temp
%均匀空间撒点

if(ESTEND_FLAG)
    load('pose_map.mat');
end

% 将起点和终点根据坐标放入pose_map中
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_temp = floor(Point_i(2)/xy_grid)+1;
j_temp = floor(Point_i(1)/xy_grid)+1;
start_index = [i_temp,j_temp,1];
pose_map(i_temp,j_temp,1)=pose_map(i_temp,j_temp,1).setPose(Point_i,[i_temp,j_temp,1]);
pose_map(i_temp,j_temp,1)=pose_map(i_temp,j_temp,1).setType('S');
% start point 所在索引的其它point就是无效点'W'
for k_temp = 2:aix_resolution
    pose_map(i_temp,j_temp,k_temp)=pose_map(i_temp,j_temp,k_temp).setType('W');
end
i_temp = floor(Point_f(2)/xy_grid)+1;
j_temp = floor(Point_f(1)/xy_grid)+1;
goal_index = [i_temp,j_temp,1];
pose_map(i_temp,j_temp,1)=pose_map(i_temp,j_temp,1).setPose(Point_f,[i_temp,j_temp,1]);
pose_map(i_temp,j_temp,1)=pose_map(i_temp,j_temp,1).setType('G');
% goal point 所在索引的其它point就是无效点'W'
for k_temp = 2:aix_resolution
    pose_map(i_temp,j_temp,k_temp)=pose_map(i_temp,j_temp,k_temp).setType('W');
end

%%%% plot sample points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_temp=1:y_points
    for j_temp=1: x_points
        for k_temp=1: aix_resolution 
        if (pose_map(i_temp,j_temp,k_temp).type == 'M')
            pose_x = pose_map(i_temp,j_temp,k_temp).pose(1);
            pose_y = pose_map(i_temp,j_temp,k_temp).pose(2);
            scatter(pose_x,pose_y,'b');
%         elseif(pose_map(i_temp,j_temp,k_temp).type == 'W')
%             pose_x = pose_map(i_temp,j_temp,k_temp).pose(1);
%             pose_y = pose_map(i_temp,j_temp,k_temp).pose(2);
%             scatter(pose_x,pose_y,'w');
        end
        end % end k_temp
    end % end j_temp
end % end i_temp
%实时刷新图像
if (PLOT_DEBUG)
    drawnow;
end
%% Step2: 搜索算法初始化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 初始化obvp
traj_cost=[];
obvpsg = SampleOBVP(Vel_factor/2,W_factor/2,Racc,dt);
obvpmid = SampleOBVP(Vel_factor,W_factor,Racc,dt);
% 初始化openlist和closelist
openlist  = OpenList();
closelist = CloseList();
goal_cost = DistanceCost(Point_i,Point_f) * heur_factor;
distance_cost = 0;
new_state=[0,0,0];
openlist = openlist.InsertNode(Point_i,start_index,Point_i,start_index,goal_cost,distance_cost,1,new_state,obvpsg);
closelist = closelist.InsertNode(Point_i,start_index);
NoPath=1;
final_index = -1;

% 多条达到终点的路线找最优的
find_goal = 0;
goal_num = 0;
GoalNodeList.node_pose     = [];  % node pose list
GoalNodeList.node_index    = []; % node index list
GoalNodeList.node_state    = []; % node state list
GoalNodeList.node_obvp     = []; % node state list
GoalNodeList.parent_pose      = [];   % parent pose list
GoalNodeList.parent_index     = [];  % parent index list
GoalNodeList.parent_listidx   = [];  % parent listidx list
GoalNodeList.heuris_cost      = [];   % h(n)
GoalNodeList.path_cost        = [];     % g(n)
GoalNodeList.whole_cost       = [];

% [obvpsg,xt,yt,qt,svf] = obvpsg.SolveMiniAccInputAcc(spi,new_state,spf);
% [obvpmid,xt,yt,qt,svf] = obvpmid.SolveMiniAccInputAcc(spi,new_state,spf);

%% Step3: A_star search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stop_flag = 0;
iter = 0;
MAX_SEARCH = 600;
while(true)
    iter = iter + 1;
    search_plot_flag = 1;
    if (mod(iter,100) == 0)
        fprintf("Having search %d times! \n",iter);
        %实时刷新图像
        if (PLOT_DEBUG)
            drawnow;
        end
    end
    extendflag = openlist.IsAllExtend();
    if(iter > MAX_SEARCH || stop_flag ==1 || extendflag)
        fprintf("END!!! search %d times! \n",iter);
        if (extendflag)
            fprintf("OpenList having extend \n");
        end
        break;
    end
    [openlist,NodePop] = openlist.MinPoP();
    parent_node = pose_map(NodePop.Index(1),NodePop.Index(2),NodePop.Index(3));
    closelist = closelist.InsertNode(parent_node.pose,parent_node.index);
    path_cost = NodePop.Path_cost;
    %%%%% 3.1 搜索节点
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (PLOT_DEBUG)
        drawnow;
    end
    extend_num = 0;
    abandon_num = 0;
    if (iter < DIR_GRID_M)
        DIR_GRID=iter;
    else
        DIR_GRID = DIR_GRID_M;
    end
    for dir_row = -DIR_GRID:1:DIR_GRID
        for dir_col = -DIR_GRID:1:DIR_GRID
            if(iter > MAX_SEARCH || stop_flag ==1 || extendflag)
                fprintf("END!!! search %d times! \n",iter);
                if (extendflag)
                    fprintf("OpenList having extend \n");
                end
                break;
            end
            %<----------------->%
            search_plot_flag = 1;
            idx = NodePop.Index;
            idx = idx + [dir_row,dir_col,0];
            %如果是自身索引,超出pose_map边界,就跳过,,
            if((dir_row == 0 && dir_col == 0) || idx(1) <= 0 || idx(1) > MAX_MAP_Y...
                    ||   idx(2) <= 0  || idx(2) > MAX_MAP_X)
                %自身node不会extend
                continue;
            end
            for dir_deep = 1:1:aix_resolution
                search_plot_flag = 1;
                idx(3) = dir_deep;
                new_node = pose_map(idx(1),idx(2),idx(3));
                %如果搜索到终点
                if (new_node.type == 'G')
                    find_goal = find_goal+1;
                    fprintf("find the %d goal! let us see! \n",find_goal);
                    % 终点不需要计算方向角度;终点角度固定
%                     [new_node,deflection]=new_node.getDeflection(NodePop.Pose,NodePop.Index);
                    %                     xy_cost = DistanceCost(parent_node.pose,new_node.pose);
                    spi = NodePop.Pose/RATIO;spi(3)=spi(3)*RATIO;
                    spf = new_node.pose/RATIO;     spf(3)=spf(3)*RATIO;
                    svi = NodePop.State;
                    [obvpsg,xt,yt,qt,svf] = obvpsg.SolveBVP(spi,svi,spf,[0,0,0]);
                    new_state = svf;
                    xy_cost       = obvpsg.path_length*RATIO;
                    theta_cost    = AngleCost(parent_node.pose,new_node.pose,c_angle);
                    % sdf cost ########################################
                    distance_cost = xy_cost + theta_cost;
                    if (SDF_COST)
                        distance_cost = getSDFcost(sdf_factor,dist_th,parent_node.sdf,new_node.sdf,distance_cost,xy_cost);
                    end
%                     distance_cost = xy_cost;
                    goal_cost     = DistanceCost(new_node.pose,Point_f) * heur_factor;
%                     goal_cost     = DistanceCost(new_node.pose,Point_f)*2;
                    traj_cost = [traj_cost;xy_cost,theta_cost,xy_cost+theta_cost,distance_cost,goal_cost,path_cost];
                    collision_flag = TrajectoryCheck(xt(:,1)*RATIO,yt(:,1)*RATIO,img);
                    %                     collision_flag = CollisionCheck(new_node.pose,parent_node.pose,img);
%                     plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'r--');
%                     plot(parent_node.pose(1),parent_node.pose(2),'o','Color','m','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','m');
                    if (collision_flag) % 如果没有碰撞就正常插入
%                         openlist = openlist.InsertNode(new_node.pose,idx,parent_node.pose,NodePop.Index,goal_cost,path_cost+distance_cost,NodePop.listidx,new_state,obvpsg);
                        plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'m--');
%                         if(PLOT_DEBUG && search_plot_flag)
%                             plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'c--');
%                             %plot([parent_node.pose(1),new_node.pose(1)],[parent_node.pose(2),new_node.pose(2)],'Color','c','LineStyle','--');
%                             search_plot_flag = 0;
% %                             drawnow;
%                         end
                            NoPath = 0;
                            goal_num = goal_num + 1;
                            fprintf("find the %d goal! happy! \n",goal_num);
                            GoalNodeList = insertGoalNode(new_node.pose,idx,parent_node.pose,NodePop.Index,goal_cost,path_cost+distance_cost,NodePop.listidx,new_state,obvpsg,GoalNodeList);
                            if (goal_num >= GoalNum)
                                stop_flag = 1;
                            end
%                         NoPath = 0;
%                         stop_flag = 1;
                        break;
                    else
                        continue;
                    end
                end % end if node.type=G
                %如果该索引是无效节点,或者已经在CloseList,就跳过该循环
                if(new_node.type == 'W' || closelist.isCloseList(idx))
                    abandon_num = abandon_num + 1;
                    continue;
                end
                %%%#####################################################%%%
                %如果该索引的yaw角超过最大偏差也跳过该循环
                % 没有搜索到终点就正常进行extend 计算
                % TODO list:使用pose_map信息进行externd加速
                % 这部分计算最为耗时，BVP or OBVP or其它kinodynamic计算都在
                % 这里，还有collision check问题，利用pose_map中的point已经
                % 存储的信息，之前计算已经externd计算过的内容直接使用，来加速
                % extend的时间
                % distance_cost = DistanceCost(parent_node.pose,new_node.pose);
                if (parent_node.type == 'S')
                    [new_node,deflection]=new_node.getDeflection(NodePop.Pose,NodePop.Index);
                    spi = NodePop.Pose/RATIO;   spi(3)=spi(3)*RATIO;
                    spf = new_node.pose/RATIO;  spf(3)=spf(3)*RATIO;
                    svi = NodePop.State;
                    [obvpsg,xt,yt,qt,svf] = obvpsg.SolveMiniAccInputAcc(spi,svi,spf);
                    new_state = svf;
                    theta_cost    = AngleCost(parent_node.pose,new_node.pose,c_angle);
                    xy_cost       = obvpsg.path_length*RATIO;
                    %<----------------->%
                    % sdf cost ########################################
                    distance_cost = xy_cost + theta_cost;
                    if (SDF_COST)
                        distance_cost = getSDFcost(sdf_factor,dist_th,parent_node.sdf,new_node.sdf,distance_cost,xy_cost);
                    end
%                     distance_cost = xy_cost;
                    goal_cost     = DistanceCost(new_node.pose,Point_f) * heur_factor;
                    traj_cost = [traj_cost;xy_cost,theta_cost,xy_cost+theta_cost,distance_cost,goal_cost,path_cost];
%                     goal_cost     = DistanceCost(new_node.pose,Point_f)*2;
                    collision_flag = TrajectoryCheck(xt(:,1)*RATIO,yt(:,1)*RATIO,img);
                    if (collision_flag) % 如果没有碰撞就正常插入
                        openlist = openlist.InsertNode(new_node.pose,idx,parent_node.pose,NodePop.Index,goal_cost,path_cost+distance_cost,NodePop.listidx,new_state,obvpsg);
                        plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'c--');
%                         if(PLOT_DEBUG && search_plot_flag)
%                             plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'c--');
%                             %                         plot([parent_node.pose(1),new_node.pose(1)],[parent_node.pose(2),new_node.pose(2)],'Color','c','LineStyle','--');
%                             search_plot_flag = 0;
% %                             drawnow;
%                         end
                    end
                    if (~collision_flag)
                        %                     plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'r-');
                    end
                else
                    [new_node,deflection]=new_node.getDeflection(NodePop.Pose,NodePop.Index);
                    spi = NodePop.Pose/RATIO;   spi(3)=spi(3)*RATIO;
                    spf = new_node.pose/RATIO;  spf(3)=spf(3)*RATIO;
                    svi = NodePop.State;
                    [obvpmid,xt,yt,qt,svf] = obvpmid.SolveMiniAccInputAcc(spi,svi,spf);
                    new_state = svf;
                    theta_cost    = AngleCost(parent_node.pose,new_node.pose,c_angle);
                    xy_cost       = obvpmid.path_length*RATIO;
                    % sdf cost ########################################
                    distance_cost = xy_cost + theta_cost;
%                     distance_cost = xy_cost;
                    if (SDF_COST)
                        distance_cost = getSDFcost(sdf_factor,dist_th,parent_node.sdf,new_node.sdf,distance_cost,xy_cost);
                    end
                    goal_cost     = DistanceCost(new_node.pose,Point_f) * heur_factor;
                    traj_cost = [traj_cost;xy_cost,theta_cost,xy_cost+theta_cost,distance_cost,goal_cost,path_cost];
%                     goal_cost     = DistanceCost(new_node.pose,Point_f)*2;
                    collision_flag = TrajectoryCheck(xt(:,1)*RATIO,yt(:,1)*RATIO,img);
                    if (collision_flag) % 如果没有碰撞就正常插入
                        openlist = openlist.InsertNode(new_node.pose,idx,parent_node.pose,NodePop.Index,goal_cost,path_cost+distance_cost,NodePop.listidx,new_state,obvpmid);
                        plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'c--','LineWidth',0.2);
%                         if(PLOT_DEBUG && search_plot_flag)
%                             plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'c--');
%                             %                         plot([parent_node.pose(1),new_node.pose(1)],[parent_node.pose(2),new_node.pose(2)],'Color','c','LineStyle','--');
%                             search_plot_flag = 0;
% %                             drawnow;
%                         end
                    end
                    if (~collision_flag)
                        %                     plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'r-');
                    end
                end
                %                 parent_node = parent_node.addAdjacent(idx,distance_cost,new_node.pose,collision_flag);
                %                 new_node    = new_node.addAdjacent(NodePop.Index,distance_cost,NodePop.Pose,collision_flag);
                %把修改后的node再放入pose_map中
                if (SEARCH_DEBUG)
                    extend_num = extend_num + 1;
                end
                pose_map(idx(1),idx(2),idx(3)) = new_node;
            end % end theta extend
        end % end y direction extend
    end % end z direction extend
    % 将拓展后的parent_node再次放入pose_map
    pose_map(NodePop.Index(1),NodePop.Index(2),NodePop.Index(3)) = parent_node;
    if (SEARCH_DEBUG)
        fprintf("Search extend_num = %d; Abandon_num = %d! \n",extend_num,abandon_num);
        [fn_list_length,~] = size(openlist.fn_list);
        fprintf("openlist.fn_list length %d \n",fn_list_length);
    end

end % end A_star search loop

fprintf("Having search %d times! \n",iter);
%% Step4: Get Path 
% 将角度变化处理得更加光滑
path = [];
path_index = [];
path_sdf = [];
if NoPath
    return
else
    n = 1;

%     path_index(n,1) = openlist.list_lenth;
%     final_index = openlist.parent_listidx(openlist.list_lenth);
    
    [~,mingoalindex] = min(GoalNodeList.whole_cost);
    path_index(n,1) = mingoalindex;
    final_index = GoalNodeList.parent_listidx(mingoalindex);

    path(n,:) = Point_f;
    [temp_sdf,~,~] = initsdfmap.getDistAndGrad(Point_f(1),Point_f(2));
    path_sdf(n,:) = [temp_sdf,0,0,0];
    while(final_index > 1)
        n = n +1;
        pose = openlist.node_pose(final_index,:);
        pose_index = openlist.node_index(final_index,:);
        path_index(n) = final_index;
        final_index = openlist.parent_listidx(final_index);
        %             fprintf('the %d final_index is : %d \n',n,final_index);
        path(n,:) =pose;
        [temp_sdf,~,~] = initsdfmap.getDistAndGrad(Point_f(1),Point_f(2));
        path_sdf(n,:) = [temp_sdf,pose_index];
    end
end

path(n+1,:) = Point_i;
path_index(n+1) = 1;

[path_length,~] = size(path);
% for idx=1:path_length-1
for idx=2:path_length-1
    plot([path(idx,1),path(idx+1,1)],[path(idx,2),path(idx+1,2)],'--','Color','g','LineWidth',1);
end
% detl = 30;
% % for idx=1:path_length
% for idx=2:path_length
%     plot([path(idx,1),path(idx,1)+cos(path(idx,3))*detl],[path(idx,2),path(idx,2)+sin(path(idx,3))*detl],'-','Color','g','LineWidth',2);
% end

%% Step5: Show the path and direction

[path_length,~]=size(path);

path=flip(path);

% angleset=zeros(path_length,1);
% 
% for idx = 2:path_length
%     angleset(idx)=angle(path(idx-1,1),path(idx-1,2),path(idx,1),path(idx,2));
% end
%##########################################################################
%%%%% 偏转角平滑处理
% for idx = 1:path_length-2
%     path(idx+1,3)=(angle(path(idx,1),path(idx,2),path(idx+1,1),path(idx+1,2))+angle(path(idx+1,1),path(idx+1,2),path(idx+2,1),path(idx+2,2)))/2;
% end

% path_m=zeros(path_length,1);
% for idx = 1:path_length-2
%     temp_angle=AngleDelta(path(idx+1,3),path(idx+2,3))/2+path(idx+1,3);
%     if (abs(temp_angle)>pi)
%         if (temp_angle < -pi)
%             temp_angle = 2*pi+temp_angle;
%         else
%             temp_angle = -2*pi+temp_angle;
%         end
%     end
%     path_m(idx+1)=temp_angle;
% end
% 计算角度差值
angle_diff = AngleDelta(path(2:path_length-1, 3), path(3:path_length, 3)) / 2 + path(2:path_length-1, 3);
angle_diff_raw = AngleDelta(path(2:path_length-1, 3), path(3:path_length, 3)) / 2 + path(2:path_length-1, 3);
angle_delta = AngleDelta(path(2:path_length-1, 3), path(3:path_length, 3));
% 处理角度超过 [-pi, pi] 范围的情况
for i=1:length(angle_diff)
    if (abs(angle_diff(i)) > pi)
        if (angle_diff(i) > 0) % -0 ~ -180
            angle_diff(i) = angle_diff(i) - 2*pi;
        else % +0 ~ +180
            angle_diff(i) = angle_diff(i) + 2*pi;
        end
    end
end
% 更新 path_m 数组
path_m = [0; angle_diff; 0];

path_m(1)=path(1,3);
path_m(end)=path(end,3);
path(:,3)=path_m;
[path_length,~] = size(path);
detl = 30;
% 计算坐标
end_x = path(:, 1) + cos(path(:, 3)) * detl;
end_y = path(:, 2) + sin(path(:, 3)) * detl;
% 绘制路径线段
plot([path(:, 1) end_x]', [path(:, 2) end_y]', '-', 'Color', 'k', 'LineWidth', 2);
% for idx=1:path_length
%     plot([path(idx,1),path(idx,1)+cos(path(idx,3))*detl],[path(idx,2),path(idx,2)+sin(path(idx,3))*detl],'-','Color','k','LineWidth',2);
% end

% [xt,yt,qt] = obvpsg.st();
% figure(fp);
% plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'b-');
% 
% for idx=1:path_length
idx=1;
path_obvp = [];
PRMX_n = [];
PRMY_n = [];
PRMQ_n = [];
PRMX_dn = [];
PRMY_dn = [];
PRMQ_dn = [];
PRMX_ddn = [];
PRMY_ddn = [];
PRMQ_ddn = [];
obvptemp = GoalNodeList.node_obvp(path_index(idx));
    [xt,yt,qt] = obvptemp.st();
    figure(fp);
    plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'b-','LineWidth',2);
    PRMX_n = [xt(:,1);PRMX_n];
    PRMY_n = [yt(:,1);PRMY_n];
    PRMQ_n = [qt(:,1);PRMQ_n];
    PRMX_dn = [xt(:,2);PRMX_dn];
    PRMY_dn = [yt(:,2);PRMY_dn];
    PRMQ_dn = [qt(:,2);PRMQ_dn];
    PRMX_ddn = [xt(:,3);PRMX_ddn];
    PRMY_ddn = [yt(:,3);PRMY_ddn];
    PRMQ_ddn = [qt(:,3);PRMQ_ddn];
    path_obvp = [path_obvp;obvptemp];
for idx=2:path_length-1
    obvptemp = openlist.node_obvp(path_index(idx));
    [xt,yt,qt] = obvptemp.st();
    PRMX_n = [xt(:,1);PRMX_n];
    PRMY_n = [yt(:,1);PRMY_n];
    PRMQ_n = [qt(:,1);PRMQ_n];
    PRMX_dn = [xt(:,2);PRMX_dn];
    PRMY_dn = [yt(:,2);PRMY_dn];
    PRMQ_dn = [qt(:,2);PRMQ_dn];
    PRMX_ddn = [xt(:,3);PRMX_ddn];
    PRMY_ddn = [yt(:,3);PRMY_ddn];
    PRMQ_ddn = [qt(:,3);PRMQ_ddn];
    figure(fp);
    plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'b-','LineWidth',2);
    path_obvp = [path_obvp;obvptemp];
end
axis on;
grid on;
set(gca,'GridColor',[0.5 0.5 0.5],'GridLineStyle','--','GridAlpha',0.2,'XMinorGrid','off','YMinorGrid','off','XTick',[0:100:800],'YTick',[0:100:800]);
% save('pathnode.mat','path','path_obvp')
save('pathnode.mat','path','path_obvp','PRMX_n','PRMY_n','PRMQ_n','PRMX_dn',...
    'PRMY_dn','PRMQ_dn','PRMX_ddn','PRMY_ddn','PRMQ_ddn')
%% Step6: show the trajectory
[path_length,~]=size(path);

% path=flip(path);

statex=zeros(path_length,4);
statey=zeros(path_length,4);
stateq=zeros(path_length,4);
Tarray=zeros(path_length,1);

statev = zeros(path_length,3);
% 初始速度都是零
statev(1,:) = [0,0,0];


% angle = @(x1,y1,x2,y2) atan2((y2-y1),(x2-x1));
% 
% angleset=zeros(path_length,1);
% 
% for idx = 2:path_length
% %     path(idx+1,3)=(angle(path(idx,1),path(idx,2),path(idx+1,1),path(idx+1,2))+angle(path(idx+1,1),path(idx+1,2),path(idx+2,1),path(idx+2,2)))/2;
%     angleset(idx)=angle(path(idx-1,1),path(idx-1,2),path(idx,1),path(idx,2));
% end
% for idx = 1:path_length-2
%     path(idx+1,3)=(angle(path(idx,1),path(idx,2),path(idx+1,1),path(idx+1,2))+angle(path(idx+1,1),path(idx+1,2),path(idx+2,1),path(idx+2,2)))/2;
% end
% 
% [path_length,~] = size(path);
% for idx=1:path_length-1
%     plot([path(idx,1),path(idx+1,1)],[path(idx,2),path(idx+1,2)],'c--');%,'LineWidth',2
% end
% detl = 20;
% for idx=1:path_length
%     plot([path(idx,1),path(idx,1)+cos(path(idx,3))*detl],[path(idx,2),path(idx,2)+sin(path(idx,3))*detl],'-','Color','b','LineWidth',1);
% %     plot([path(idx,1),path(idx,1)+cos(angleset(idx))*detl],[path(idx,2),path(idx,2)+sin(angleset(idx))*detl],'-','Color','b','LineWidth',1);
% end

%%%%%%%%第一次启动速度不够,所以时间加倍
% obvp=SampleOBVP(Vel_factor/2,W_factor/2,Racc,dt);
% idx=1;
% spi = path(idx,:)/RATIO;
% spf = path(idx+1,:)/RATIO;
% spi(3)=spi(3)*RATIO;
% spf(3)=spf(3)*RATIO;
% 
% [obvp,xt,yt,qt,svf] = obvp.SolveMiniAccInputAcc(spi,statev(idx,:),spf);
% statev(idx+1,:) = svf;
% Tarray(idx) = obvp.refer_T;
% figure(fp);
% plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'--','Color','r','LineWidth',2);
% 
% if (OBVP_DEBUG)
%     figure(fv);
%     tf = sum(Tarray(1:idx));
%     ti = tf-Tarray(idx);
%     t = ti:dt:tf+dt;
%     plot(t,xt(:,2),'-','Color','b');
%     plot(t,yt(:,2),'-','Color','r');
%     figure(fa);
%     plot(t,xt(:,3),'-','Color','b');
%     plot(t,yt(:,3),'-','Color','r');
%     figure(fqp);
%     plot(t,rad2deg(qt(:,1)),'-','Color','r');
%     figure(fqv)
%     plot(t,qt(:,2),'-','Color','b');
%     figure(fqa);
%     plot(t,qt(:,3),'-','Color','b');
% end
% 
% obvp=SampleOBVP(Vel_factor,W_factor,Racc,dt);
% for idx=2:path_length-2
%     spi = path(idx,:)/RATIO;
%     spf = path(idx+1,:)/RATIO;
%     spi(3)=spi(3)*RATIO;
%     spf(3)=spf(3)*RATIO;
%     [obvp,xt,yt,qt,svf] = obvp.SolveMiniAccInputAcc(spi,statev(idx,:),spf);
%     statev(idx+1,:) = svf;
%     Tarray(idx) = obvp.refer_T;
%     figure(fp);
%     plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'--','Color','r','LineWidth',2);
%     if (OBVP_DEBUG)
%         figure(fv);
%         tf = sum(Tarray(1:idx));
%         ti = tf-Tarray(idx);
%         t = ti:dt:tf+dt;
%         plot(t,xt(:,2),'-','Color','b');
%         plot(t,yt(:,2),'-','Color','r');
%         figure(fa);
%         plot(t,xt(:,3),'-','Color','b');
%         plot(t,yt(:,3),'-','Color','r');
%         figure(fqp);
%         plot(t,rad2deg(qt(:,1)),'-','Color','r');
%         figure(fqv)
%         plot(t,qt(:,2),'-','Color','b');
%         figure(fqa);
%         plot(t,qt(:,3),'-','Color','b');
%     end
% end
% 
% idx = path_length-1;
% obvp=SampleOBVP(Vel_factor/2,W_factor/2,Racc,dt);
% spi = path(idx,:)/RATIO;
% spf = path(idx+1,:)/RATIO;
% spi(3)=spi(3)*RATIO;
% spf(3)=spf(3)*RATIO;
% [obvp,xt,yt,qt,svf] = obvp.SolveBVP(spi,statev(idx,:),spf,[0,0,0]);
% statev(idx+1,:) = svf;
% Tarray(idx) = obvp.refer_T;
% figure(fp);
% plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'--','Color','r','LineWidth',2);
% if (OBVP_DEBUG)
%     figure(fv);
%     tf = sum(Tarray(1:idx));
%     ti = tf-Tarray(idx);
%     t = ti:dt:tf+dt;
%     plot(t,xt(:,2),'-','Color','b');
%     plot(t,yt(:,2),'-','Color','r');
%     figure(fa);
%     plot(t,xt(:,3),'-','Color','b');
%     plot(t,yt(:,3),'-','Color','r');
%     figure(fqp);
%     plot(t,rad2deg(qt(:,1)),'-','Color','r');
%     figure(fqv)
%     plot(t,qt(:,2),'-','Color','b');
%     figure(fqa);
%     plot(t,qt(:,3),'-','Color','b');
% end
% 
% hold off
% fprintf("OVER !!! \n");

if(~ESTEND_FLAG)
    save('pose_map.mat','pose_map');
    saveas(gcf,'lazykinoPRM.fig');
end

% filename = "E:\datas\Swift\Optimization\TRAJ_COST.csv";
% % traj_cost = [traj_cost;xy_cost,theta_cost,xy_cost+theta_cost,distance_cost,goal_cost,path_cost];
% TRAJ_DATA = traj_cost;
% TRAJ_DATA = round(TRAJ_DATA,6);
% writematrix(TRAJ_DATA,filename);


function sdfcost = getSDFcost(sdf_factor,dist_th,parsdf,newsdf,distance_cost,xy_cost)
parsdf = dist_th - parsdf;
if (parsdf < 0)
    parsdf = 0;
end
newsdf = dist_th - newsdf;
if (newsdf < 0)
    newsdf = 0;
end
dist_factor = (newsdf + parsdf)/(2*dist_th);
sdfcost = dist_factor * distance_cost * sdf_factor + distance_cost;
end