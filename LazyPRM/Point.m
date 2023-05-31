%***************************************
%@Author: Wentao Zhang
%@Date: 2023-3-7
%@E-mail: zwt190315@163.com
%@Reference: ##########
%@Problems: 
%@Description:
%@TODO：##########
%***************************************
classdef Point
    properties
        type;       % point type
        % Start -> S;Goal -> G;Midpoint -> M;Invalid -> W;Extend -> E;
        pose;       % pose of point  [x,y,θ]
        index;      % pose_map index [rows,cols,deepth]
        extandmap;  % point extand map 2->obsticle;1->free;0->not judge
        mapgrid;    % extandmap size = mapgrid
        deepth;     % sample points deepth
        adjacent_points; % adjacent points whole 
        adjacent_indexs; % adjacent point index for quick indexing
        adjacent_poses;  % adjecent point pose
        adjacent_costs; % adjacent point path costs
        deflection_angle; % deflection angle between patent node and child node
        min_mapidx; % min cost map index
        min_cost;   % min cost
        judeg_length; % the adjacent point having beening extend
        sdf;        % Euclidean signed distance fields Information
    end
    
    methods
        %% Construct function
%         function P = Point(Pose, Type)
        function P = Point()
            %定义节点的坐标和类型
            P.type = 'W';
            P.pose = [0,0,0];
            P.min_mapidx = [0,0,0];
            %初始最小cost=Naf
            P.min_cost = 99999;
            %已拓展的邻近点为0
            P.judeg_length = 0;
        end
        %% Init Point
        function P = InitPoint(P,Adjacent_index)
            % 初始化部分节点的连接关系
            extmap_idx = getidx(Adjacent_index); %索引偏移
            P.extandmap(extmap_idx) = 0;
            %已经判断的邻近点数重置
            P.judeg_length = P.judeg_length-1;
            P.type = 'M';
        end
        %% set P grid
        function P = setGrid(P,Mapgrid,Deepth)
            % 原本是初始化的事儿,但是不会用MATLAB初始化类数组
            P.mapgrid = Mapgrid;
            P.deepth  = Deepth;
            P.adjacent_points = Mapgrid^2*Deepth;
        end
        %% set P sdf
        function P = setSDF(P,sdf)
            % 设定节点的sdf
            P.sdf = sdf;
        end
        %% set Pose function
        function P = setPose(P,Pose,Index) 
            %设定节点的Pose
            P.pose = Pose;
            P.index = Index;
        end
        %% set Type function
        function P = setType(P,Type) 
            %设定节点的Type
            P.type = Type;
        end
        %% add adjacent node function
        function P = addAdjacent(P,Adjacent_index,Adjacent_cost,Adjacent_pose,obsticle_flag) 
            %添加节点的临近节点
            extmap_idx =getidx(Adjacent_index); %九宫格索引偏移
            %<--------------------------DEBUG---------------------------->%
%           fprintf("extmap_idx= %d , %d , %d \n",extmap_idx(1),extmap_idx(2),extmap_idx(3))
            
            if (obsticle_flag)
                P.extandmap(extmap_idx) = 2;
            else
                P.extandmap(extmap_idx) = 1;
                P.adjacent_costs(extmap_idx)  = Adjacent_cost;
                P.adjacent_indexs(extmap_idx) = Adjacent_index;
                P.adjacent_poses(extmap_idx) = Adjacent_pose;
                % 判断最小cost 无碰撞的点添加才有意义
                if(Adjacent_cost < P.min_cost)
                    P.min_mapidx = extmap_idx;
                end
            end
            
            %记录已经判断的点数
            P.judeg_length = P.judeg_length + 1;
        end
        %% get deflection angle
        function [P,angle] = getDeflection(P,Pose,Index)
            % return deflection angle from Pose to obj.pose
            angle = atan2((P.pose(2)-Pose(2)),(P.pose(1)-Pose(1)));
%             extmap_idx =getidx(P,Index); %九宫格索引偏移
%             P.deflection_angle(extmap_idx)=angle;
            P.pose(3)=angle;
        end
    end
end

function idx = getidx(P,Adjacent_index)
    extmap_idx = Adjacent_index-P.index + [P.mapgrid,P.mapgrid,0]; %九宫格索引偏移
    idx=extmap_idx(1)*P.mapgrid+extmap_idx(2);
    idx = idx * P.deepth + Adjacent_index(3);
end