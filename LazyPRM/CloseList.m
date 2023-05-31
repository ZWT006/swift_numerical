classdef CloseList
    %CLOSELIST 图搜索算法的CloseList
    %   用于保存已经搜索过的节点
    
    properties
        list_length;    % close list length
        node_pose;      % node pose list
        node_index;     % node index list
    end
    
    methods
        %% 构造函数
        function obj = CloseList()
            %CLOSELIST 构造此类的实例
            %   此处显示详细说明
            obj.list_length = 0;
        end
        %% 插入新节点
        function obj = InsertNode(obj,Pose,Index)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            obj.list_length = obj.list_length + 1;
            obj.node_pose = [obj.node_pose;Pose];
            obj.node_index = [obj.node_index;Index];
        end
        %% 插入新节点
        function flag = isCloseList(obj,Index)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            flag = false;
            for idx = 1:obj.list_length
               if (Index(1)==obj.node_index(idx,1) && ...
                   Index(2)==obj.node_index(idx,2) && ...
                   Index(3)==obj.node_index(idx,3) )
                   flag = true;
                   return;
               end
            end
        end
    end
end

