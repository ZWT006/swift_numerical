classdef OpenList
    %OPENLIST 图搜索算法的OpenList类
    %   用于保存等待搜寻的节点
    properties
        list_lenth;    % open list length
        node_pose;  % node pose list
        node_index; % node index list
        node_state; % node state [p v a j]
        node_obvp;  % obvp array
        parent_pose;   % parent pose list
        parent_index;  % parent index list
        parent_listidx; % parent node index in openlist for for quick indexing
        heuris_cost;   % h(n)
        path_cost;     % g(n)
        whole_cost;    % f(n)
        extend_flag;   % haven been poped ?
        fn_list;    % f(n)的排序表，降序排列 |f_cost|list_index|
        DEBUG;      %debug flag
    end
    
    methods
        %% Construct function
        function OP = OpenList()
            %定义节点的坐标和类型
            OP.list_lenth = 0;
        end
        
        %% if openlist all extend?
        function flag = IsAllExtend(OP)
            %判断OpenList是否全部由已经extend yes -> true; no -> false;
            [list_length,~] = size(OP.fn_list);
            if(list_length == 0)
                flag = true;
            else
                flag = false;
            end
        end
        %% if node is in openlist? 
        function list_index = isOpenList(OP,Current_index)
            %判断该点是否已经在index_list中 是 -> idx; 否 -> 0;
            list_index = 0;
            [fn_length,~] = size(OP.node_index);
            for idx = 1:fn_length
                if (OP.node_index(idx,1) == Current_index(1) && ...
                    OP.node_index(idx,2) == Current_index(2) && ...
                    OP.node_index(idx,3) == Current_index(3))
                    list_index = idx;
                    return;
                end
            end
        end
%          %% Pop Index pose
%          function Pose = PosePoP(OP,index)
%             % Pop the index pose
%             % 去除min fn node index and pop it
%            Pose = OP.node_pose(index);
%         end
        %% Pop min fn node
        function [OP,NodePop] = MinPoP(OP)
            % Pop the min fn node
            % 去除min fn node index and pop it
            min_idx = OP.fn_list(1,2);
            OP.fn_list(1,:) = [];
            NodePop.Pose  = OP.node_pose(min_idx,:);
            NodePop.Index = OP.node_index(min_idx,:);
            NodePop.State = OP.node_state(min_idx,:);
            NodePop.Path_cost = OP.path_cost(min_idx,:);
            NodePop.listidx = min_idx;
        end
        %% Insert node
        function OP=InsertNode(OP,Current_pose,Current_index,Parent_pose,Parent_index,Heris_cost,Path_cost,Parent_listidx,Current_state,Current_obvp)
            %OpenList中添加新节点
            % Step1:判断是否已经有这个节点
            list_index = isOpenList(OP,Current_index);
            Whole_Cost = Heris_cost + Path_cost;
            % 如果没有这个节点 就将该节点插入
            if(list_index == 0)
                OP.list_lenth = OP.list_lenth + 1;
                OP.node_pose     = [OP.node_pose;Current_pose];  % node pose list
                OP.node_index    = [OP.node_index;Current_index]; % node index list
                OP.node_state    = [OP.node_state;Current_state]; % node state list
                OP.node_obvp     = [OP.node_obvp;Current_obvp]; % node state list
                OP.parent_pose      = [OP.parent_pose;Parent_pose];   % parent pose list
                OP.parent_index     = [OP.parent_index;Parent_index];  % parent index list
                OP.parent_listidx   = [OP.parent_listidx;Parent_listidx];  % parent listidx list
                OP.heuris_cost      = [OP.heuris_cost;Heris_cost];   % h(n)
                OP.path_cost        = [OP.path_cost;Path_cost];     % g(n)
                OP.whole_cost       = [OP.whole_cost;Whole_Cost];  %f(n)
                [listlength_old,~] = size(OP.fn_list);
                OP = InsertFnlist(OP,Whole_Cost,OP.list_lenth);
                [listlength_now,~] = size(OP.fn_list);
                if (listlength_now == listlength_old)
                    fprintf("ERROR: OP.InsertNode GG \n");
                    save OpenList OP.fn_list OP.node_index;
                    fprintf("LOG: can not insert fn_list, the length not change \n");
                end
            else
                % 如果现有的node的fn比现在节点的大，那就把它更新了
                if (OP.whole_cost(list_index) > Whole_Cost)
                    % 删除该index的fn_list
                    [listlength_old,~] = size(OP.fn_list);

                    OP = DeletFnlist(OP,list_index);
                    
                    [listlength_now,~] = size(OP.fn_list);
                    if (listlength_now == listlength_old)
                        fprintf("ERROR: OP.InsertNode GG \n");
                        save OpenList OP.fn_list OP.node_index;
                        fprintf("LOG: can not delet fn_list data, the length not change \n");
                    end
                    % 重新添加
                    OP = InsertFnlist(OP,Whole_Cost,list_index);
                    % 更新该点的信息
                    OP.node_state(list_index,:)       = Current_state;
                    OP.node_pose(list_index,:)        = Current_pose;
                    OP.node_obvp(list_index)          = Current_obvp;   % parent pose list
                    OP.parent_pose(list_index,:)      = Parent_pose;   % parent pose list
                    OP.parent_index(list_index,:)     = Parent_index;  % parent index list
                    OP.parent_listidx(list_index,:)   = Parent_listidx;  % parent listidx list
                    OP.heuris_cost(list_index,:)      = Heris_cost;   % h(n)
                    OP.path_cost(list_index,:)        = Path_cost;     % g(n)
                    OP.whole_cost(list_index,:)       = Whole_Cost;  %f(n)
                end
            end
        end
    end
end

%% min cost list insert
function OP = InsertFnlist(OP,Whole_Cost,ListIdenx)
[list_length,~] = size(OP.fn_list);
for idx=1:list_length
    if (OP.fn_list(idx,1) >= Whole_Cost)
        % 此时fn_list中不可能有一样 listidx的点，因为已经进行过判断了
        if (OP.fn_list(idx,2) == ListIdenx)
            fprintf("ERROR: OP.InsertFnlist GG \n");
            save OpenList OP.fn_list OP.node_index;
            fprintf("LOG: fn_list have error index %d \n",ListIdenx);
        end
        % 如果Whole_Cost比较小，就插入该列表
        OP.fn_list =  [OP.fn_list(1:(idx-1),:);Whole_Cost,ListIdenx;OP.fn_list(idx:end,:)];
        return;
    end
end
% 如果Whole_Cost最大，就加在列表末尾
OP.fn_list =  [OP.fn_list;Whole_Cost,ListIdenx];
end
%% min cost list delet special one
function OP = DeletFnlist(OP,ListIdenx)
% fn_list中某个节点的fn_cost要更新，重排该列表
[list_length,~] = size(OP.fn_list);
for idx=1:list_length
    if (OP.fn_list(idx,2) == ListIdenx)
        OP.fn_list(idx,:) = []; % 清除该点赋空值，队列
        return;
    end
end
end