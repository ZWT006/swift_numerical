function GoalNodeList = insertGoalNode(Current_pose,Current_index,Parent_pose,Parent_index,Heris_cost,Path_cost,Parent_listidx,Current_state,Current_obvp,GoalNodeList)
%INSERTGOALNODE 此处显示有关此函数的摘要
%   此处显示详细说明
Whole_Cost = Heris_cost + Path_cost;
GoalNodeList.node_pose     = [GoalNodeList.node_pose;Current_pose];  % node pose list
GoalNodeList.node_index    = [GoalNodeList.node_index;Current_index]; % node index list
GoalNodeList.node_state    = [GoalNodeList.node_state;Current_state]; % node state list
GoalNodeList.node_obvp     = [GoalNodeList.node_obvp;Current_obvp]; % node state list
GoalNodeList.parent_pose      = [GoalNodeList.parent_pose;Parent_pose];   % parent pose list
GoalNodeList.parent_index     = [GoalNodeList.parent_index;Parent_index];  % parent index list
GoalNodeList.parent_listidx   = [GoalNodeList.parent_listidx;Parent_listidx];  % parent listidx list
GoalNodeList.heuris_cost      = [GoalNodeList.heuris_cost;Heris_cost];   % h(n)
GoalNodeList.path_cost        = [GoalNodeList.path_cost;Path_cost];     % g(n)
GoalNodeList.whole_cost       = [GoalNodeList.whole_cost;Whole_Cost];  %f(n)
end

