function feasible = TrajectoryCheck(xtraj,ytraj,map)
%COLLISIONCHECK 检查轨迹线上是否有障碍 no -> true; yes -> false
%输入为xy轨迹和地图信息，若轨迹线不会经过障碍物，则返回true, 碰到障碍物则为返回false
feasible=true;
traj_length = length(xtraj);
% 就是一种遍历比较，判断离散的轨迹点是否经过障碍
% 这个比较方式跟地图和轨迹的分辨率有关,暂时就这样直接比较离散的空间点;
% TODO: 完善轨迹点的collisioncheck
for idx=1:traj_length
    posCheck = [xtraj(idx),ytraj(idx)];
    %ceil 向上取整; floor向下取整
    if ~(feasiblePoint(ceil(posCheck),map) && feasiblePoint(floor(posCheck),map) && ...
            feasiblePoint([ceil(posCheck(1)) floor(posCheck(2))],map) && feasiblePoint([floor(posCheck(1)) ceil(posCheck(2))],map))
        feasible=false;break;
    end
end


    function feasible=feasiblePoint(point,map)
        feasible=true;
        if ~(point(1)>=1 && point(1)<=size(map,2) && ... % x in map
                point(2)>=1 && point(2)<=size(map,1) && ... % y in map
                map(point(2),point(1))==255) % x,y is Free(是空白点)
            feasible=false;
        end
    end
end
