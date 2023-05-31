function feasible = CollisionCheck(start,target,map)
%COLLISIONCHECK 检查两点之间是否有障碍 no -> true; yes -> false
%输入为两节点坐标和地图信息，若两节点直线连接不会经过障碍物，则返回true, 碰到障碍物则为返回false
feasible=true;
startPose = [start(1),start(2)];
goalPose  = [target(1),target(2)];
%两个节点的方向
dir=atan2(goalPose(1)-startPose(1),goalPose(2)-startPose(2));
%就是一种遍历比较，以一个较小的增量从0扩展到r长度，检查这条直线上离散点是否是障碍物
for r=0:0.5:sqrt(sum((startPose-goalPose).^2))
    posCheck = startPose + r.*[sin(dir) cos(dir)];
    %ceil 向上取整; floor向下取整
    if ~(feasiblePoint(ceil(posCheck),map) && feasiblePoint(floor(posCheck),map) && ...
            feasiblePoint([ceil(posCheck(1)) floor(posCheck(2))],map) && feasiblePoint([floor(posCheck(1)) ceil(posCheck(2))],map))
        feasible=false;break;
    end
end
%比较末尾的点
if ~feasiblePoint([floor(goalPose(1)),ceil(goalPose(2))],map), feasible=false; end

    function feasible=feasiblePoint(point,map)
        feasible=true;
        if ~(point(1)>=1 && point(1)<=size(map,2) && ... % x in map
                point(2)>=1 && point(2)<=size(map,1) && ... % y in map
                map(point(2),point(1))==255) % x,y is Free(是空白点)
            feasible=false;
        end
    end
end
