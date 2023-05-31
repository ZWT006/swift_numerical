function delta = AngleDelta(start,goal)
%ANGLEDELTA 计算-180~180的角度差 start + delta = goal;
%   计算角度差,考虑方向
goal = mod(goal,2*pi);
if (goal < 0)
    goal = goal + 2*pi;
end % 按周期转化为0~2pi
start = mod(start,2*pi);
if (start < 0)
    start = start + 2*pi;
end % % 按周期转化为0~2pi
delta = goal - start;
% 有方向的转换
if (abs(delta) > pi)
    if (delta > 0) % -0 ~ -180
        delta = delta - 2*pi;
    else % +0 ~ +180
        delta = delta + 2*pi;
    end
end
end

