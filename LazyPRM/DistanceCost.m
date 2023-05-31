function cost = DistanceCost(start, goal)
    %计算两点之间的二范数 角度使用一个c_angle权重
    % xy_grid = 40;YAW_BIAS = pi/2;
    c_angle = 20;
    angle = abs(start(3)-goal(3));
    % 角度查的一个转化问题
    if (angle > pi)
        angle = start(3) + 2*pi - goal(3);
    end
    cost = sqrt((start(1) - goal(1))^2 + (start(2) - goal(2))^2 + c_angle*angle^2);
end

