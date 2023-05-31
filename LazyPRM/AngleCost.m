function cost = AngleCost(start, goal,c_angle)
    %计算两个node直接的转角cost 角度使用一个c_angle权重
    % xy_grid = 40;YAW_BIAS = pi/2;
%     c_angle = 30;
    angle = AngleDelta(start(3),goal(3));
    cost = c_angle*angle^2;
end

