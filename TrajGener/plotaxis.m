% 定义坐标范围
xrange = [-5, 5];
yrange = [-5, 5];
zrange = [-5, 5];

% 绘制 x 轴
plot3([xrange(1), xrange(2)], [0, 0], [0, 0], 'r', 'LineWidth', 1.5);
hold on;

% 绘制 y 轴
plot3([0, 0], [yrange(1), yrange(2)], [0, 0], 'g', 'LineWidth', 1.5);

% 绘制 z 轴
plot3([0, 0], [0, 0], [zrange(1), zrange(2)], 'b', 'LineWidth', 1.5);

% 设置坐标范围和比例尺
xlim(xrange);
ylim(yrange);
zlim(zrange);
daspect([1, 1, 1]); % 等比例缩放

% 添加坐标轴标签
xlabel('x');
ylabel('y');
zlabel('z');
grid on
% 添加图例
legend('x-axis', 'y-axis', 'z-axis');
