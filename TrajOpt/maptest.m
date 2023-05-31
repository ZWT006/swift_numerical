% 读取地图
close all
clear
map = imread('F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\map\map5.png');
cmap = colormap;
cmap(1,:)=0;
colormap(cmap);

% imshow(map);
% 将地图转换为二值图像
bw_map = im2gray(map);
% bw_map = im2bw(map);
bw_map = imbinarize(bw_map);
bw_map = ~bw_map;
% 计算距离变换
dist_transform = bwdist(bw_map);

% 将距离变换值符号化
% sdf_map = -1 * (bw_map - 0.5) .* dist_transform + (bw_map - 0.5) .* dist_transform;

% 显示SDF距离场地图
% imshow(dist_transform, [0 100]);

% 保存SDF距离场地图
% imwrite(sdf_map, 'sdf_map.png');


D1 = bwdist(bw_map,'euclidean');
D2 = bwdist(bw_map,'cityblock');
D3 = bwdist(bw_map,'chessboard');
D4 = bwdist(bw_map,'quasi-euclidean');

RGB1 = repmat(rescale(D1), [1 1 3]);
RGB2 = repmat(rescale(D2), [1 1 3]);
RGB3 = repmat(rescale(D3), [1 1 3]);
RGB4 = repmat(rescale(D4), [1 1 3]);

% figure
% subplot(2,2,1), imshow(RGB1), title('Euclidean')
% hold on, imcontour(D1)
% subplot(2,2,2), imshow(RGB2), title('Cityblock')
% hold on, imcontour(D2)
% subplot(2,2,3), imshow(RGB3), title('Chessboard')
% hold on, imcontour(D3)
% subplot(2,2,4), imshow(RGB4), title('Quasi-Euclidean')
% hold on, imcontour(D4)

% figure
% 将灰度图归一化到[0,1]区间
I_norm = mat2gray(dist_transform);

% 生成黄绿色渐进色彩映射表
% cmap = parula(8);
colormap(cool);
cmap = colormap;
% colormap(cmap);
% 将灰度图映射到黄绿色渐进色彩映射表上，得到对应的RGB颜色值
% 将颜色梯度翻转
% rgb = interp1(flip(linspace(0, 1, 256)), cmap, I_norm);
rgb = interp1(linspace(0, 1, 256), cmap, I_norm);
% imshow(rgb);
% imshow(map);
% 将RGB颜色值与灰度图叠加，得到最终的彩色图
set(gca,'Ydir','normal');
imshow(rgb.*repmat(~bw_map,[1 1 3]))
title('Euclidean')
hold on
% imcontour(D1,'c-')
imcontour(D1)


bw_map = im2gray(map);
% bw_map = im2bw(map);
bw_map = imbinarize(bw_map);
bw_map = ~bw_map;
% 计算距离变换
dist_transform = bwdist(bw_map);
D1 = bwdist(bw_map,'euclidean');
% 将灰度图归一化到[0,1]区间
I_norm = mat2gray(dist_transform);
% 生成黄绿色渐进色彩映射表
% cmap = parula(8);
colormap(cool);
cmap = colormap;
% 将灰度图映射到黄绿色渐进色彩映射表上，得到对应的RGB颜色值
% 将颜色梯度翻转
rgb = interp1(linspace(0, 1, 256), cmap, I_norm);
% 将RGB颜色值与灰度图叠加，得到最终的彩色图
set(gca,'Ydir','normal');
imshow(rgb.*repmat(~bw_map,[1 1 3]))
title('Euclidean')
hold on
% imcontour(D1,'c-')
imcontour(D1)