classdef sdfMap
    %SDFMAP sdf 地图信息相关功能函数
    %   此处显示详细说明
    
    properties
        map;            % 原始地图
        bw_map;         % 二值图
        D1;             % 等距线
        rgb;            % RGB图像
        I_norm;         % 归一化 I_norm
        dist_transform; % edf 地图
        resolution;     % 分辨率是
        RATION;         % 分辨率倒数
        d_th; % 距离代价的阈值
        cols;
        rows;
    end
    
    methods
        function obj = sdfMap(map)
            %SDFMAP 
            %   矩阵初始化
            obj.map = map;
            bw_map = im2gray(map);
            bw_map = imbinarize(bw_map);
            bw_map = ~bw_map;
            dist_transform = bwdist(bw_map);
            D1 = bwdist(bw_map,'euclidean');
            I_norm = mat2gray(dist_transform);
            cmap = bone; % autumn; % cool winter bone
            rgb = interp1(linspace(0, 1, 256), cmap, I_norm);
            obj.bw_map = bw_map;
            obj.I_norm = I_norm;
            obj.dist_transform = dist_transform;
            obj.D1 = D1;
            obj.rgb = rgb;
            obj.resolution = 0.01; % 分辨率是 1 pixel = 1 cm
            obj.RATION = 100;
            [obj.cols,obj.rows,~]=size(dist_transform);
        end
        
        function showSDFMap(obj,fp)
            figure(fp);
            imshow(~obj.bw_map);
            hold on;
            im = imshow(obj.rgb.*repmat(~obj.bw_map,[1 1 3]));
            im.AlphaData = 0.2;
            imcontour(obj.D1);
%             title('Euclidean SDF Map');
            axis on
            grid on
            axis equal
            set(gca,'Ydir','normal');
        end

        function [dist,xgrad,ygrad] = getDistAndGrad(obj,xpos,ypos)
            %getDistAndGrad 获取[xpos,ypos] 点的SDF信息和Grad信息
            xpos = xpos * obj.RATION; ypos = ypos * obj.RATION;
            xpos = floor(xpos); ypos = floor(ypos);
            % 这里比较地图大小的BUG x,y好像反了
            if (xpos<=0) 
                xpos=2;end
            if (ypos<=0) 
                ypos=2;end
            if (xpos+1>=obj.rows) 
%                 fprintf("xpos=%d\n",xpos);
                xpos=obj.rows-2;end
            if (ypos+1>=obj.cols) 
%                 fprintf("ypos%d\n",ypos);
                ypos=obj.cols-2;end
%             fprintf("xpos=%d,ypos=%d\n",xpos,ypos);
            distx=zeros(2,1);disty=zeros(2,1);
            if(isnan(xpos) || isnan(ypos))
                dist = 0;
            else
                dist = obj.dist_transform(ypos+1,xpos+1);
                distx(1) = obj.dist_transform(ypos+1,xpos);
                distx(2) = obj.dist_transform(ypos+1,xpos+2);
                disty(1) = obj.dist_transform(ypos,xpos+1);
                disty(2) = obj.dist_transform(ypos+2,xpos+1);
            end
%             if (dist) < 0.0001
%                 dist = 99999;
%             end
            %%%% fix grad bug
            xgrad = sum((dist-distx(1) + distx(2)-dist))/2.0;
            ygrad = sum((dist-disty(1) + disty(2)-dist))/2.0;
            % BUG 好像这里求梯度还得除以分辨率?
            % BUG 这里的梯度方向好像不对呀??? 是不是应该顺序依次相减?
            % 这个距离场的梯度绝对是错的
            dist = dist * obj.resolution; % 分辨率转化
        end
    end
end