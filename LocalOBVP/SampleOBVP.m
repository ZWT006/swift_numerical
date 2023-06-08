classdef SampleOBVP
    %SAMPLEOBVP 采样阶段的快速OBVP计算 解析形式
    %   此处显示详细说明
    
    properties
        vel_factor; % 线速度因子;用于计算参考时间
        w_factor;   % 角速度因子
        refer_T;    % 参考时间
        xc;yc;qc;   % polynomial coefficients for [x,y,q]
        spi;svi;sai;spf;
        costJ;
        R;          % [x,y,q]的权重系数
        dt;         % 离散时间
        path_length;% path length integration
    end
    
    methods
        %% construction function
        function obj = SampleOBVP(Vel_factor,W_factor,R,dt)
            %SAMPLEOBVP 定义比例因子和cost_function的权重
            %   
            obj.vel_factor = Vel_factor;
            obj.w_factor   = W_factor;
            obj.R = R;
            obj.dt = dt;
        end
        %% SolveBVP
        function [obj,xt,yt,qt,svf] = SolveBVP(obj,spi,svi,spf,svf)
            %Solve 计算两点之间的轨迹 fix start;final
            % 1)计算轨迹点的多项式系数;
            % 2)轨迹参考时间
            % 3)轨迹的综合cost
            idx=1;idy=2;idq=3; %state的索引
            % step1: 角度处理
            deltaq = AngleDelta(spi(idq),spf(idq));
            spf(idq) = spi(idq) + deltaq;
            % step2: 求参考时间,线速度时间和角速度时间取较大值
            T = max([sqrt((spi(idx)-spf(idx))^2+(spi(idy)-spf(idy))^2)/obj.vel_factor,abs(deltaq)/obj.w_factor]);
            % step3: obvp 计算
%             si=[spi(idx),svi(idx),0];sf=[spf(idx),svf(idx),0];
            si=[spi(idx),svi(idx),0];sf=[spf(idx),svf(idx),0];
            c0 = si(1);
            c1 = si(2);
            c2 = si(3)/2;
            c3 = -((-20*sf(1) + 20*si(1) - sf(3)*T^2 + 3*si(3)*T^2 + 8*T*sf(2) + 12*T*si(2))/(2*T^3));
            c4 = -((30*sf(1) - 30*si(1) + 2*sf(3)*T^2 - 3*si(3)*T^2 - 14*T*sf(2) - 16*T*si(2))/(2*T^4));
            c5 = -((-12*sf(1) + 12*si(1) - sf(3)*T^2 + si(3)*T^2 + 6*T*sf(2) + 6*T*si(2))/(2*T^5));
            obj.xc=[c5 c4 c3 c2 c1 c0];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            si=[spi(idy),svi(idy),0];sf=[spf(idy),svf(idy),0];
            c0 = si(1);
            c1 = si(2);
            c2 = si(3)/2;
            c3 = -((-20*sf(1) + 20*si(1) - sf(3)*T^2 + 3*si(3)*T^2 + 8*T*sf(2) + 12*T*si(2))/(2*T^3));
            c4 = -((30*sf(1) - 30*si(1) + 2*sf(3)*T^2 - 3*si(3)*T^2 - 14*T*sf(2) - 16*T*si(2))/(2*T^4));
            c5 = -((-12*sf(1) + 12*si(1) - sf(3)*T^2 + si(3)*T^2 + 6*T*sf(2) + 6*T*si(2))/(2*T^5));
            obj.yc=[c5 c4 c3 c2 c1 c0];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            si=[spi(idq),svi(idq),0];sf=[spf(idq),svf(idq),0];
            c0 = si(1);
            c1 = si(2);
            c2 = si(3)/2;
            c3 = -((-20*sf(1) + 20*si(1) - sf(3)*T^2 + 3*si(3)*T^2 + 8*T*sf(2) + 12*T*si(2))/(2*T^3));
            c4 = -((30*sf(1) - 30*si(1) + 2*sf(3)*T^2 - 3*si(3)*T^2 - 14*T*sf(2) - 16*T*si(2))/(2*T^4));
            c5 = -((-12*sf(1) + 12*si(1) - sf(3)*T^2 + si(3)*T^2 + 6*T*sf(2) + 6*T*si(2))/(2*T^5));
            obj.qc=[c5 c4 c3 c2 c1 c0];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.costJ=0;
            obj.refer_T = T;
            xt=f_st(obj.xc,obj.dt,T);
            yt=f_st(obj.yc,obj.dt,T);
            qt=f_st(obj.qc,obj.dt,T);
            % 返回svf 分清节点
            svf = [xt(end,2),yt(end,2),qt(end,2)];
            % path length integration
            obj.path_length = path_ds(obj.xc,obj.yc,obj.dt,T);
            % 
            obj.spf = spf;
            obj.spi = spi;
            obj.svi = svi;
        end
        %% SolveMiniAccInputAcc
        function [obj,xt,yt,qt,svf] = SolveMiniAccInputAcc(obj,spi,svi,spf)
            %Solve 计算两点之间的轨迹 Mini Acc; Input Acc
            % 1)计算轨迹点的多项式系数;
            % 2)轨迹参考时间
            % 3)轨迹的综合cost
            idx=1;idy=2;idq=3; %state的索引
            % step1: 角度处理
            deltaq = AngleDelta(spi(idq),spf(idq));
            spf(idq) = spi(idq) + deltaq;
            % step2: 求参考时间,线速度时间和角速度时间取较大值
            T = max([sqrt((spi(idx)-spf(idx))^2+(spi(idy)-spf(idy))^2)/obj.vel_factor,abs(deltaq)/obj.w_factor]);
            % step3: obvp 计算
            alphax=3*(svi(idx)*T+spi(idx)-spf(idx))/T^3;
            alphay=3*(svi(idy)*T+spi(idy)-spf(idy))/T^3;
            % 对这里的*obj.R抱有疑问,不确定是不是会影响结果准确性
            alphaq=3*(svi(idq)*T+spi(idq)-spf(idq))/T^3*obj.R;
            obj.costJ=alphax/3*T^3+alphay/3*T^3+alphaq/3*T^3*obj.R;
            obj.xc=[alphax/6;-alphax/2*T;svi(idx);spi(idx)];
            obj.yc=[alphay/6;-alphay/2*T;svi(idy);spi(idy)];
            obj.qc=[alphaq/6;-alphaq/2*T;svi(idq);spi(idq)];
            obj.refer_T = T;
            xt=f_st(obj.xc,obj.dt,T);
            yt=f_st(obj.yc,obj.dt,T);
            qt=f_st(obj.qc,obj.dt,T);
            % 返回svf 分清节点
            svf = [xt(end,2),yt(end,2),qt(end,2)];
            % path length integration
            obj.path_length = path_ds(obj.xc,obj.yc,obj.dt,T);
            % 
            obj.spf = spf;
            obj.spi = spi;
            obj.svi = svi;
        end
        %%  
        function [obj,xt,yt,qt,svf,saf] = SolveMiniJerkInputJerk(obj,spi,svi,sai,spf)
            %Solve 计算两点之间的轨迹 Mini Acc; Input Acc
            % 1)计算轨迹点的多项式系数;
            % 2)轨迹参考时间
            % 3)轨迹的综合cost
            idx=1;idy=2;idq=3; %state的索引
            % step1: 角度处理
            deltaq = AngleDelta(spi(idq),spf(idq));
            spf(idq) = spi(idq) + deltaq;
            % step2: 求参考时间,线速度时间和角速度时间取较大值
            T = max([sqrt((spi(idx)-spf(idx))^2+(spi(idy)-spf(idy))^2)/obj.vel_factor,abs(deltaq)/obj.w_factor]);
            % step3: obvp 计算
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 20*(pf-pi-vi*temp_T-ai*temp_T^2/2)/temp_T^5
            % calculate alpha
            alphax=20*(-sai(idx)*T^2/2-svi(idx)*T-spi(idx)+spf(idx))/T^5;
            alphay=20*(-sai(idy)*T^2/2-svi(idy)*T-spi(idy)+spf(idy))/T^5;
            alphaq=20*(-sai(idq)*T^2/2-svi(idq)*T-spi(idq)+spf(idq))/T^5;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (T^7*alpha^2)/56 + (T^4*alpha*ao)/4 + T*ao^2
            % calculate cost 利用a的多项式积分进行计算
            costax = (T^7*alphax^2)/56 + (T^4*alphax*sai(idx))/4 + T*sai(idx)^2;
            costay = (T^7*alphay^2)/56 + (T^4*alphay*sai(idy))/4 + T*sai(idy)^2;
            costaq = ((T^7*alphaq^2)/56 + (T^4*alphaq*sai(idq))/4 + T*sai(idq)^2)*obj.R;
            obj.costJ = costax + costay + costaq;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % alpha/120 -alpha/24*obj.T alpha/12*obj.T^2 si(3)/2 si(2) si(1);
            %  calculate polynomial coefficients
            obj.xc=[alphax/120;-alphax/24*T;alphax/12*T^2;sai(idx)/2;svi(idx);spi(idx)];
            obj.yc=[alphay/120;-alphay/24*T;alphay/12*T^2;sai(idy)/2;svi(idy);spi(idy)];
            obj.qc=[alphaq/120;-alphaq/24*T;alphaq/12*T^2;sai(idq)/2;svi(idq);spi(idq)];
            obj.refer_T = T;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % polynomial calculate
            xt=f_st(obj.xc,obj.dt,T);
            yt=f_st(obj.yc,obj.dt,T);
            qt=f_st(obj.qc,obj.dt,T);
            % 返回svf 分清节点
            svf = [xt(end,2),yt(end,2),qt(end,2)];
            saf = [xt(end,3),yt(end,3),qt(end,3)];
            
            % path length integration
            obj.path_length = path_ds(obj.xc,obj.yc,obj.dt,T);

            % 
            obj.spf = spf;
            obj.spi = spi;
            obj.svi = svi;
            obj.sai = sai;
        end
        %%
        function [xt,yt,qt] = st(obj)
            % polynomial calculate
            xt=f_st(obj.xc,obj.dt,obj.refer_T);
            yt=f_st(obj.yc,obj.dt,obj.refer_T);
            qt=f_st(obj.qc,obj.dt,obj.refer_T);
        end
    end
end

function st = f_st(psc,dt,T)
    t=0:dt:T+dt;
    st(:,1) = polyval(psc,t);
    vsc=polyder(psc);
    st(:,2) = polyval(vsc,t);
    asc=polyder(vsc);
    st(:,3) = polyval(asc,t);
    jsc=polyder(asc);
    st(:,4) = polyval(jsc,t);
end

function path_length = path_ds(pxc,pyc,ds,T)
    t=0:ds:T+ds;
    x(:,1) = polyval(pxc,t);
    y(:,1) = polyval(pyc,t);
    path_length = sum(sqrt(diff(x).^2 + diff(y).^2));
%     path_length = 0;
%     for idx=1:length(t)-1
%         path_length = path_length + sqrt((x(idx)-x(idx+1))^2+(y(idx)-y(idx+1))^2);
%     end
end
