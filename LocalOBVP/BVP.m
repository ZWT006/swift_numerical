classdef BVP
    %BVP 求解多项式BVP问题
    %   根据终末条件求解多项式轨迹
    
    properties
        T; % 轨迹的周期
        c5;c4;c3;c2;c1;c0; %5th order polynomial trajectory
        % f(t) = c5*t^5 + c4*t^4 + c3*t^3 + c2*t^2 + c1*t +c0;
        % f'(t) = 5*c5*t^4 + 4*c4*t^3 + 3*c3*t^2 + 2*c2*t + c1;
        % f''(t) = 20*c5*t^3 + 12*c4*t^2 + 6*c3*t + 2*c2;
    end
    
    methods
        %%
        function obj = BVP()
            %BVP 构造此类啥都不需要
            %   
        end
        %%
        function obj = Solve(obj,si,sf,T)
            %Solve 求解多项式系数
            % 根据解析解公式直接求解
            % c5 -> -((-12 pf + 12 pi - af T^2 + ai T^2 + 6 T vf + 6 T vi)/(2 T^5)), 
            % c4 -> -((30 pf - 30 pi + 2 af T^2 - 3 ai T^2 - 14 T vf - 16 T vi)/(2 T^4)), 
            % c3 -> -((-20 pf + 20 pi - af T^2 + 3 ai T^2 + 8 T vf + 12 T vi)/(2 T^3)), 
            % c2 -> ai/2, 
            % c1 -> vi, 
            % c0 -> pi
            obj.T = T;
            obj.c0 = si(1);
            obj.c1 = si(2);
            obj.c2 = si(3)/2;
            obj.c3 = -((-20*sf(1) + 20*si(1) - sf(3)*T^2 + 3*si(3)*T^2 + 8*T*sf(2) + 12*T*si(2))/(2*T^3));
            obj.c4 = -((30*sf(1) - 30*si(1) + 2*sf(3)*T^2 - 3*si(3)*T^2 - 14*T*sf(2) - 16*T*si(2))/(2*T^4));
            obj.c5 = -((-12*sf(1) + 12*si(1) - sf(3)*T^2 + si(3)*T^2 + 6*T*sf(2) + 6*T*si(2))/(2*T^5));
        end
        %% st
        function [pst,vst,ast]= ft(obj,dt)
            t=0:dt:obj.T;
            pst =obj.c5*t.^5 +obj.c4*t.^4 +obj.c3*t.^3 +obj.c2*t.^2 +obj.c1*t + obj.c0;
            vst =5*obj.c5*t.^4 + 4*obj.c4*t.^3 + 3*obj.c3*t.^2 + 2*obj.c2*t + obj.c1 ;
            ast =20*obj.c5*t.^3 + 12*obj.c4*t.^2 + 6*obj.c3*t + 2*obj.c2 ;
        end
    end
end

