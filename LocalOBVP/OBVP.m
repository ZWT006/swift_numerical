classdef OBVP
    %OBVP 求解多项式OBVP问题
    %   
    
    properties
        T; % 轨迹的周期
        opt_T; %计算出来的最优时间
        c5;c4;c3;c2;c1;c0; %5th order polynomial trajectory
        alpha;beta;gamma;
        J; % cost function
        % f(t) = c5*t^5 + c4*t^4 + c3*t^3 + c2*t^2 + c1*t +c0;
        % f'(t) = 5*c5*t^4 + 4*c4*t^3 + 3*c3*t^2 + 2*c2*t + c1;
        % f''(t) = 20*c5*t^3 + 12*c4*t^2 + 6*c3*t + 2*c2;
    end
    
    methods
        %%
        function obj = OBVP()
            %OBVP 构造此类啥都不需要
            %   
        end
        %% 两端固定时间固定的Optimal Control求解
        function obj = FixSolve(obj,si,sf,T)
            %Solve 求解多项式系数 fix init; fix final; fix T
            % 根据解析解公式直接求解
            % delta_p = pf - pi - vi*T - 1/2*ai*T^2;
            % delta_v = vf - vi - ai*T
            % delta_a = af - ai;
            pi=si(1);vi=si(2);ai=si(3);
            pf=sf(1);vf=sf(2);af=sf(3);
            delta_p = pf - pi - vi*T - 1/2*ai*T^2;
            delta_v = vf - vi - ai*T;
            delta_a = af - ai;
            MA =[720,-360*T,60*T^2;
                -360*T,168*T^2,-24*T^3;
                60*T^2,-24*T^3,3*T^4];
            aby=1/T^5*MA*[delta_p;delta_v;delta_a];
            obj.alpha = aby(1);
            obj.beta  = aby(2);
            obj.gamma = aby(3);
            obj.J = obj.gamma^2 + obj.beta*obj.gamma*obj.T + ...
                1/3*(obj.beta^2+obj.alpha*obj.gamma)*obj.T^2 + ...
                1/4*obj.alpha*obj.beta*obj.T^3 + 1/20*obj.alpha^2*T^4;
            obj.T = T;
            obj.c0 = si(1);
            obj.c1 = si(2);
            obj.c2 = si(3)/2;
            obj.c3 = aby(3)/6;
            obj.c4 = aby(2)/24;
            obj.c5 = aby(1)/120;
        end
        %% init固定 final 任意 时间不定的Optimal Control求解
        function [obj,flag]= FreeFinalFreeTimeSolve(obj,si,sf)
            %Solve 求解多项式系数 fix init; free final; free T
            % 根据解析解公式直接求解
            flag = true;
            pi=si(1);vi=si(2);ai=si(3);
            pf=sf(1);
            pfi = pf - pi;
            ply_J1 = [-1/2*ai,-vi,pfi];
            ply_J2 = [ai,4*vi,-6*pfi];
%             solve_t = [roots(ply_J1);roots(ply_J2)];
            solve_t = roots(ply_J2);
            L_solve = length(solve_t);
            realt = zeros(4,1);
            ndx=1;
            for idx=1:L_solve
                if isreal(solve_t(idx))
                    if solve_t(idx)>0
                        realt(ndx)=solve_t(idx);
                        ndx=ndx+1;
                    end
                end
            end
            if (ndx==1)
                % 该方程没有根
                flag = false;
                return;
            else
                temp_alpha=zeros(ndx-1,1);
                temp_costJ=zeros(ndx-1,1);
                for idx=1:ndx-1
                    temp_T = realt(idx);
                    temp_alpha(idx) = 20*(pf-pi-vi*temp_T-ai*temp_T^2/2)/temp_T^5;
                    temp_costJ(idx) = temp_alpha(idx)^2*temp_T^4/20*temp_T;
                end
            end
            [obj.J,idx] = min(temp_costJ);

            obj.T = realt(idx);
            alpha = temp_alpha(idx);
            obj.c0 = si(1);
            obj.c1 = si(2);
            obj.c2 = si(3)/2;
            obj.c3 = alpha/12*obj.T^2;
            obj.c4 = -alpha/24*obj.T;
            obj.c5 = alpha/120;
        end
        %% pst,vst,ast,jst
        function [pst,vst,ast,jst]= ft(obj,dt)
            %返回
            t=0:dt:obj.T;
            pst =obj.c5*t.^5 +obj.c4*t.^4 +obj.c3*t.^3 +obj.c2*t.^2 +obj.c1*t + obj.c0;
            vst =5*obj.c5*t.^4 + 4*obj.c4*t.^3 + 3*obj.c3*t.^2 + 2*obj.c2*t + obj.c1 ;
            ast =20*obj.c5*t.^3 + 12*obj.c4*t.^2 + 6*obj.c3*t + 2*obj.c2 ;
            jst =60*obj.c5*t.^2 + 24*obj.c4*t + 6*obj.c3;
        end
    end
end

