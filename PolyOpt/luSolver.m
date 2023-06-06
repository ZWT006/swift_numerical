classdef luSolver
    %LUSOLVER 用于线性方程组的求解
    %   Ax=b
    
    properties
        Aeq; % A 矩阵
        beq; % b 向量
        lowerBw = 50;
        upperBw = 50;
        LU;
        x;
    end
    
    methods
        %##################################################################
        function obj = luSolver(Aeq,beq,lowerBw,upperBw)
            %LUSOLVER 
            obj.Aeq = Aeq;
            obj.beq = beq;
            obj.lowerBw = lowerBw;
            obj.upperBw = upperBw;
        end
        %##################################################################
        function obj = factorizeLU(obj)
            n = size(obj.Aeq,1);
            LU = obj.Aeq;
            
            for k = 1:n-1
                iM = min(k+obj.lowerBw,n);
                cVl = LU(k,k);
                for i = k+1:iM
                    if LU(i,k) ~= 0
                        LU(i,k) = LU(i,k)/cVl;
%                         LU(i,k+1:n) = LU(i,k+1:n) - LU(i,k)*LU(k,k+1:n);
%                         LU(i,k) = 0;
                    end
                end
                jM = min(k+obj.upperBw,n);
                for j = k+1:jM
                    cVl = LU(k,j);
                    if cVl ~= 0
                        for i = k+1:iM
                            if LU(i,k) ~= 0
                                LU(i,j) = LU(i,j) - LU(i,k)*cVl;
                            end
                        end
                    end
                end
            end
            obj.LU = LU;
        end
        %##################################################################
        function obj = Solve(obj)
            n = size(obj.Aeq,1);
            m = size(obj.beq,2);
            x = obj.beq;
            % 前向消元
            for j = 1:n-1
                iM = min(j+obj.lowerBw,n);
                for i = j+1:iM
                    if obj.LU(i,j) ~= 0
                        x(i,:) = x(i,:) - obj.LU(i,j)*x(j,:);
                    end
                end
            end
            
            % 后向代入
            for j = n:-1:1
                x(j,:) = x(j,:)/obj.LU(j,j);
                iM = max(1,j-obj.upperBw);
                for i = iM:j-1
                    if obj.LU(i,j) ~= 0
                        x(i,:) = x(i,:) - obj.LU(i,j)*x(j,:);
                    end
                end
            end
            obj.x = x;
        end
    end
end

