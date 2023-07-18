classdef BandedSystem
    %BANDEDSYSTEM 此处显示有关此类的摘要
    %   此处显示详细说明
    
properties
        N
        lowerBw
        upperBw
        ptrData
end
    
methods
    function obj = BandedSystem(n, p, q)
        obj.N = n;
        obj.lowerBw = p;
        obj.upperBw = q;
        actualSize = obj.N * (obj.lowerBw + obj.upperBw + 1);
        obj.ptrData = zeros(1, actualSize);
    end
    
    function obj = delete(obj)
        % No need to implement explicit destructor in MATLAB
    end
    
    function obj = reset(obj)
        obj.ptrData = zeros(1, obj.N * (obj.lowerBw + obj.upperBw + 1));
    end
    
    function value = subsref(obj, s)
        % Overload MATLAB's subsref to enable accessing elements using parenthesis
        if strcmp(s(1).type, '()')
            i = s(1).subs{1};
            j = s(1).subs{2};
            value = obj.ptrData((i - j + obj.upperBw) * obj.N + j);
        else
            value = builtin('subsref', obj, s); % Return default value for other cases
        end
    end
    
    function obj = subsasgn(obj, s, value)
        % Overload MATLAB's subsasgn to enable assigning values using parenthesis
        if strcmp(s(1).type, '()')
            i = s(1).subs{1};
            j = s(1).subs{2};
            obj.ptrData((i - j + obj.upperBw) * obj.N + j) = value;
        end
    end

    function obj = factorizeLU(obj)
        for k = 0:(obj.N - 2)
            iM = min(k + obj.lowerBw, obj.N - 1);
            cVl = obj(k + 1, k + 1);
            for i = (k + 2):iM
                if obj(i, k + 1) ~= 0.0
                    obj(i, k + 1) = obj(i, k + 1) / cVl;
                end
            end
            jM = min(k + obj.upperBw, obj.N - 1);
            for j = (k + 2):jM
                cVl = obj(k + 1, j);
                if cVl ~= 0.0
                    for i = (k + 2):iM
                        if obj(i, k + 1) ~= 0.0
                            obj(i, j) = obj(i, j) - (obj(i, k + 1) * cVl);
                        end
                    end
                end
            end
        end
    end

    function b = solve(obj, b)
        for j = 0:(obj.N - 1)
            iM = min(j + obj.lowerBw, obj.N - 1);
            for i = (j + 1):iM
                if obj(i, j + 1) ~= 0.0
                    b(i) = b(i) - obj(i, j + 1) * b(j + 1);
                end
            end
        end
        for j = (obj.N - 1):-1:0
            b(j + 1) = b(j + 1) / obj(j + 1, j + 1);
            iM = max(0, j - obj.upperBw);
            for i = iM:(j - 1)
                if obj(i + 1, j + 1) ~= 0.0
                    b(i + 1) = b(i + 1) - obj(i + 1, j + 1) * b(j + 1);
                end
            end
        end
    end

    function b = solveAdj(obj, b)
        for j = 0:(obj.N - 1)
            iM = min(j + obj.upperBw, obj.N - 1);
            b(j + 1) = b(j + 1) / obj(j + 1, j + 1);
            for i = (j + 1):iM
                if obj(j + 1, i + 1) ~= 0.0
                    b(i + 1) = b(i + 1) - obj(j + 1, i + 1) * b(j + 1);
                end
            end
        end
        for j = (obj.N - 1):-1:0
            iM = max(0, j - obj.lowerBw);
            for i = iM:(j - 1)
                if obj(j + 1, i + 1) ~= 0.0
                    b(i + 1) = b(i + 1) - obj(j + 1, i + 1) * b(j + 1);
                end
            end
        end
    end
   end
end

