function [cost,grad] = CostFunc(coeffs,segpoly)
%COSTFUNC 计算Trajectory Optimization所有cost的CostFunction
%   此处显示详细说明

DEBUG_PRINT = segpoly.DEBUG_PRINT;
DEBUG_PLOT  = segpoly.DEBUG_PLOT;


[smoCost,smograd]=smoothCost(coeffs,segpoly);
[obsCost,obsgrad]=obstacleCost(coeffs,segpoly);
[dynCost,dyngrad]=dynamicCost(coeffs,segpoly);
[timCost,timgrad]=timeCost(coeffs,segpoly);
[ovaCost,ovagrad]=ovalCost(coeffs,segpoly);

lambda_smooth   = segpoly.lambda_smooth;
lambda_obstacle = segpoly.lambda_obstacle;
lambda_dynamic  = segpoly.lambda_dynamic;
lambda_time     = segpoly.lambda_time;
lambda_oval     = segpoly.lambda_oval;


if (DEBUG_PRINT)
    if(lambda_smooth ~= 0)
%         [smoCost,smograd]=smoothCost(coeffs,segpoly);
%         cost = lambda_smooth*smoCost;
%         grad = lambda_smooth*smograd;
        fprintf("smoCost = %8.6f; ",lambda_smooth*smoCost);
    end
    if(lambda_obstacle ~= 0)
%         [obsCost,obsgrad]=obstacleCost(coeffs,segpoly);
%         cost = lambda_obstacle*obsCost + cost;
%         grad = lambda_obstacle*obsgrad + grad;
        fprintf("obsCost = %8.6f; ",lambda_obstacle*obsCost);
    end
    if(lambda_dynamic ~= 0)
%         [dynCost,dyngrad]=dynamicCost(coeffs,segpoly);
%         cost = lambda_dynamic*dynCost + cost;
%         grad = lambda_dynamic*dyngrad + grad;
        fprintf("dynCost = %8.6f; ",lambda_dynamic*dynCost);
    end
    if(lambda_oval ~= 0)
%         [ovaCost,ovagrad]=ovalCost(coeffs,segpoly);
%         cost = lambda_oval*ovaCost + cost;
%         grad = lambda_oval*ovagrad + grad;
        fprintf("ovaCost = %8.6f; ",lambda_oval*ovaCost);
    end
    if(lambda_time ~= 0)
%         [timCost,timgrad]=timeCost(coeffs,segpoly);
%         cost = lambda_time*timCost + cost;
%         grad = lambda_time*timgrad + grad;
        fprintf("timCost = %8.6f; ",lambda_time*timCost);
    end
    fprintf("\n");
% else
%     [smoCost,smograd]=smoothCost(coeffs,segpoly);
%     [obsCost,obsgrad]=obstacleCost(coeffs,segpoly);
%     [dynCost,dyngrad]=dynamicCost(coeffs,segpoly);
%     [timCost,timgrad]=timeCost(coeffs,segpoly);
%     [ovaCost,ovagrad]=ovalCost(coeffs,segpoly);
%     cost = lambda_smooth*smoCost + lambda_obstacle*obsCost + lambda_dynamic*dynCost + lambda_time*timCost + lambda_oval*ovaCost;
%     grad = lambda_smooth*smograd + lambda_obstacle*obsgrad + lambda_dynamic*dyngrad + lambda_time*timgrad + lambda_oval*ovagrad;
end

% filename = "E:\datas\Swift\Optimization\NLoptCostGradNoTime.csv";
% TRAJ_DATA = [[smoCost;0;smograd], [obsCost;0;obsgrad], [dynCost;0;dyngrad], [timCost;0;timgrad], [ovaCost;0;ovagrad]];
% TRAJ_DATA = round(TRAJ_DATA,4);
% writematrix(TRAJ_DATA,filename);

if(DEBUG_PLOT)
global iter costArray;
iter = iter + 1;
costArray = [costArray;[lambda_smooth*smoCost,lambda_obstacle*obsCost,lambda_dynamic*dynCost,lambda_time*timCost,lambda_oval*ovaCost]];
end
% fprintf("smoCost = %12.6f; obsCost = %12.6f\n",lambda_smooth*smoCost,lambda_obstacle*obsCost);
% fprintf("smoCost = %8.6f; obsCost = %8.6f; dynCost = %8.6f\n",lambda_smooth*smoCost,lambda_obstacle*obsCost,lambda_dynamic*dynCost);

% cost = lambda_smooth*smoCost + lambda_obstacle*obsCost;
% grad = lambda_smooth*smograd + lambda_obstacle*obsgrad;

% cost = lambda_smooth*smoCost + lambda_time*timCost + lambda_oval*ovaCost;
% grad = lambda_smooth*smograd + lambda_time*timgrad + lambda_oval*ovagrad;

cost = lambda_smooth*smoCost + lambda_obstacle*obsCost + lambda_dynamic*dynCost + lambda_time*timCost + lambda_oval*ovaCost;
grad = lambda_smooth*smograd + lambda_obstacle*obsgrad + lambda_dynamic*dyngrad + lambda_time*timgrad + lambda_oval*ovagrad;

% cost = lambda_smooth*smoCost;
% grad = lambda_smooth*smograd;

% cost = lambda_obstacle*obsCost;
% grad = lambda_obstacle*obsgrad;

end