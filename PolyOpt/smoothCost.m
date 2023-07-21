function [cost,grad] = smoothCost(coeffs,segpoly)
%SMOOTHCOST 此处显示有关此函数的摘要
% equation
% \min_{s(t),T} \mathcal{J} _{\mathrm{s}}=\int_0^{T_{\mathrm{j}}}{u_{\mathrm{j}}}(t)^{\mathrm{T}_{\mathrm{j}}}Ru_{\mathrm{j}}(t)dt+\rho T
% \\
% s^{(k)}(t)=u(t),s^{(k)}(t)=\sum_{i=k}^n{c_{\mathrm{ji}}t^i}
% \\
% \frac{\partial \mathcal{J} _{\mathrm{s}}}{\partial c_i}=2\left( \int_0^{T_i}{\beta}(t)\beta (t)^{\mathrm{T}}\mathrm{d}t \right) c_i
% \\
% \frac{\partial \mathcal{J} _{\mathrm{s}}}{\partial T_i}=c_{i}^{\mathrm{T}}\beta \left( T_i \right) \beta \left( T_i \right) ^{\mathrm{T}}c_i+\rho 
% 推导公式见：HUST_MS\1-Aggressive Navigation\documents\TOequation0.nb
% SmoothCostatT =  576 c4^2 T + 2880 c4 c5 T^2 + 4800 c5^2 T^3 + 5760 c4 c6 T^3 + 
%                  21600 c5 c6 T^4 + 10080 c4 c7 T^4 + 25920 c6^2 T^5 + 
%                  40320 c5 c7 T^5 + 100800 c6 c7 T^6 + 100800 c7^2 T^7
% SmoothCostgradcatT = {0, 0, 0, 0, 1152 c4 t + 2880 c5 t^2 + 5760 c6 t^3 + 10080 c7 t^4, 
%                      2880 c4 t^2 + 9600 c5 t^3 + 21600 c6 t^4 + 40320 c7 t^5, 
%                      5760 c4 t^3 + 21600 c5 t^4 + 51840 c6 t^5 + 100800 c7 t^6, 
%                      10080 c4 t^4 + 40320 c5 t^5 + 100800 c6 t^6 + 201600 c7 t^7}
% SmoothCostgradt =    576 c4^2 + 5760 c4 c5 T + 14400 c5^2 T^2 + 17280 c4 c6 T^2 + 
%                      86400 c5 c6 T^3 + 40320 c4 c7 T^3 + 129600 c6^2 T^4 + 
%                      201600 c5 c7 T^4 + 604800 c6 c7 T^5 + 705600 c7^2 T^6
n_seg   = segpoly.seg;
n_order = segpoly.norder;
n_dim   = segpoly.Dim;

% 时间最优选项
TimeOptimal = segpoly.TimeOptimal;

% 等式约束降维选项
ReduceOptimalValue = segpoly.ReduceOptimalValue;

% //Debug: 保存数据
if isfield(segpoly,'SAVE_DATA')
    SAVE_DATA = segpoly.SAVE_DATA;
else
    SAVE_DATA = false;
end

if (TimeOptimal)
    ts = coeffs(end-n_seg + 1:end); % 最优的n_seg个变量是tau
    coeffs(end-n_seg + 1:end) = [];
    ts = exp(ts);
else
    ts = segpoly.T;
end
% disp("Iter ts datas:");
% disp(ts');

if (ReduceOptimalValue)
    segpoly.T = ts;
    coeffs = segpoly.traj.SolveCoeffs(coeffs,segpoly);
end

gdsmo = [];
gdsmocoeff = [];

cost = 0;
grad = zeros(1,n_order*n_dim*n_seg);
gradt = zeros(1,n_seg);
for j = 1: n_seg
    % 取出 nodr长度的系数
    xc = coeffs((j-1)*n_order*n_dim + 1:(j-1)*n_order*n_dim + n_order);
    xcost = getSegPolyCost(xc,ts(j));
    xgrad = getSegPolyGardatc(xc,ts(j));
    yc = coeffs((j-1)*n_order*n_dim + n_order + 1:(j-1)*n_order*n_dim + n_order*2);
    ycost = getSegPolyCost(yc,ts(j));
    ygrad = getSegPolyGardatc(yc,ts(j));
    qc = coeffs((j-1)*n_order*n_dim + n_order*2 + 1:(j-1)*n_order*n_dim + n_order*3);
    qcost = getSegPolyCost(qc,ts(j));
    qgrad = getSegPolyGardatc(qc,ts(j));
    cost = cost + xcost + ycost + segpoly.R * qcost;
    grad((j-1)*n_order*n_dim + 1:(j-1)*n_order*n_dim + n_order*3) = [xgrad,ygrad,segpoly.R * qgrad];
    xgradt = getSegPolyGardatT(xc,ts(j));
    ygradt = getSegPolyGardatT(yc,ts(j));
    qgradt = getSegPolyGardatT(qc,ts(j));
    gradt(j) = (xgradt + ygradt + segpoly.R * qgradt)*ts(j); % 对数的导数还是自身
    
    gdsmocoeff = [gdsmocoeff;xc';yc';qc'];
    
    t = ts(j);
    c=xc;
    c(1)=[];
    xcgd = [576*c(4)^2 , 5760*c(4)*c(5)*t , 14400*c(5)^2*t^2 , 17280*c(4)*c(6)*t^2 , ...
         86400*c(5)*c(6)*t^3 , 40320*c(4)*c(7)*t^3 , 129600*c(6)^2*t^4 , ...
         201600*c(5)*c(7)*t^4 , 604800*c(6)*c(7)*t^5 , 705600*c(7)^2*t^6];
    c=yc;
    c(1)=[];
   ycgd = [576*c(4)^2 , 5760*c(4)*c(5)*t , 14400*c(5)^2*t^2 , 17280*c(4)*c(6)*t^2 , ...
         86400*c(5)*c(6)*t^3 , 40320*c(4)*c(7)*t^3 , 129600*c(6)^2*t^4 , ...
         201600*c(5)*c(7)*t^4 , 604800*c(6)*c(7)*t^5 , 705600*c(7)^2*t^6];
    c=qc;
    c(1)=[];
    qcgd = [576*c(4)^2 , 5760*c(4)*c(5)*t , 14400*c(5)^2*t^2 , 17280*c(4)*c(6)*t^2 , ...
         86400*c(5)*c(6)*t^3 , 40320*c(4)*c(7)*t^3 , 129600*c(6)^2*t^4 , ...
         201600*c(5)*c(7)*t^4 , 604800*c(6)*c(7)*t^5 , 705600*c(7)^2*t^6];
    xyqgd = [xcgd,ycgd,qcgd];
    gdsmo =[gdsmo, xyqgd'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% save debug data
if (SAVE_DATA)
    filename = "E:\datas\Swift\Debug\MATLABgdsmo.csv";
    TRAJ_DATA = gdsmo;
    TRAJ_DATA = round(TRAJ_DATA,8);
    writematrix(TRAJ_DATA,filename);
    filename = "E:\datas\Swift\Debug\MATLABgdsmocoeff.csv";
    TRAJ_DATA = gdsmocoeff;
    TRAJ_DATA = round(TRAJ_DATA,8);
    writematrix(TRAJ_DATA,filename);
    fprintf("save gdsmo debug datas, 上天保佑 \n");
end

% for j = 1: n_seg
%     % 取出 nodr长度的系数
%     xc = coeffs((j-1)*n_order*n_dim + 1:(j-1)*n_order*n_dim + n_order);
%     xgrad = getSegPolyGardatc(xc,ts(j));
%     yc = coeffs((j-1)*n_order*n_dim + n_order + 1:(j-1)*n_order*n_dim + n_order*2);
%     ygrad = getSegPolyGardatc(yc,ts(j));
%     qc = coeffs((j-1)*n_order*n_dim + n_order*2 + 1:(j-1)*n_order*n_dim + n_order*3);
%     qgrad = getSegPolyGardatc(qc,ts(j));
%     grad((j-1)*n_order*n_dim + 1:(j-1)*n_order*n_dim + n_order*3) = [xgrad,ygrad,segpoly.R * qgrad];
%     xgradt = getSegPolyGardatT(xc,ts(j));
%     ygradt = getSegPolyGardatT(yc,ts(j));
%     qgradt = getSegPolyGardatT(qc,ts(j));
%     gradt(j) = (xgradt + ygradt + segpoly.R * qgradt)*ts(j); % 对数的导数还是自身
% end

if(ReduceOptimalValue)
    grad = segpoly.traj.getReduceOptVelue(grad);
end
if(TimeOptimal)
    grad = [grad,gradt];
end
grad=grad';
end

%#########################################################################%
function cost = getSegPolyCost(c,t)
c(1)=[]; %MATLAB从1开始索引,所以c(7)应该是符号公式中的第8个系数
cost = 576 *c(4)^2*t + 2880*c(4)*c(5)*t^2 + 4800*c(5)^2*t^3 + 5760*c(4)*c(6)*t^3 + ... 
       21600*c(5)*c(6)*t^4 + 10080*c(4)*c(7)*t^4 + 25920*c(6)^2*t^5 + ... 
       40320*c(5)*c(7)*t^5 + 100800*c(6)*c(7)*t^6 + 100800*c(7)^2*t^7 ;
end
function gard = getSegPolyGardatc(c,t)
c(1)=[];
gard =  [0, 0, 0, 0, 1152*c(4)*t + 2880*c(5)*t^2 + 5760*c(6)*t^3 + 10080*c(7)*t^4,... 
         2880*c(4)*t^2 + 9600*c(5)*t^3 + 21600*c(6)*t^4 + 40320*c(7)*t^5, ...
         5760*c(4)*t^3 + 21600*c(5)*t^4 + 51840*c(6)*t^5 + 100800*c(7)*t^6, ...
         10080*c(4)*t^4 + 40320*c(5)*t^5 + 100800*c(6)*t^6 + 201600*c(7)*t^7];
end
function gard = getSegPolyGardatT(c,t)
c(1)=[];
gard =   576*c(4)^2 + 5760*c(4)*c(5)*t + 14400*c(5)^2*t^2 + 17280*c(4)*c(6)*t^2 + ...
         86400*c(5)*c(6)*t^3 + 40320*c(4)*c(7)*t^3 + 129600*c(6)^2*t^4 + ...
         201600*c(5)*c(7)*t^4 + 604800*c(6)*c(7)*t^5 + 705600*c(7)^2*t^6;
end

% %#########################################################################%
% %%%% acceleration cost
% %%%% getSegPolyCost
% 4 * c(2)^2 * t + 12 * c(2) * c(3) * t^2 + 12 * c(3)^2 * t^3 + 16 * c(2) * c(4) * t^3 + 36 * c(3) * c(4) * t^4 +
% 20 * c(2) * c(5) * t^4 + (144 * c(4)^2 * t^5) / 5 + 48 * c(3) * c(5) * t^5 + 24 * c(2) * c(6) * t^5 +
% 80 * c(4) * c(5) * t^6 + 60 * c(3) * c(6) * t^6 + 28 * c(2) * c(7) * t^6 + (400 * c(5)^2 * t^7) / 7 +
% 720/7 * c(4) * c(6) * t^7 + 72 * c(3) * c(7) * t^7 + 150 * c(5) * c(6) * t^8 + 126 * c(4) * c(7) * t^8 +
% 100 * c(6)^2 * t^9 + 560/3 * c(5) * c(7) * t^9 + 252 * c(6) * c(7) * t^10 + (1764 * c(7)^2 * t^11) / 11
% %%%% getSegPolyGardatc
% [
% 0,
% 0,
% 8 * c(2) * t + 12 * c(3) * t^2 + 16 * c(4) * t^3 + 20 * c(5) * t^4 + 24 * c(6) * t^5 + 28 * c(7) * t^6,
% 12 * c(2) * t^2 + 24 * c(3) * t^3 + 36 * c(4) * t^4 + 48 * c(5) * t^5 + 60 * c(6) * t^6 + 72 * c(7) * t^7,
% 16 * c(2) * t^3 + 36 * c(3) * t^4 + (288 * c(4) * t^5) / 5 + 80 * c(5) * t^6 + (720 * c(6) * t^7) / 7 + 126 * c(7) * t^8,
% 20 * c(2) * t^4 + 48 * c(3) * t^5 + 80 * c(4) * t^6 + (800 * c(5) * t^7) / 7 + 150 * c(6) * t^8 + (560 * c(7) * t^9) / 3,
% 24 * c(2) * t^5 + 60 * c(3) * t^6 + (720 * c(4) * t^7) / 7 + 150 * c(5) * t^8 + 200 * c(6) * t^9 + 252 * c(7) * t^10,
% 28 * c(2) * t^6 + 72 * c(3) * t^7 + 126 * c(4) * t^8 + (560 * c(5) * t^9) / 3 + 252 * c(6) * t^10 + (3528 * c(7) * t^11) / 11
% ]
% %%%% getSegPolyGardatT
% 4 * c(2)^2 + 24 * c(2) * c(3) * t + 36 * c(3)^2 * t^2 + 48 * c(2) * c(4) * t^2 + 144 * c(3) * c(4) * t^3 +
% 80 * c(2) * c(5) * t^3 + 144 * c(4)^2 * t^4 + 240 * c(3) * c(5) * t^4 + 120 * c(2) * c(6) * t^4 +
% 480 * c(4) * c(5) * t^5 + 360 * c(3) * c(6) * t^5 + 168 * c(2) * c(7) * t^5 + 400 * c(5)^2 * t^6 +
% 720 * c(4) * c(6) * t^6 + 504 * c(3) * c(7) * t^6 + 1200 * c(5) * c(6) * t^7 + 1008 * c(4) * c(7) * t^7 +
% 900 * c(6)^2 * t^8 + 1680 * c(5) * c(7) * t^8 + 2520 * c(6) * c(7) * t^9 + 1764 * c(7)^2 * t^10
% 
% %#########################################################################%
% %%%% jerk cost
% %%%% getSegPolyCost
% 36 * c(3)^2 * t + 144 * c(3) * c(4) * t^2 + 192 * c(4)^2 * t^3 + 240 * c(3) * c(5) * t^3 +
% 720 * c(4) * c(5) * t^4 + 360 * c(3) * c(6) * t^4 + 720 * c(5)^2 * t^5 + 1152 * c(4) * c(6) * t^5 +
% 504 * c(3) * c(7) * t^5 + 2400 * c(5) * c(6) * t^6 + 1680* c(4) * c(7) * t^6 +
% (14400 * c(6)^2 * t^7) / 7 + 3600 * c(5) * c(7) * t^7 + 6300 * c(6) * c(7) * t^8 + 4900 * c(7)^2 * t^9
% %%%% getSegPolyGardatc
% [
% 0,
% 0,
% 0,
% 72 * c(3) * t + 144 * c(4) * t^2 + 240 * c(5) * t^3 + 360 * c(6) * t^4 + 504 * c(7) * t^5,
% 144 * c(3) * t^2 + 384 * c(4) * t^3 + 720 * c(5) * t^4 + 1152 * c(6) * t^5 + 1680 * c(7) * t^6,
% 240 * c(3) * t^3 + 720 * c(4) * t^4 + 1440 * c(5) * t^5 + 2400 * c(6) * t^6 + 3600 * c(7) * t^7,
% 360 * c(3) * t^4 + 1152 * c(4) * t^5 + 2400 * c(5) * t^6 + (28800 * c(6) * t^7) / 7 + 6300 * c(7) * t^8,
% 504 * c(3) * t^5 + 1680 * c(4) * t^6 + 3600 * c(5) * t^7 + 6300 * c(6) * t^8 + 9800 * c(7) * t^9
% ]
% %%%% getSegPolyGardatT
% 36 * c(3)^2 + 288 * c(3) * c(4) * t + 576 * c(4)^2 * t^2 + 720 * c(3) * c(5) * t^2 +
% 2880 * c(4) * c(5) * t^3 + 1440 * c(3) * c(6) * t^3 + 3600 * c(5)^2 * t^4 + 5760 * c(4) * c(6) * t^4 +
% 2520 * c(3) * c(7) * t^4 + 14400 * c(5) * c(6) * t^5 + 10080 * c(4) * c(7) * t^5 +
% 14400 * c(6)^2 * t^6 + 25200 * c(5) * c(7) * t^6 + 50400 * c(6) * c(7) * t^7 + 44100 * c(7)^2 * t^8
