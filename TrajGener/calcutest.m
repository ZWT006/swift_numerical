close all;clear;clc;


n_inputorder = 4;
n_order = 6;
n_seg = 2;
ts=ones(n_seg,1);
ts(1)=2;
start_cond=ones(n_inputorder,1);

n_all_poly = n_seg*(n_order+1);
coeff_n = n_order+1;

Aeq_start = zeros(n_inputorder, n_all_poly);
% STEP 2.1: write expression of Aeq_start and beq_start
polyCoeff =  getCoeffCons(0,n_order,n_inputorder);
bcoeff=BerneteinCoeff(n_order);
bezierCoeff = polyCoeff*bcoeff;
t = ts(1);

% scal_t=ones(n_inputorder,1);
% for idx=1:n_inputorder
% scal_t(idx) = t^(idx-1);
% end



bezierCoeff_m = bezierCoeff./scal_t;
Aeq_start(1:n_inputorder,1:coeff_n) = bezierCoeff_m;
beq_start =  start_cond';