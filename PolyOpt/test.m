%% C++ 代码快速计算函数对比测试
LOOPNUM = 80000;

% 定义范围
t = linspace(-10, 2, 1000);  % 曲线的范围，可以根据需要进行调整

% 计算 expC2 和 MATLAB 的 exp 对应点的值
expC2Start = tic;
for i=1:LOOPNUM
    expC2_values = expC2(t);
end
expC2End = toc(expC2Start);
expStart = tic;
for i=1:LOOPNUM
    exp_values = exp(t);
end
expEnd = toc(expStart);

% 计算 logC2 和 MATLAB 的 log 对应点的值
logC2Start = tic;
for i=1:LOOPNUM
    logC2_values = logC2(t);
end
logC2End = toc(logC2Start);
logStart = tic;
for i=1:LOOPNUM
    log_values = log(t);
end
logEnd = toc(logStart);

fprintf("expC2 time %f, exp time %f \n",expC2End,expEnd);
fprintf("logC2 time %f, log time %f \n",logC2End,logEnd);

% 绘制对比曲线
figure;
subplot(2, 1, 1);
plot(t, expC2_values, 'r', t, exp_values, 'b');
title('expC2 vs exp');
legend('expC2', 'exp');

subplot(2, 1, 2);
plot(t, logC2_values, 'r', t, log_values, 'b');
title('logC2 vs log');
legend('logC2', 'log');



%% BandedSystem
lowerBw = 10;
upperBw = 8;
Mat = BandedSystem(100,lowerBw,upperBw);
for i=1:100
    for j=1:100
        if ((i-j)<= lowerBw && (j-i)<= upperBw)
            Mat(i,j) = i+j;
        end
    end
end

fprintf("Mat(i,j) = %2.2f \n",Mat(99,98));
fprintf("Mat(i,j) = %2.2f \n",Mat(70,60));
fprintf("Mat(i,j) = %2.2f \n",Mat(75,70));
fprintf("Mat(i,j) = %2.2f \n",Mat(2,3));
fprintf("Mat(i,j) = %2.2f \n",Mat(80,3));

function result = expC2(t)
    if t > 0.0
        result = ((0.5 .* t + 1.0) .* t + 1.0);
    else
        result = 1.0 ./ ((0.5 .* t - 1.0) .* t + 1.0);
    end
end

function result = logC2(T)
    if T > 1.0
        result = sqrt(2.0 .* T - 1.0) - 1.0;
    else
        result = 1.0 - sqrt(2.0 ./ T - 1.0);
    end
end