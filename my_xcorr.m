function [rho, lags] = my_xcorr(x, t, max_lag)
% 输入：
% x - 一维时序数据
% t - 时序信号对应的时间
% max_lag - 最大滞后时间步数

    % 检查输入序列均值是否为0，不为零则处理：
    mean_value = mean(x);
    if mean_value > max(abs(x)) * 1e-6
        x = x - mean_value;
    end
    
    len = length(x);
    
    max_lag = min(max_lag, len-1); % 最大延迟不超过数据序列长度
    if max_lag == 0
        max_lag = len-1;
    end
    
    variance = mean(x .* x); % 计算方差（均值形式）
    total_variance = variance * len; % 计算总平方和

    rho = zeros(1, max_lag+1);
    lags = zeros(1, max_lag+1);
    rho(1) = 1; % 滞后0的自相关为1
    for i = 1:max_lag
        % 计算分子：求和而不是求平均
        temp = sum(x(1:len-i) .* x(1+i:end));
        rho(i+1) = temp / total_variance; % 除以总平方和
        lags(i+1) = t(i+1) - t(1);
    end
end