function time_scale = calculate_timescale(all_data, all_time, method)
% 计算粒子速度时间尺度
% 输入:
%   all_data - 粒子数据，其中第1列为流向分量
%   all_time - 时间序列
%   method - 方法类型: 'fit' 或 'integral'
% 输出:
%   time_scale - 计算得到的时间尺度

    % 提取流向分量并去中心化
    x = all_data(:,1); % 流向分量
    x = x - mean(x);
    t = all_time;
    
    % 计算自相关函数
    [rho, lags] = my_xcorr(x, t, 00000);
    
    % 根据指定方法计算时间尺度
    switch method
        case 'fit'
            % 方法1: 拟合自相关函数
            [~, t_l, ~, ~, ~] = lesutils.fit_exp_model(lags, rho, 0);
            time_scale = -1 / t_l; % 拟合系数为负倒数
            
        case 'integral'
            % 方法2: 积分到自相关函数过零点
            [~, ~, zero_cross_timescale, ~] = integrate_from_start(lags, rho);
            time_scale = zero_cross_timescale;
            
        otherwise
            error('不支持的method类型。请使用 ''fit'' 或 ''integral''');
    end
end
