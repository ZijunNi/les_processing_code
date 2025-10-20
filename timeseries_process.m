% 本代码主要用于处理时序信息（包括Euler和Lagrange）
%% Lagrange时间尺度计算
close all
clear A B
[A,B] = filter_cell(particleData);

%%
check_data_quality(A);
%%
calculate_ensemble_timescale(A);


%%
[R_avg, tau, T_l, T_l_u, T_l_v, T_l_w, R_u, R_v, R_w] = lesutils.calculate_velocity_timescale( ...
    particleData{32*16+19}.data(:,1), ...
    particleData{32*16+19}.data(:,5:7));

zero_cross_idx_u = find(R_u < 0, 1);
if isempty(zero_cross_idx_u), zero_cross_idx_u = numel(R_u); end
zero_cross_idx_v = find(R_v < 0, 1);
if isempty(zero_cross_idx_v), zero_cross_idx_v = numel(R_v); end
zero_cross_idx_w = find(R_w < 0, 1);
if isempty(zero_cross_idx_w), zero_cross_idx_w = numel(R_w); end
zero_cross_idx_avg = find(R_avg < 0, 1);
if isempty(zero_cross_idx_avg), zero_cross_idx_avg = numel(R_avg); end

figure('Color', 'white');
plot(tau, R_u, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Streamwise (u)'); hold on;
plot(tau, R_v, 'g--', 'LineWidth', 1.2, 'DisplayName', 'Wall-normal (v)');
plot(tau, R_w, 'b--', 'LineWidth', 1.2, 'DisplayName', 'Spanwise (w)');
plot(tau, R_avg, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Average');
plot(tau(1:zero_cross_idx_avg), R_avg(1:zero_cross_idx_avg), 'm', 'LineWidth', 2.0, ...
     'DisplayName', 'Integral Range');
yline(0, 'k-');
xlabel('Time $\tau$', 'Interpreter', 'latex');
ylabel('Autocorrelation $R(\tau)$', 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex');
grid on;

fprintf('\n拉格朗日时间尺度结果:\n');
fprintf('u方向: T_L_u = %.4f (积分上限 τ_u = %.4f)\n', T_L_u, tau(zero_cross_idx_u));
fprintf('v方向: T_L_v = %.4f (积分上限 τ_v = %.4f)\n', T_L_v, tau(zero_cross_idx_v));
fprintf('w方向: T_L_w = %.4f (积分上限 τ_w = %.4f)\n', T_L_w, tau(zero_cross_idx_w));
fprintf('平均: T_L = %.4f (积分上限 τ = %.4f)\n', T_l, tau(zero_cross_idx_avg));


function [R_ensemble, tau, T_int_ensemble, T_int_u, T_int_v, T_int_w] = calculate_ensemble_timescale(all_data)
    % 计算多个采样点的集成积分时间尺度
    %
    % 输入:
    %   all_data - 元胞数组，每个元素是一个N×7矩阵 [时间, x, y, z, u, v, w]
    %   dt - 采样时间间隔
    %
    % 输出:
    %   R_ensemble - 集成平均自相关函数
    %   tau - 时间延迟向量
    %   T_int_ensemble - 集成平均积分时间尺度
    %   T_int_u - u方向的集成平均积分时间尺度
    %   T_int_v - v方向的集成平均积分时间尺度
    %   T_int_w - w方向的集成平均积分时间尺度



    nPoints = numel(all_data);
    nTime = size(all_data{1}, 1);
    maxlag = min([10000, floor(nTime / 2), nTime - 1]);
    numLags = maxlag + 1;

    all_R_u = zeros(nPoints, numLags);
    all_R_v = zeros(nPoints, numLags);
    all_R_w = zeros(nPoints, numLags);
    tau = [];

    fprintf('计算%d个采样点的自相关函数...\n', nPoints);
    for i = 1:nPoints
        data = all_data{i};
        time = data(:,1);
        
        % 提取速度分量
        u = data(:,5);
        v = data(:,6);
        w = data(:,7);
        
        % 计算速度脉动
        u_prime = u - mean(u);
        v_prime = v - mean(v);
        w_prime = w - mean(w);

        [R_u, lags] = my_xcorr(u_prime, time, maxlag);
        [R_v, ~   ] = my_xcorr(v_prime, time, maxlag);
        [R_w, ~   ] = my_xcorr(w_prime, time, maxlag);

        all_R_u(i, :) = R_u;
        all_R_v(i, :) = R_v;
        all_R_w(i, :) = R_w;

        if isempty(tau)
            tau = lags(:);
        end
        
        % 显示进度
        if mod(i, ceil(nPoints/10)) == 0
            fprintf('已完成 %.0f%%\n', (i/nPoints)*100);
        end
    end

    R_u_ensemble = mean(all_R_u, 1).';
    R_v_ensemble = mean(all_R_v, 1).';
    R_w_ensemble = mean(all_R_w, 1).';
    R_ensemble = (R_u_ensemble + R_v_ensemble + R_w_ensemble) / 3;

    % 计算各方向的集成积分时间尺度
    % u方向
    zero_cross_idx_u = find(R_u_ensemble < 0, 1);
    if isempty(zero_cross_idx_u)
        zero_cross_idx_u = length(R_u_ensemble);
        warning('u方向自相关函数未过零点，积分到最大延迟');
    end
    T_int_u = trapz(tau(1:zero_cross_idx_u), R_u_ensemble(1:zero_cross_idx_u));
    
    % v方向
    zero_cross_idx_v = find(R_v_ensemble < 0, 1);
    if isempty(zero_cross_idx_v)
        zero_cross_idx_v = length(R_v_ensemble);
        warning('v方向自相关函数未过零点，积分到最大延迟');
    end
    T_int_v = trapz(tau(1:zero_cross_idx_v), R_v_ensemble(1:zero_cross_idx_v));
    
    % w方向
    zero_cross_idx_w = find(R_w_ensemble < 0, 1);
    if isempty(zero_cross_idx_w)
        zero_cross_idx_w = length(R_w_ensemble);
        warning('w方向自相关函数未过零点，积分到最大延迟');
    end
    T_int_w = trapz(tau(1:zero_cross_idx_w), R_w_ensemble(1:zero_cross_idx_w));
    
    % 计算总体集成积分时间尺度
    zero_cross_idx = find(R_ensemble < 0, 1);
    if isempty(zero_cross_idx)
        zero_cross_idx = length(R_ensemble);
        warning('平均自相关函数未过零点，积分到最大延迟');
    end
    T_int_ensemble = trapz(tau(1:zero_cross_idx), R_ensemble(1:zero_cross_idx));
    
    % 可视化结果
    create_visualizations(tau, R_u_ensemble, R_v_ensemble, R_w_ensemble, R_ensemble, ...
                         zero_cross_idx_u, zero_cross_idx_v, zero_cross_idx_w, zero_cross_idx, ...
                         T_int_u, T_int_v, T_int_w, T_int_ensemble);
    
    % 输出结果
    fprintf('\n集成积分时间尺度结果 (基于%d个采样点):\n', nPoints);
    fprintf('u方向: T_int_u = %.4f (积分上限 τ_u = %.4f)\n', T_int_u, tau(zero_cross_idx_u));
    fprintf('v方向: T_int_v = %.4f (积分上限 τ_v = %.4f)\n', T_int_v, tau(zero_cross_idx_v));
    fprintf('w方向: T_int_w = %.4f (积分上限 τ_w = %.4f)\n', T_int_w, tau(zero_cross_idx_w));
    fprintf('总体平均: T_int_ensemble = %.4f (积分上限 τ = %.4f)\n', T_int_ensemble, tau(zero_cross_idx));
end

function create_visualizations(tau, R_u, R_v, R_w, R_avg, ...
                              idx_u, idx_v, idx_w, idx_avg, ...
                              T_u, T_v, T_w, T_avg)
    % 创建可视化图表
    
    % 自相关函数图
    figure('Position', [100, 100, 800, 600], 'Color', 'white')
    
    % 绘制各方向自相关函数
    plot(tau, R_u, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Streamwise (u)')
    hold on
    plot(tau, R_v, 'g--', 'LineWidth', 1.2, 'DisplayName', 'Wall-normal (v)')
    plot(tau, R_w, 'b--', 'LineWidth', 1.2, 'DisplayName', 'Spanwise (w)')
    plot(tau, R_avg, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Ensemble Average')
    
    % 标记积分区域
    plot(tau(1:idx_u), R_u(1:idx_u), 'r', 'LineWidth', 3, 'DisplayName', sprintf('u Integral (T=%.3f)', T_u))
    plot(tau(1:idx_v), R_v(1:idx_v), 'g', 'LineWidth', 3, 'DisplayName', sprintf('v Integral (T=%.3f)', T_v))
    plot(tau(1:idx_w), R_w(1:idx_w), 'b', 'LineWidth', 3, 'DisplayName', sprintf('w Integral (T=%.3f)', T_w))
    plot(tau(1:idx_avg), R_avg(1:idx_avg), 'm', 'LineWidth', 3, 'DisplayName', sprintf('Avg Integral (T=%.3f)', T_avg))
    
    % 添加零线
    yline(0, 'k-', 'LineWidth', 1)
    
    % 添加垂直线标记积分上限
    xline(tau(idx_u), 'r--', 'LineWidth', 1)
    xline(tau(idx_v), 'g--', 'LineWidth', 1)
    xline(tau(idx_w), 'b--', 'LineWidth', 1)
    xline(tau(idx_avg), 'm--', 'LineWidth', 2)
    
    % 设置图表属性
    xlabel('Time Lag $\tau$', 'FontSize', 14, 'Interpreter', 'latex')
    ylabel('Autocorrelation Function $R(\tau)$', 'FontSize', 14, 'Interpreter', 'latex')
    title('Ensemble Autocorrelation Functions', 'FontSize', 16)
    legend('Location', 'best', 'FontSize', 10)
    grid on
    set(gca, 'FontSize', 12)
    
    % 设置坐标轴范围
    xlim([0, min([tau(idx_u)*2, tau(idx_v)*2, tau(idx_w)*2, tau(idx_avg)*2, tau(end)])])
    
    % 添加统计信息文本框
    annotation('textbox', [0.15, 0.15, 0.3, 0.2], 'String', ...
        {sprintf('T_{int,u} = %.3f', T_u), ...
         sprintf('T_{int,v} = %.3f', T_v), ...
         sprintf('T_{int,w} = %.3f', T_w), ...
         sprintf('T_{int,avg} = %.3f', T_avg)}, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FontSize', 10);
end


function [A, B] = filter_cell(cellArray)
    % 定义特定值
    dataSpecificValue = 5e-4;  % data矩阵第一行第二个元素需要等于的值
    idSpecificValue = 0;   % id需要等于的值（不需要筛选则填0）
    
    A = {};
    B = [];
    
    for i = 1:length(cellArray)
        currentCell = cellArray{i};
        id = currentCell.id;
        data = currentCell.data;
        
        % 检查data条件：第一行第二个元素等于特定值
        % dataCondition = false;
        % if size(data, 1) >= 1 && size(data, 2) >= 2
            dataCondition = abs(data(1, 3) - dataSpecificValue)<1e-4;
        % end
        
        % 检查id条件：id等于特定值
        idCondition = (id == idSpecificValue);
        
        % 如果任一条件满足，则提取数据
        if dataCondition || idCondition
            A{end+1} = data;
            B(end+1) = id;
        end
    end
end

% 在计算自相关函数前，进行数据质量检查
function check_data_quality(all_data)
    figure('Position', [100, 100, 1200, 800]);
    
    for i = 1:min(9, length(all_data)) % 检查前9个点
        data = all_data{i};
        u = data(:,5);
        
        subplot(3, 3, i);
        
        % 检查平稳性
        plot(u);
        title(sprintf('Point %d: Velocity Time Series', i));
        xlabel('Time step');
        ylabel('Velocity');
        
        % 计算并显示基本统计量
        fprintf('Point %d: mean=%.3f, std=%.3f, skewness=%.3f, kurtosis=%.3f\n', ...
                i, mean(u), std(u), skewness(u), kurtosis(u));
    end
    
    % 检查所有点的统计特性分布
    all_means = cellfun(@(x) mean(x(:,5)), all_data);
    all_stds = cellfun(@(x) std(x(:,5)), all_data);
    
    figure;
    subplot(1, 2, 1);
    histogram(all_means);
    title('Distribution of Mean Velocity across Points');
    
    subplot(1, 2, 2);
    histogram(all_stds);
    title('Distribution of Velocity STD across Points');
end
