
load('./processed_data/fluid_lagrange_7_particles_2025-10-04_09-05-17.mat')% 采用二阶Lagrange插值方法
analyzeParticleGroupLagrangianTimescale(particleData, 21, 'analysis_results.mat');

load('./processed_data/fluid_lagrange_2025-10-18_09-13-34.mat')% 采用二阶B样条插值方法
analyzeParticleGroupLagrangianTimescale(particleData, 32, 'analysis_results.mat');

% % 批量分析多个粒子群
% groupNumbers = [1, 3, 5, 7];
% for i = 1:length(groupNumbers)
%     fprintf('\n\n分析粒子群 %d...\n', groupNumbers(i));
%     analyzeParticleGroupLagrangianTimescale(particleData, groupNumbers(i), 'batch_analysis.mat');
% end

function analyzeParticleGroupLagrangianTimescale(particleData, groupNumber, filename)
    % 分析粒子群的Lagrangian时间尺度
    %
    % 输入参数：
    %   particleData - 结构体数组，包含data和id字段
    %   groupNumber - 粒子群编号（0-63）
    %   filename - 保存结果的文件名
    
    tic;
    
    % 步骤1：计算粒子群的平均位置和速度，并获取单个粒子数据
    fprintf('步骤1: 计算粒子群 %d 的平均速度和位置...\n', groupNumber);
    [avgPosition, avgVelocity, individualVelocities, time] = calculateParticleGroupAveragesForTimescale(particleData, groupNumber);
    
    % 步骤2：提取流向分量（x方向速度）并进行处理
    fprintf('步骤2: 处理流向速度分量...\n');
    x = avgVelocity(:, 2);  % 第二列是平均Vx（流向分量）
    t = avgVelocity(:, 1);  % 第一列是时间
    
    % 去均值
    mean_x = mean(x);
    fprintf('平均流向速度: %.6f\n', mean_x);
    x = x - mean_x;
    
    % 步骤3：计算自相关函数
    fprintf('步骤3: 计算自相关函数...\n');
    [rho, lags] = my_xcorr(x, t, 0000);
    
    % 步骤4：积分自相关函数
    fprintf('步骤4: 积分自相关函数...\n');
    [x_inte, inte, zero_cross_timescale, zero_cross_point] = integrate_from_start(lags, rho);
    
    num_sample = 1;
    
    % 步骤5：拟合指数模型
    fprintf('步骤5: 拟合指数模型...\n');
    [~, t_l, r_squared, x_fit, y_fit] = lesutils.fit_exp_model(lags, rho, 0);
    t_l = -1 / t_l;  % 拟合系数为负倒数
    
    % 获取平均位置（使用最终时刻的位置）
    loc = avgPosition(end, 2:4);
    
    % 保存结果
    all_loc(num_sample, 1:3) = loc;
    all_t_l(num_sample) = t_l;
    all_r_squared(num_sample) = r_squared;
    
    % 步骤6：绘图
    fprintf('步骤6: 生成结果图表...\n');
    if(1)
        figure('Position', [100, 100, 1200, 500]);
        sgtitle(['Particle Group ', num2str(groupNumber), ' Lagrangian Timescale at Final Position (', ...
                num2str(loc(1), '%.4f'), ',', num2str(loc(2), '%.4f'), ',', ...
                num2str(loc(3), '%.4f'), ')']);
        
        % 第一个子图：自相关函数和拟合
        subplot(1, 3, 1);
        plot(lags, rho, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Autocorrelation Function');
        hold on;
        plot(x_fit, y_fit, 'r-.', 'LineWidth', 1.5, ...
             'DisplayName', ['Exponential Fit $\rho(s)=\exp(s/T_L)$' newline ' with $T_L = $', num2str(t_l, '%.4f')]);
        grid on;
        xlabel('Lag');
        ylabel('Autocorrelation');
        legend('Interpreter', 'latex', 'Location', 'best');
        title('Autocorrelation Function and Fit');
        
        % 第二个子图：自相关函数积分
        subplot(1, 3, 2);
        plot(x_inte, inte, 'b-', 'LineWidth', 1.5);
        hold on;
        line([0.5, x_inte(end)], [t_l, t_l], 'LineStyle', '--', 'Color', 'red', 'LineWidth', 1.5);
        line([0.5, x_inte(end)], [zero_cross_timescale, zero_cross_timescale], ...
             'LineStyle', '--', 'Color', 'green', 'LineWidth', 1.5);
        text(1, t_l + 0.005, ['$T_L = $', num2str(t_l, '%.4f')], ...
             'Interpreter', 'latex', 'Color', 'red', 'FontSize', 10);
        text(1, zero_cross_timescale + 0.005, ['$T_L = $', num2str(zero_cross_timescale, '%.4f')], ...
             'Interpreter', 'latex', 'Color', 'green', 'FontSize', 10);
        grid on;
        xlabel('Integration Limit');
        ylabel('Integral Value');
        title('Integration of Autocorrelation Function');
        legend('Integration', 'Fit Timescale', 'Zero-cross Timescale', 'Location', 'best');
        
        % 第三个子图：平均速度时间序列和单个粒子速度
        subplot(1, 3, 3);
        
        % 首先绘制所有单个粒子的速度（低透明度细线）
        numParticles = size(individualVelocities, 2);
        colors = lines(7); % 使用不同的颜色区分粒子
        
        for i = 1:numParticles
            % 计算单个粒子的平均速度
            particle_mean_vx = mean(individualVelocities(:, i));
            % 绘制去均值后的单个粒子速度
            plot(time, individualVelocities(:, i) - particle_mean_vx, ...
                 'Color', [colors(mod(i-1,7)+1, :), 0.2], ... % 低透明度
                 'LineWidth', 0.5, ...
                 'HandleVisibility', 'off'); % 不在图例中显示单个粒子
        end
        
        % 然后绘制平均速度（粗黑线）
        plot(t, x, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Group Mean Velocity (fluctuation)');
        hold on;
        
        % 绘制零线
        plot([t(1), t(end)], [0, 0], 'r--', 'LineWidth', 1.5, ...
             'DisplayName', 'Zero Level');
        
        grid on;
        xlabel('Time');
        ylabel('Streamwise Velocity Fluctuation (Vx'')');
        title('Streamwise Velocity Time Series');
        legend('show', 'Location', 'best');
        
        % 添加文本说明单个粒子数量
        text(0.02, 0.98, sprintf('Individual particles: %d', numParticles), ...
             'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'white', 'FontSize', 10);
    end
    
    % 步骤7：保存结果
    fprintf('步骤7: 保存结果...\n');
    save(['timescale_group', num2str(groupNumber), '_', filename], ...
         "all_r_squared", "all_t_l", "all_loc", "avgPosition", "avgVelocity", "individualVelocities");
    
    % 输出总结信息
    fprintf('\n=== 分析完成 ===\n');
    fprintf('粒子群编号: %d\n', groupNumber);
    fprintf('粒子数量: %d\n', numParticles);
    fprintf('最终位置: (%.4f, %.4f, %.4f)\n', loc);
    fprintf('Lagrangian时间尺度 (拟合): %.6f\n', t_l);
    fprintf('Lagrangian时间尺度 (零交叉): %.6f\n', zero_cross_timescale);
    fprintf('拟合优度 R²: %.6f\n', r_squared);
    fprintf('平均流向速度: %.6f\n', mean_x);
    
    toc;
end

function [avgPosition, avgVelocity, individualVelocities, time] = calculateParticleGroupAveragesForTimescale(particleData, groupNumber)
    % 专门为时间尺度分析计算粒子群平均位置和速度的辅助函数
    % 现在返回单个粒子的速度数据
    
    % 查找该粒子群中的所有粒子
    groupParticles = [];
    particleCount = 0;
    
    for i = 1:length(particleData)
        % 将id补0到5位数并提取粒子群编号
        id_str = sprintf('%05d', particleData{i}.id);
        particleGroup = str2double(id_str(1:2));
        
        if particleGroup == groupNumber
            particleCount = particleCount + 1;
            groupParticles(particleCount).data = particleData{i}.data;
        end
    end
    
    % 检查是否找到粒子
    if particleCount == 0
        error('未找到粒子群 %d 的任何粒子', groupNumber);
    end
    
    fprintf('  找到 %d 个粒子\n', particleCount);
    
    % 获取时间序列（假设所有粒子的时间序列相同）
    time = groupParticles(1).data(:, 1);
    numTimeSteps = length(time);
    
    % 初始化累加矩阵
    positionSum = zeros(numTimeSteps, 3);  % [x, y, z]
    velocitySum = zeros(numTimeSteps, 3);  % [vx, vy, vz]
    
    % 初始化单个粒子速度存储矩阵
    individualVelocities = zeros(numTimeSteps, particleCount);
    
    % 累加所有粒子的位置和速度
    validParticles = 0;
    for i = 1:particleCount
        particleDataMatrix = groupParticles(i).data;
        
        % 确保时间步数一致
        if size(particleDataMatrix, 1) == numTimeSteps
            positionSum = positionSum + particleDataMatrix(:, 2:4);
            velocitySum = velocitySum + particleDataMatrix(:, 5:7);
            
            % 存储单个粒子的流向速度
            individualVelocities(:, i) = particleDataMatrix(:, 5); % 第5列是vx
            
            validParticles = validParticles + 1;
        else
            fprintf('  警告: 粒子 %d 的时间步数不一致，跳过\n', i);
        end
    end
    
    % 只保留有效粒子的数据
    individualVelocities = individualVelocities(:, 1:validParticles);
    
    % 计算平均值
    avgPosition = [time, positionSum / validParticles];
    avgVelocity = [time, velocitySum / validParticles];
    
    fprintf('  有效粒子数: %d\n', validParticles);
    fprintf('  时间步数: %d\n', numTimeSteps);
    fprintf('  时间范围: %.2f ~ %.2f\n', time(1), time(end));
end
