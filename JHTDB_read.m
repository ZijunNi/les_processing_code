%% 本代码用于读取并合并JHTDB数据库中下载的单点时序数据
clear,clc
% 设置子目录路径
subdir = './data/JHTDB_data'; % 请替换为您的子目录名称

% 获取所有mat文件
files = dir(fullfile(subdir, '*_to_*.mat'));
nFiles = length(files);

% 提取文件名中的数字并排序
fileNumbers = zeros(nFiles, 1);
for i = 1:nFiles
    [~, name, ~] = fileparts(files(i).name);
    numbers = sscanf(name, '%d_to_%d');
    fileNumbers(i) = numbers(1);
end
[sortedNumbers, sortedIdx] = sort(fileNumbers);
files = files(sortedIdx);

% 初始化变量
all_time = [];
all_data = [];
prev_end_time = NaN;
prev_end_data = NaN;

% 循环处理每个文件
for i = 1:nFiles
    % 加载数据
    data = load(fullfile(subdir, files(i).name));
    
    % 获取变量（假设时间变量名为't'，数据变量名为'v'）
    t = data.times_plot;
    v = data.result;
    
    % 检查数据维度并压缩
    v = squeeze(v);
    
    % 检查连续性（从第二个文件开始）
    if i > 1
        time_diff = t(1) - prev_end_time;
        data_diff = max(abs(v(1, :) - prev_end_data), [], 'all');
        
        if time_diff > 1e-10 || data_diff > 1e-10
            warning('文件 %s 与上一个文件不连续', files(i).name);
            disp(data_diff);
        else
            % 去除重复的时间点和数据
            t(1) = [];
            v(1, :) = [];
        end
    end
    
    % 连接数据
    all_time = [all_time; t(:)];
    all_data = [all_data; v];
    
    % 保存当前文件的末尾值用于下一个文件的检查
    prev_end_time = t(end);
    prev_end_data = v(end, :);
    % pause
end

% 保存结果
save('combined_data.mat', 'all_time', 'all_data', '-v7.3');
disp('数据已保存到 combined_data.mat');

%% 计算单点Euler时间尺度
clear,clc
load("combined_data.mat")
[R_avg, tau, T_L, T_L_u, T_L_v, T_L_w, R_u, R_v, R_w] = lesutils.calculate_velocity_timescale(all_time, all_data);

zero_cross_idx_u = find(R_u < 0, 1);
if isempty(zero_cross_idx_u), zero_cross_idx_u = numel(R_u); end
zero_cross_idx_v = find(R_v < 0, 1);
if isempty(zero_cross_idx_v), zero_cross_idx_v = numel(R_v); end
zero_cross_idx_w = find(R_w < 0, 1);
if isempty(zero_cross_idx_w), zero_cross_idx_w = numel(R_w); end
zero_cross_idx = find(R_avg < 0, 1);
if isempty(zero_cross_idx), zero_cross_idx = numel(R_avg); end

figure('Color', 'white')
plot(tau, R_u, 'r--', 'LineWidth', 1.2)
hold on
plot(tau, R_v, 'g--', 'LineWidth', 1.2)
plot(tau, R_w, 'b--', 'LineWidth', 1.2)
plot(tau, R_avg, 'k-', 'LineWidth', 2.5)
plot(tau(1:zero_cross_idx), R_avg(1:zero_cross_idx), 'm', 'LineWidth', 2)
yline(0, 'k-')
xlabel('Time $t$', 'FontSize', 12, 'Interpreter','latex')
ylabel('Autocorrelation Function of Velocity $R(\tau)$', 'FontSize', 12, 'Interpreter','latex')
legend({'Streamwise', 'Wall-normal', 'Spanwise', 'Average','Integral Range'}, ...
       'Location', 'best', 'FontSize', 10, 'Interpreter','latex')
grid on
set(gca, 'FontSize', 11)
       
fprintf('\n拉格朗日时间尺度结果:\n');
fprintf('u方向: T_L_u = %.4f (积分上限 τ_u = %.4f)\n', T_L_u, tau(zero_cross_idx_u));
fprintf('v方向: T_L_v = %.4f (积分上限 τ_v = %.4f)\n', T_L_v, tau(zero_cross_idx_v));
fprintf('w方向: T_L_w = %.4f (积分上限 τ_w = %.4f)\n', T_L_w, tau(zero_cross_idx_w));
fprintf('平均: T_L = %.4f (积分上限 τ = %.4f)\n', T_L, tau(zero_cross_idx));
% xlim([0 1])
