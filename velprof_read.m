% 本代码读取particle_solver输出的"velprof-"文件并对其进行相关处理
% 指定目录
data = read_velprof_files('./data/coarse_drid_test');


%% 绘制平均速度剖面
if ~isempty(data)
    n = 40;
    % 访问第一个文件的平均速度数据
    first_file_mean_vel = data(n).mean_velocity;
    
    % 获取所有文件的数字标识
    file_numbers = [data.file_number];
    
    % 绘制某个文件的平均速度剖面
    figure;
    plot(data(n).coordinates, data(n).mean_velocity, 'b-', 'LineWidth', 2);
    hold on
    plot(data(n).coordinates, flip(data(n).mean_velocity), 'r-', 'LineWidth', 2);
    hold off
    xlabel('平均速度');
    ylabel('坐标');
    xlim([0 180])
    set(gca,'XScale','log')
    title(sprintf('文件 %d 的平均速度剖面', data(n).file_number));
end

%% 绘制RMS速度剖面
if ~isempty(data)
    n = 41;
    % 访问第一个文件的平均速度数据
    first_file_mean_vel = data(n).rms_velocity;
    
    % 获取所有文件的数字标识
    file_numbers = [data.file_number];
    
    % 绘制某个文件的平均速度剖面
    figure;
    plot(data(n).coordinates, data(n).rms_velocity, 'b-', 'LineWidth', 2);
    hold on
    plot(data(n).coordinates, flip(data(n).rms_velocity), 'r-', 'LineWidth', 2);
    hold off
    xlabel('平均速度');
    ylabel('坐标');
    xlim([0 180])
    % set(gca,'XScale','log')
    title(sprintf('文件 %d 的RMS速度剖面', data(n).file_number));
end




function velprof_data = read_velprof_files(directory_path)
% 读取特定目录下以velprof-开头的文件并整理为结构体数组
% 输入:
%   directory_path - 文件所在目录路径
% 输出:
%   velprof_data - 结构体数组，包含文件名数字和五列数据

    % 如果未提供目录路径，使用当前目录
    if nargin < 1
        directory_path = '.';
    end
    
    % 获取所有以velprof-开头的文件
    file_pattern = fullfile(directory_path, 'velprof-*');
    file_list = dir(file_pattern);
    
    if isempty(file_list)
        error('在目录 %s 中未找到以"velprof-"开头的文件', directory_path);
    end
    
    % 提取文件名中的数字并排序
    file_info = [];
    for i = 1:length(file_list)
        filename = file_list(i).name;
        
        % 提取velprof-后面的数字部分
        num_str = regexp(filename, 'velprof-(\d+)', 'tokens');
        if ~isempty(num_str)
            file_num = str2double(num_str{1}{1});
            file_info = [file_info; struct('filename', filename, 'number', file_num, 'index', i)];
        end
    end
    
    % 按数字大小排序
    [~, sort_idx] = sort([file_info.number]);
    file_info = file_info(sort_idx);
    
    % 初始化结构体数组
    velprof_data = struct();
    
    % 读取每个文件并存储数据
    for i = 1:length(file_info)
        file_path = fullfile(directory_path, file_info(i).filename);
        
        try
            % 读取五列数据
            data = load(file_path);
            
            % 验证数据格式
            if size(data, 2) ~= 5
                warning('文件 %s 不是五列数据，跳过该文件', file_info(i).filename);
                continue;
            end
            
            % 存储到结构体
            velprof_data(i).file_number = file_info(i).number;
            velprof_data(i).filename = file_info(i).filename;
            velprof_data(i).coordinates = data(:, 1);        % 第一列：坐标
            velprof_data(i).mean_velocity = data(:, 2);      % 第二列：平均速度
            velprof_data(i).rms_velocity = data(:, 3);       % 第三列：均方根速度
            velprof_data(i).vertical_mean_velocity = data(:, 4);  % 第四列：垂向平均速度
            velprof_data(i).spanwise_mean_velocity = data(:, 5);  % 第五列：展向平均速度
            
            fprintf('成功读取文件: %s (数字: %d)\n', file_info(i).filename, file_info(i).number);
            
        catch ME
            warning('读取文件 %s 时出错: %s', file_info(i).filename, ME.message);
        end
    end
    
    % 如果所有文件读取失败，返回空结构体
    if isempty(fieldnames(velprof_data))
        velprof_data = [];
        warning('未能成功读取任何文件');
    end
    
    fprintf('\n总共成功处理 %d 个文件\n', length(velprof_data));
end