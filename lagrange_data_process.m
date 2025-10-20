% clear,clc
% % 主函数：调用数据读取和物理验证
% folderPath = 'F:\本机文档\科研\local_test\particle_solver\runscript\lagrange_data'; 
folderPath = 'F:\本机文档\data\bspline_test'; 
% folderPath = 'F:\本机文档\科研\local_test\particle_solver\runscript\lagrange_data';

% 读取数据
particleData = load_lagrange_data(folderPath,'',1,1);
particleData(cellfun(@isempty,particleData))=[];
fprintf('成功读取 %d 个时间步数据\n', length(particleData));
% 生成格式化的时间字符串
timeStr = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM-SS');
fileName = ['processed_data\fluid_lagrange_', timeStr, '.mat'];
save(fileName,"particleData")


%% 物理验证
clc
[meanErr, maxErr] = verify_lagrange_data(particleData,1e-4);

%% MPI验证
[x,y,z] = plot_lagrange_velocity(particleData,true,'xy_u');%第一个分量为绘图的横坐标
hold on
axis("equal")
xlim([0*pi 4*pi])
% ylim([0*pi 2*pi])
ylim([0 2])
xx = 0:0.01:4*pi;
yy = 2*ones(size(xx));
line(xx,yy,'Color','black','linewidth',2)
hold on
yy = 0*ones(size(xx));
line(xx,yy,'Color','black','linewidth',2)
% 设置 x 轴刻度和标签
xticks(0:pi:4*pi);
xticklabels({ '0','\pi','2\pi','3\pi','4\pi'});
hold off
% yticks(0:pi/2:2*pi);
% yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});

%% 流体粒子速度剖面

load("grid_test.mat")
bin_edges = grid_test.data(:,1);
[bin_centers,mean_vel,rms_vel] = calculate_mean_velocity(particleData,bin_edges);
% 绘制折线图
figure;
plot(bin_centers*180, mean_vel/0.0643994, 'o-', 'MarkerSize', 6,'DisplayName','$\overline{u^+}$');
% 设置图形属性
xlabel('$y^+$',Interpreter='latex');
ylabel('Mean Streamwise Velocity of Fluid Particles $\overline{u^+}$',Interpreter='latex');
set(gca,'XScale','log')
xlim([0, 180]);
hold on
plot(bin_edges(1:25)*180,bin_edges(1:25)*180,'-.',LineWidth=2,DisplayName='Linear Law, $u^+ = y^+$')
plot(bin_edges*180,log(bin_edges*180)/0.4+5,'-.',LineWidth=2,DisplayName='Log Law, $u^+ = \log(y^+)/0.4+5$')
legend(Interpreter="latex")
daspect([1/(19) 1/((180)) 1])
%% 流体粒子的均方根剖面
figure;
plot(bin_centers*180,rms_vel/0.0643994,'o-', 'MarkerSize', 6,'DisplayName','$\sqrt{\overline{(u^\prime)^2}}/u_\tau$');
% legend(Interpreter="latex")
xlabel('$y^+$',Interpreter='latex');
ylabel('$\sqrt{\overline{(u^\prime)^2}}/u_\tau$',Interpreter='latex');
xlim([0, 180]);
daspect([1/(3) 1/((180)) 1])

%% 临时验证区
clc;
all_data = particleData;
for i = 1:length(all_data)
    arr = all_data{i}.data;
    id = all_data{i}.id;
    set_vel = fix(id/1000)/64;
    err_count = 0;
    for j = 1:size(arr,1)
        err = abs(arr(j,5)-set_vel);
        if(err>1/64)
            err_count = err_count + 1;
            % fprintf('编号为%d的粒子速度有误，在第%d行，其理论速度为%.4f，实际速度为%.4f。\n',id,j,set_vel,arr(j,5));
        end
    end
    if(err_count==size(arr,1))
           fprintf('编号为%d的粒子速度有误，在第%d行，其理论速度为%.4f，实际速度为%.4f。\n',id,j,set_vel,arr(j,5));
    end
end
%% 函数部分


function particleData = load_lagrange_data(folderPath,file_prefix,sample_begin,sample_step)
    % 读取特定文件夹下所有fluid_*.txt文件并排序存储
    % sample_step - 为采样样本间隔

    fileList = dir(fullfile(folderPath, [file_prefix,'*.txt']));
    if isempty(fileList)
        error('未找到符合命名规则fluid_*.txt的文件');
    end
    
    % 提取文件编号并排序
    fileNumbers = zeros(1, length(fileList));
    for i = 1:length(fileList)
        [~, name] = fileparts(fileList(i).name);
        numStr = regexp(name, ['(?<=',file_prefix,'_)\d+'], 'match', 'once');
        fileNumbers(i) = str2double(numStr);
    end
    [~, sortedIdx] = sort(fileNumbers);
    sortedFiles = fileList(sortedIdx);
    disp('排序完成。');
    % 读取文件数据到元胞数组
    particleData = cell(1, length(sortedFiles));
    for i = sample_begin:sample_step:length(sortedFiles)
        filePath = fullfile(folderPath, sortedFiles(i).name);
        part_data.data = load(filePath);
        part_data.id = fileNumbers(i);
        particleData{i} = part_data;
        disp(['已经读取',num2str(i),'/',num2str(length(sortedFiles)),'个。'])
    end
end

function [isValid, invalidIndices] = verify_lagrange_data(all_data, tolerance)
    %   [isValid, invalidIndices] = validatePhysics(data, tolerance)
    %
    % 输入:
    %   data      - n×7 矩阵，列顺序: [时间, x, y, z, vx, vy, vz]
    %   tolerance - 允许的误差容限（标量）
    %
    % 输出:
    %   isValid       - 是否所有数据点均通过验证 (true/false)
    %   invalidIndices - 未通过验证的数据行索引（从第二行开始）
    %
    % 边界条件:
    %   x方向: 周期边界 (0 ~ 4*pi)
    %   z方向: 周期边界 (0 ~ 2*pi)
    %   y方向: 反射边界 (0 和 2)
    
    % 检查数据点数量
for j = 1:length(all_data)
    data = all_data{j}.data;
    if size(data, 1) < 2
        error('至少需要两个数据点才能进行验证');
    end
    
    % 定义周期边界范围
    Lx = 4 * pi;  % x方向周期
    Lz = 2 * pi;  % z方向周期
    y_min = 0;    % y方向下边界
    y_max = 2;    % y方向上边界
    
    % 提取时间、位置和速度
    times = data(:, 1);
    positions = data(:, 2:4);
    velocities = data(:, 5:7);
    
    % 计算时间间隔 (Δt = t_i - t_{i-1})
    dt = diff(times);
    
    % 计算实际位移变化 (考虑边界条件)
    actualDisplacement = zeros(size(positions, 1)-1, 3);
    
    for i = 2:size(positions, 1)
        % 前一时间点位置
        prev_pos = positions(i-1, :);
        
        % 当前时间点位置
        curr_pos = positions(i, :);
        
        % 1. x方向：周期边界处理
        dx = curr_pos(1) - prev_pos(1);
        dx_corrected = dx - round(dx/Lx) * Lx;  % 最小位移差
        
        % 2. z方向：周期边界处理
        dz = curr_pos(3) - prev_pos(3);
        dz_corrected = dz - round(dz/Lz) * Lz;  % 最小位移差
        
        % 3. y方向：反射边界处理
        % 由于反射边界，实际位移可能经过多次反射
        % 这里使用直线位移（反射边界下的实际位移变化等于位置差）
        y_pos = prev_pos(2) + velocities(i,2)*dt(i-1);
        if(y_pos<y_min)
            dy = - curr_pos(2) - prev_pos(2);
        elseif(y_pos>y_max)
            dy = y_max*2 - curr_pos(2) - prev_pos(2);
        else
            dy = curr_pos(2) - prev_pos(2);
        end
        
        actualDisplacement(i-1, :) = [dx_corrected, dy, dz_corrected];
    end
    
    % 计算理论位移变化 (Δp_theoretical = v_{i-1} * Δt)
    theoreticalDisplacement = velocities(1:end-1, :) .* dt;
    
    % 计算误差
    error = actualDisplacement - theoreticalDisplacement;
    
    % 标记超标点（从第二行开始验证）
    invalidMask = any(abs(error) > tolerance, 2);
    
    % 输出结果
    invalidIndices = find(invalidMask) + 1;  % +1 对应原始数据行号
    isValid = isempty(invalidIndices);

    % 显示验证结果
    if isValid
        % disp('所有相邻时间点间的位移变化验证通过');
    else
        fprintf('编号为 %d 的粒子的%d个数据点未通过验证，相关信息如下：\n',all_data{j}.id,numel(invalidIndices));
        % for i =1:length(invalidIndices)
        %     fprintf('| 行号 |  时间  |                 位置                 |      速度     |\n')
        %     fprintf(' %d  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g\n',invalidIndices-1,data(invalidIndices-1,:))
        %     fprintf(' %d  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g\n',invalidIndices,data(invalidIndices,:))
        %     fprintf(' %d  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g\n',invalidIndices+1,data(invalidIndices+1,:))
        % end
        fprintf('最大误差: [X:%.4g, Y:%.4g, Z:%.4g]\n\n',max(abs(error)));
    end
end
end

function [x,y,z] = plot_lagrange_velocity(cellArray, plot_switch, component_str)
% PLOT_LAGRANGE_VELOCITY 绘制粒子轨迹并用颜色表示速度分量
%   [x,y,z] = plot_lagrange_velocity(cellArray, plot_switch, component_str)
%   输入: 
%     cellArray - 元胞数组，每个元素为N×7矩阵
%                矩阵列: [?, X, Y, Z, Vx, Vy, Vz]
%                其中第2-4列为位置(X,Y,Z)，第5-7列为速度分量(Vx,Vy,Vz)
%     plot_switch - 绘图开关 (1=绘图, 0=不绘图)
%     component_str - 分量选择字符串 (格式: '平面_分量', 如 'xy_u')
%                 平面: 'xy', 'xz', 'yz'
%                 分量: 'u'(Vx), 'v'(Vy), 'w'(Vz)
%
%   输出: 二维散点图，点颜色表示速度分量值
%         当输入只有一个粒子时，额外绘制粒子轨迹线

% 设置默认分量
if nargin < 3
    component_str = 'xy_u'; % 默认在xy平面绘制u分量
end

% 检查输入有效性
if ~iscell(cellArray)
    error('输入必须是元胞数组');
end

% 解析分量选择字符串
if length(component_str) ~= 4 || component_str(3) ~= '_'
    error('分量字符串格式错误，应为类似 "xy_u" 格式');
end

plane = component_str(1:2);
comp_char = component_str(4);

% 确定位置分量索引
switch plane
    case 'xy'
        x_col = 2; % X位置
        y_col = 3; % Y位置
        x_label = '$x$';
        y_label = '$y$';
    case 'xz'
        x_col = 2; % X位置
        y_col = 4; % Z位置
        x_label = '$x$';
        y_label = '$z$';
    case 'yz'
        x_col = 3; % Y位置
        y_col = 4; % Z位置
        x_label = '$y$';
        y_label = '$z$';
    case 'yx'
        x_col = 3; % Y位置
        y_col = 2; % Z位置
        x_label = '$y$';
        y_label = '$x$';
    case 'zx'
        x_col = 4; % Y位置
        y_col = 2; % Z位置
        x_label = '$z$';
        y_label = '$x$';
    case 'zy'
        x_col = 4; % Y位置
        y_col = 3; % Z位置
        x_label = '$z$';
        y_label = '$y$';
    otherwise
        error('不支持的平面选择，请选择 "xy", "xz", 或 "yz"');
end

% 确定速度分量索引
switch comp_char
    case 'u'
        v_col = 5; % Vx
        c_label = '流向速度u_x';
    case 'v'
        v_col = 6; % Vy
        c_label = '流向速度u_y';
    case 'w'
        v_col = 7; % Vz
        c_label = '流向速度u_z';
    otherwise
        error('不支持的速度分量，请选择 "u", "v", 或 "w"');
end

% 计算总数据点数
numPoints = 0;
for i = 1:numel(cellArray)
    if size(cellArray{i}.data, 2) ~= 7
        error('每个元胞元素必须是7列矩阵');
    end
    numPoints = numPoints + size(cellArray{i}.data, 1);
end

% 预分配内存
allXY = zeros(numPoints, 2);  % 存储绘图位置 [X, Y] 或 [X,Z] 等
allVel = zeros(numPoints, 1); % 存储速度分量值

% 提取并合并所有数据
idx = 1;
for i = 1:numel(cellArray)
    currentMatrix = cellArray{i}.data;
    n = size(currentMatrix, 1);
    
    % 提取位置和速度分量
    allXY(idx:idx+n-1, :) = [currentMatrix(:, x_col), currentMatrix(:, y_col)];
    allVel(idx:idx+n-1) = currentMatrix(:, v_col);
    
    idx = idx + n;
end

x = allXY(:,1);
y = allXY(:,2);
z = allVel;

% 绘图部分
if plot_switch
    figure;
    hold on;
    
    % 如果只有一个粒子，绘制轨迹线
    if numel(cellArray) == 1 && size(cellArray{1}.data, 1) > 1
        singleData = cellArray{1}.data;
        plot(singleData(:, x_col), singleData(:, y_col), ...
            'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
    
    % 绘制散点图
    sc = scatter(x, y, 12, z, 'filled');
    
    % 设置图形属性
    colormap(jet);
    colorbar;
    xlabel(x_label, 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(y_label, 'Interpreter', 'latex', 'FontSize', 12);
    grid on;
    
    % 设置颜色条
    c = colorbar;
    c.Label.String = c_label;
    c.Label.FontSize = 12;
    
    % 设置颜色轴范围
    velRange = [min(allVel), max(allVel)];
    if velRange(1) ~= velRange(2)
        caxis(velRange);
    end
    
    % 添加图例说明（如果是单粒子轨迹）
    if numel(cellArray) == 1
        annotation('textbox', [0.8, 0.05, 0.15, 0.05], ...
            'String', '轨迹线', ...
            'Color', 'k', ...
            'EdgeColor', 'none', ...
            'BackgroundColor', 'none', ...
            'FontSize', 10);
    end
    
    % hold off;
end
end

function [bin_centers,mean_velocities,rms_fluctuations] = calculate_mean_velocity(cell_data,bin_edges)
    % 定义分区边界 (0~2.1 以包含所有数据)
    num_bins = length(bin_edges) - 1;
    
    % 初始化累加器
    sum_velocity = zeros(num_bins, 1);
    sum_velocity_sq = zeros(num_bins, 1);  % 新增：速度平方的累加
    count_points = zeros(num_bins, 1);
    
    % 第一遍遍历：计算每个区间的平均速度
    for k = 1:length(cell_data)
        data_matrix = cell_data{k}.data;
        
        % 提取第三列(y坐标)和第五列(速度)
        y_values = data_matrix(:, 3);
        velocities = data_matrix(:, 5);
        
        % 分配数据点到区间
        bin_indices = discretize(y_values, bin_edges);
        
        % 遍历每个区间
        for bin_idx = 1:num_bins
            % 找出当前区间的数据点
            mask = (bin_indices == bin_idx);
            current_velocities = velocities(mask);
            
            % 累加速度和计数
            if ~isempty(current_velocities)
                sum_velocity(bin_idx) = sum_velocity(bin_idx) + sum(current_velocities);
                sum_velocity_sq(bin_idx) = sum_velocity_sq(bin_idx) + sum(current_velocities.^2); % 累加速度平方
                count_points(bin_idx) = count_points(bin_idx) + length(current_velocities);
            end
        end
    end
    
    % 计算每个区间的平均速度
    mean_velocities = sum_velocity ./ count_points;
    
    % 计算速度脉动的均方根 (RMS)
    % 公式：RMS = sqrt(mean(v^2) - mean(v)^2)
    mean_velocity_sq = sum_velocity_sq ./ count_points;
    rms_fluctuations = sqrt(mean_velocity_sq - mean_velocities.^2);
    
    % 计算区间中心点作为x坐标
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;    
end

% 合并数据集（注意需要编号一一对应，按particleData_1先particleData_2后的顺序排列
function mergedParticleData = mergeParticleData(particleData_1, particleData_2)
    % 去除两个变量中的空元胞
    particleData_1 = particleData_1(~cellfun('isempty', particleData_1));
    particleData_2 = particleData_2(~cellfun('isempty', particleData_2));
    
    % 提取所有 id
    ids1 = cellfun(@(x) x.id, particleData_1);
    ids2 = cellfun(@(x) x.id, particleData_2);
    
    % 检查 id 是否一一对应
    if numel(ids1) ~= numel(ids2) || ~all(sort(ids1) == sort(ids2))
        error('两个 particleData 的 id 不匹配，无法合并');
    end
    
    % 为 particleData_2 创建 id 到索引的映射
    idMap = containers.Map(ids2, 1:numel(ids2));
    
    % 预分配合并后的元胞数组
    mergedParticleData = cell(size(particleData_1));
    
    % 遍历 particleData_1 中的每个元素
    for i = 1:numel(particleData_1)
        currentId = particleData_1{i}.id;
        
        % 检查 particleData_2 中是否存在相同 id
        if ~isKey(idMap, currentId)
            error('在 particleData_2 中找不到 id=%d 的数据', currentId);
        end
        
        % 获取两个数据集
        data1 = particleData_1{i}.data;
        data2 = particleData_2{idMap(currentId)}.data;
        
        % 检查列数是否一致
        if size(data1, 2) ~= size(data2, 2)
            error('id=%d 的数据列数不一致: %d vs %d', currentId, size(data1,2), size(data2,2));
        end
        
        % 垂直合并数据
        mergedData = [data1; data2];
        
        % 创建合并后的结构体
        mergedStruct = struct('id', currentId, 'data', mergedData);
        mergedParticleData{i} = mergedStruct;
    end
end