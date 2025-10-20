function plotParticleGroupTrajectories(particleData, groupNumbers)
    % 绘制指定粒子群的中心粒子轨迹及与其他粒子的连线
    %
    % 输入参数：
    %   particleData - 结构体数组，包含data和id字段
    %   groupNumbers - 粒子群编号（可以是单个数字或数组），范围0-63
    
    % 创建图形窗口
    figure;
    hold on;
    grid on;
    
    % 设置图形属性
    xlabel('X');
    ylabel('Y'); 
    zlabel('Z');
    title(['粒子群 ', num2str(groupNumbers), ' 的轨迹与连线']);
    
    % 设置固定的坐标轴范围
    xlim([0, 4*pi]);
    ylim([0, 2*pi]);
    zlim([0, 2]);
    
    % 确保groupNumbers是数组形式
    if isscalar(groupNumbers)
        groupNumbers = [groupNumbers];
    end
    
    % 颜色映射，为不同的粒子群分配不同颜色
    colors = lines(length(groupNumbers));
    
    % 遍历所有指定的粒子群
    for g = 1:length(groupNumbers)
        currentGroup = groupNumbers(g);
        
        % 查找该粒子群中的所有粒子
        groupParticles = [];
        for i = 1:length(particleData)
            % 将id补0到5位数并提取粒子群编号和粒子编号
            id_str = sprintf('%05d', particleData{i}.id);
            particleGroup = str2double(id_str(1:2));
            particleId = str2double(id_str(3:5));
            
            if particleGroup == currentGroup
                groupParticles(particleId).data = particleData{i}.data;
                groupParticles(particleId).id = particleId;
            end
        end
        
        % 检查是否找到中心粒子（编号为1）
        if length(groupParticles) < 1 || isempty(groupParticles(1).data)
            fprintf('警告: 未找到粒子群 %d 的中心粒子\n', currentGroup);
            continue;
        end
        
        % 获取中心粒子的轨迹数据
        centerData = groupParticles(1).data;
        centerPositions = centerData(:, 2:4);  % 中心粒子位置 (x, y, z)
        
        % 绘制中心粒子的轨迹
        plot3(centerPositions(:, 1), centerPositions(:, 2), centerPositions(:, 3), ...
              'Color', colors(g, :), 'LineWidth', 2.5, ...
              'DisplayName', ['粒子群 ', num2str(currentGroup), ' 中心']);
        
        % 标记中心粒子的起点和终点
        plot3(centerPositions(1, 1), centerPositions(1, 2), centerPositions(1, 3), ...
              'o', 'Color', colors(g, :), 'MarkerSize', 8, ...
              'MarkerFaceColor', colors(g, :), 'HandleVisibility', 'off');
        plot3(centerPositions(end, 1), centerPositions(end, 2), centerPositions(end, 3), ...
              's', 'Color', colors(g, :), 'MarkerSize', 8, ...
              'MarkerFaceColor', colors(g, :), 'HandleVisibility', 'off');
        
        % 绘制与其他六个粒子的连线
        % 假设其他粒子的编号为2-7（上下左右前后）
        neighborIds = 2:7;
        for n = neighborIds
            if length(groupParticles) >= n && ~isempty(groupParticles(n).data)
                neighborData = groupParticles(n).data;
                neighborPositions = neighborData(:, 2:4);
                
                % 绘制从中心粒子到相邻粒子的连线
                % 可以选择在特定时间点绘制连线，避免过于密集
                timeStep = max(1, floor(size(centerPositions, 1) / 20)); % 约20个连线
                for t = 1:timeStep:size(centerPositions, 1)
                    % 确保时间点有效
                    if t <= size(neighborPositions, 1)
                        plot3([centerPositions(t, 1), neighborPositions(t, 1)], ...
                              [centerPositions(t, 2), neighborPositions(t, 2)], ...
                              [centerPositions(t, 3), neighborPositions(t, 3)], ...
                              'Color', 'red', 'LineStyle', '-', ...
                              'LineWidth', 0.5, 'HandleVisibility', 'off');
                    end
                end
                
                % 标记相邻粒子的位置（可选）
                plot3(neighborPositions(1, 1), neighborPositions(1, 2), neighborPositions(1, 3), ...
                      '^', 'Color', colors(g, :), 'MarkerSize', 4, ...
                      'MarkerFaceColor', colors(g, :), 'HandleVisibility', 'off');
            end
        end
        
        fprintf('粒子群 %d: 已绘制中心粒子轨迹及与相邻粒子的连线\n', currentGroup);
    end
    
    % 添加图例和调整视图
    legend('show', 'Location', 'best');
    view(3);  % 3D视图
    axis equal;
    hold off;
end