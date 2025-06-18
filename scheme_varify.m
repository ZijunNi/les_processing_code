% 本代码主要对数值格式中的随机参数正定性进行验证
clear,clc;
% 设置参数
% clear,clc;
format long
Delta_t = 5e-3;     % 固定时间步长（需小于 T_i 和 tau_p）
T_max = 10000*Delta_t;       % T_i 扫描最大值
tau_max = 10000*Delta_t;     % tau_p 扫描最大值
num_points = 5000; % 网格分辨率

% 参数说明:
% Delta_t : 固定时间步长 (需满足 Delta_t < T_i 和 Delta_t < tau_p)
% T_max   : T_i 的扫描最大值
% tau_max : tau_p 的扫描最大值
% num_points : 每个参数的采样点数

% 生成 T_i 和 tau_p 的网格 (确保大于 Delta_t)
T_values = linspace(log(1*Delta_t) , log(T_max), num_points);    % 避免 T_i <= Delta_t
tau_values = linspace(log(1*Delta_t) , log(tau_max), num_points);% 避免 tau_p <= Delta_t
[T_grid, tau_grid] = meshgrid(exp(T_values), exp(tau_values));

% 预分配结果矩阵
result_matrix = zeros(size(T_grid));

%% 计算P_22^2：

% 逐点计算表达式值
for i = 1:num_points
    for j = 1:num_points
        f = omega_2(T_grid(i,j), tau_grid(i,j), Delta_t);
        h = omega_gamma(T_grid(i,j), tau_grid(i,j), Delta_t);
        g = gamma_2(T_grid(i,j), tau_grid(i,j), Delta_t);
        result_matrix(i,j) = (f-h/g);
    end
end



fprintf('最小值为%.6e\n',min(result_matrix(:)));
fprintf('最大值为%.6e\n',max(result_matrix(:)));

% 找到最小值点的位置
[min_val, min_idx] = min(result_matrix(:));
[min_row, min_col] = ind2sub(size(result_matrix), min_idx);
T_min = T_grid(min_row, min_col);
tau_min = tau_grid(min_row, min_col);
fprintf('最小值位于 T_i = %.6e (%.2fΔT) 和 tau_p = %.6e (%.2fΔT) 的交点\n',...
        T_min, T_min/Delta_t, tau_min, tau_min/Delta_t);


% 提取负值点的坐标
[neg_rows, neg_cols] = find(result_matrix < 0);
T_neg = T_grid(sub2ind(size(T_grid), neg_rows, neg_cols));
tau_neg = tau_grid(sub2ind(size(tau_grid), neg_rows, neg_cols));



% 绘制P_22^2等高线图
figure;
    contourf(T_grid/Delta_t, tau_grid/Delta_t, result_matrix, 'LineColor', 'none');
    
    % --- 计算数据范围并设置颜色映射范围 ---
    max_abs = max(abs(result_matrix(:))); % 获取数据绝对值的最大值
    clim([-max_abs, max_abs]);           % 设置对称的颜色范围
    
    % --- 创建自定义红绿渐变色阶 ---
    num_colors = 256;                    % 色阶总数
    half = floor(num_colors/2);          % 中间分割点
    
    % 红色渐变：深红[0.5,0,0] -> 浅红[1,0.7,0.7]
    red_ramp = [
        linspace(0.5, 1, half)' , ...  % R通道
        linspace(0, 0.7, half)' , ...  % G通道
        linspace(0, 0.7, half)'   ...  % B通道
    ];
    
    % 绿色渐变：浅绿[0.7,1,0.7] -> 深绿[0,0.5,0]
    green_ramp = [
        linspace(0.7, 0, half)'  , ... % R通道
        linspace(1, 0.5, half)'  , ... % G通道
        linspace(0.7, 0, half)'    ... % B通道
    ];
    
    % 合并色阶
    custom_cmap = [red_ramp; green_ramp];
    colormap(custom_cmap);
    
    % --- 坐标轴设置 ---
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('$T_i/\Delta t$','Interpreter','latex');
    ylabel('$\tau_p/\Delta t$','Interpreter','latex');
    title('P_{22}计算格式正定性验证');
    daspect([1 1 1])
    % --- 显示颜色条 ---
    c = colorbar;
    c.Label.String = '数值大小';          % 添加颜色条标签
    c.Ticks = linspace(-max_abs, max_abs, 5); % 设置5个刻度点
    grid on;


%% 计算P_33^2：

% 本代码主要对数值格式中的随机参数正定性进行验证
clear,clc;
% 设置参数
% clear,clc;
format long
Delta_t = 5e-3;     % 固定时间步长（需小于 T_i 和 tau_p）
T_max = 10000*Delta_t;       % T_i 扫描最大值
tau_max = 10000*Delta_t;     % tau_p 扫描最大值
num_points = 5000; % 网格分辨率

% 参数说明:
% Delta_t : 固定时间步长 (需满足 Delta_t < T_i 和 Delta_t < tau_p)
% T_max   : T_i 的扫描最大值
% tau_max : tau_p 的扫描最大值
% num_points : 每个参数的采样点数

% 生成 T_i 和 tau_p 的网格 (确保大于 Delta_t)
T_values = linspace(log(1*Delta_t) , log(T_max), num_points);    % 避免 T_i <= Delta_t
tau_values = linspace(log(1*Delta_t) , log(tau_max), num_points);% 避免 tau_p <= Delta_t
[T_grid, tau_grid] = meshgrid(exp(T_values), exp(tau_values));

% 预分配结果矩阵
result_matrix = zeros(size(T_grid));


% 逐点计算表达式值
for i = 1:num_points
    for j = 1:num_points
        P_21 = omega_gamma(T_grid(i,j), tau_grid(i,j), Delta_t)/sqrt(gamma_2(T_grid(i,j), tau_grid(i,j), Delta_t));
        P_31 = gamma_Gamma(T_grid(i,j), tau_grid(i,j), Delta_t)/sqrt(gamma_2(T_grid(i,j), tau_grid(i,j), Delta_t));
        squared_P_22 = omega_2(T_grid(i,j), tau_grid(i,j), Delta_t)-omega_gamma(T_grid(i,j), tau_grid(i,j), Delta_t)^2/gamma_2(T_grid(i,j), tau_grid(i,j), Delta_t);
        squared_P_32 = (Gamma_Omega(T_grid(i,j), tau_grid(i,j), Delta_t)-P_21*P_31)^2/squared_P_22;
        squared_P_33(i,j) = Gamma_2(T_grid(i,j), tau_grid(i,j), Delta_t)-P_31^2-squared_P_32;
        result_matrix(i,j) =  squared_P_33(i,j);
    end
end





fprintf('最小值为%.6e\n',min(result_matrix(:)));
fprintf('最大值为%.6e\n',max(result_matrix(:)));

% 找到最小值点的位置
[min_val, min_idx] = min(result_matrix(:));
[min_row, min_col] = ind2sub(size(result_matrix), min_idx);
T_min = T_grid(min_row, min_col);
tau_min = tau_grid(min_row, min_col);
fprintf('最小值位于 T_i = %.6e (%.2fΔT) 和 tau_p = %.6e (%.2fΔT) 的交点\n',...
        T_min, T_min/Delta_t, tau_min, tau_min/Delta_t);


% 提取负值点的坐标
[neg_rows, neg_cols] = find(result_matrix < 0);
T_neg = T_grid(sub2ind(size(T_grid), neg_rows, neg_cols));
tau_neg = tau_grid(sub2ind(size(tau_grid), neg_rows, neg_cols));



% 绘制P_22^2等高线图
figure;
    contourf(T_grid/Delta_t, tau_grid/Delta_t, result_matrix, 'LineColor', 'none');
    
    % --- 计算数据范围并设置颜色映射范围 ---
    max_abs = max(abs(result_matrix(:))); % 获取数据绝对值的最大值
    if(isinf(max_abs))
        max_abs = max(result_matrix(:));
    end
    clim([-max_abs, max_abs]);           % 设置对称的颜色范围
    
    % --- 创建自定义红绿渐变色阶 ---
    num_colors = 256;                    % 色阶总数
    half = floor(num_colors/2);          % 中间分割点
    
    % 红色渐变：深红[0.5,0,0] -> 浅红[1,0.7,0.7]
    red_ramp = [
        linspace(0.5, 1, half)' , ...  % R通道
        linspace(0, 0.7, half)' , ...  % G通道
        linspace(0, 0.7, half)'   ...  % B通道
    ];
    
    % 绿色渐变：浅绿[0.7,1,0.7] -> 深绿[0,0.5,0]
    green_ramp = [
        linspace(0.7, 0, half)'  , ... % R通道
        linspace(1, 0.5, half)'  , ... % G通道
        linspace(0.7, 0, half)'    ... % B通道
    ];
    
    % 合并色阶
    custom_cmap = [red_ramp; green_ramp];
    colormap(custom_cmap);
    
    % --- 坐标轴设置 ---
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('$T_i/\Delta t$','Interpreter','latex');
    ylabel('$\tau_p/\Delta t$','Interpreter','latex');
    title('P_{33}计算格式正定性验证');
    daspect([1 1 1])
    % --- 显示颜色条 ---
    c = colorbar;
    c.Label.String = '数值大小';          % 添加颜色条标签
    c.Ticks = linspace(-max_abs, max_abs, 5); % 设置5个刻度点
    grid on;


% 调用绘图函数
% plot_expression_positivity(Delta_t, T_max, tau_max, num_points);
% 
% 
% figure_name = '双正.pdf';
% filename = fullfile(pwd,'../figures',figure_name);
% exportgraphics(gcf, filename, 'ContentType', 'vector');


function result_matrix = plot_expression_positivity(Delta_t, T_max, tau_max, num_points)
    % 参数说明:
    % Delta_t : 固定时间步长 (需满足 Delta_t < T_i 和 Delta_t < tau_p)
    % T_max   : T_i 的扫描最大值
    % tau_max : tau_p 的扫描最大值
    % num_points : 每个参数的采样点数
    
    % 生成 T_i 和 tau_p 的网格 (确保大于 Delta_t)
    T_values = linspace(Delta_t , T_max, num_points);    % 避免 T_i <= Delta_t
    tau_values = linspace(Delta_t , tau_max, num_points);% 避免 tau_p <= Delta_t
    [T_grid, tau_grid] = meshgrid(T_values, tau_values);
    
    % 预分配结果矩阵
    result_matrix = zeros(size(T_grid));
    
    % 逐点计算表达式值
    for i = 1:num_points
        for j = 1:num_points
            result_matrix(i,j) = omega_2(T_grid(i,j), tau_grid(i,j), Delta_t);
        end
    end
    [~,I] = min(result_matrix(:));
    fprintf('最小值为%.6e\n',min(result_matrix(:)));

    fprintf('最大值为%.6e\n',max(result_matrix(:)));
        % 提取负值点的坐标
    [neg_rows, neg_cols] = find(result_matrix < 0);
    T_neg = T_grid(sub2ind(size(T_grid), neg_rows, neg_cols));
    tau_neg = tau_grid(sub2ind(size(tau_grid), neg_rows, neg_cols));
    
    % 绘制等高线图
    figure;
    contourf(T_grid/Delta_t, tau_grid/Delta_t, result_matrix, 'LineColor', 'none');
    hold on;
    
    % --- 新增：用红色星号标记负值点 ---
    % if ~isempty(T_neg)
    %     scatter(T_neg, tau_neg, 50, 'r*', 'LineWidth', 1.5, ...
    %            'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
    % end
    
    % 标注零等高线（正负分界线）
    contour(T_grid, tau_grid, result_matrix, [0, 0], ...
                    'LineWidth', 2, 'LineColor', 'k');
    
    % --- 修改颜色映射：负区域深红色，正区域深绿色 ---
    colormap([0.9 0.4 0.4; 0.4 0.9 0.4]); % [深红; 深绿]
    clim([-1e-5, 1e-5]);    
    % 标签和标题
    xlabel('$T_i/\Delta t$',Interpreter='latex');
    ylabel('$\tau_p/\Delta t$',Interpreter='latex');
    title(sprintf('正定性等高线图'));
    colorbar('Ticks', [], 'TickLabels', {}); % 隐藏色标数值（仅显示颜色意义）
    % text(T_max*0.8, tau_max*0.9, '负区域', 'Color', [0.7 0 0], 'FontSize', 12);
    % text(T_max*0.2, tau_max*0.2, '正区域', 'Color', [0 0.5 0], 'FontSize', 12);
    
    % --- 新增图例说明 ---
    % if ~isempty(T_neg)
    %     legend('负值点', 'Location', 'northwest');
    % end
    
    grid on;
end



% 统计量计算函数

function result = omega_2(T_i, tau_p, Delta_t)
    term1 = (T_i - tau_p)^2 * Delta_t;
    term2 = (T_i^3 / 2) * (1 - exp(-2 * Delta_t / T_i));
    term3 = (tau_p^3 / 2) * (1 - exp(-2 * Delta_t / tau_p));
    term4 = -2 * T_i^2 * (T_i - tau_p) * (1 - exp(-Delta_t / T_i));%原始为-
    term5 = 2 * tau_p^2 * (T_i - tau_p) * (1 - exp(-Delta_t / tau_p));
    term6 = -2 * (T_i^2 * tau_p^2) / (T_i + tau_p) * (1 - exp(-Delta_t / T_i) * exp(-Delta_t / tau_p));%原始为-
    result = term1 + term2 + term3 + term4 + term5 + term6;

    result = result/(Delta_t^3);%归一化
end


function result = omega_gamma(T_i, tau_p, delta_t)
    % 计算复杂表达式的三个部分
    term1 = (T_i - tau_p) * (1 - exp(-delta_t / T_i));
    term2 = (T_i / 2) * (1 - exp(-2 * delta_t / T_i));
    term3 = (tau_p^2 / (T_i + tau_p)) * (1 - exp(-delta_t / T_i) * exp(-delta_t / tau_p));
    
    % 组合各项并乘以 T_i 得到最终结果
    sum_terms = term1 - term2 + term3;
    result = (T_i * sum_terms)/(delta_t^2);
end

function result = gamma_2(T_i, tau_p, delta_t)
    % 计算 <γ_i^2(t)> 对应公式 (136)
    term = 1 - exp(-2*delta_t / T_i);
    result = (T_i / 2) * term / (delta_t);
end


function cov_Gamma2 = Gamma_2(T_i, tau_p, dt)
    % 计算 <Γ_i^2(t)> 对应公式 (137)
    term1 = T_i/2 * (1 - exp(-2*dt/T_i));
    term2 = - (2*tau_p*T_i)/(T_i + tau_p) * (1 - exp(-dt/T_i)*exp(-dt/tau_p));
    term3 = tau_p/2 * (1 - exp(-2*dt/tau_p));
    cov_Gamma2 = (term1 + term2 + term3);
end

% function cov_Omega2 = Omega_2(dt, T_i, tau_p, theta_i, B_i)
%     % 计算 <Ω_i^2(t)> 对应公式 (138)
%     term1 = (T_i - tau_p)^2 * dt;
%     term2 = (T_i^3 / 2) * (1 - exp(-2*dt/T_i));
%     term3 = (tau_p^3 / 2) * (1 - exp(-2*dt/tau_p));
%     term4 = -2*T_i^2*(T_i - tau_p) * (1 - exp(-dt/T_i));
%     term5 = 2*tau_p^2*(T_i - tau_p) * (1 - exp(-dt/tau_p));
%     term6 = -2*(T_i^2 * tau_p^2)/(T_i + tau_p) * (1 - exp(-dt/T_i)*exp(-dt/tau_p));
%     normalized = term1 + term2 + term3 + term4 + term5 + term6;
%     cov_Omega2 = normalized;
% end

function cov_gamma_Gamma = gamma_Gamma(T_i, tau_p, dt)
    % 计算 <γ_i(t)Γ_i(t)> 对应公式 (139)
    term1 = 1/2 * (1 - exp(-2*dt/T_i));
    term2 = - tau_p/(T_i + tau_p) * (1 - exp(-dt/T_i)*exp(-dt/tau_p));
    cov_gamma_Gamma = T_i * (term1 + term2);
end

function cov_gamma_Omega = gamma_Omega(T_i, tau_p, dt)
    % 计算 <γ_i(t)Ω_i(t)> 对应公式 (140)
    term1 = (T_i - tau_p) * (1 - exp(-dt/T_i));
    term2 = - (T_i/2) * (1 - exp(-2*dt/T_i));
    term3 = (tau_p^2)/(T_i + tau_p) * (1 - exp(-dt/T_i)*exp(-dt/tau_p));
    cov_gamma_Omega = T_i * (term1 + term2 + term3);
end

function cov_Gamma_Omega = Gamma_Omega(T_i, tau_p, dt)
    % 计算 <Γ_i(t)Ω_i(t)> 对应公式 (141)
    term1 = (T_i - tau_p) * (T_i*(1 - exp(-dt/T_i)) - tau_p*(1 - exp(-dt/tau_p)));
    term2 = - (T_i^2/2) * (1 - exp(-2*dt/T_i));
    term3 = - (tau_p^2/2) * (1 - exp(-2*dt/tau_p));
    term4 = T_i * tau_p * (1 - exp(-dt/T_i)*exp(-dt/tau_p));
    normalized = term1 + term2 + term3 + term4;
    cov_Gamma_Omega =  normalized;
end