function [x_out, integral_values, peak_value, peak_x] = integrate_from_start(x, y)
    % 检查输入数组长度是否一致
    if length(x) ~= length(y)
        error('输入数组x和y的长度必须相同。');
    end
    
    % 确保x是单调递增的
    if any(diff(x) < 0)
        error('自变量数组x必须是单调递增的。');
    end
    
    % 使用cumtrapz进行累积梯形积分
    integral_values = cumtrapz(x, y);
    
    % 输出自变量数组
    x_out = x;
    
    % 寻找第一个极大值点
    [peak_x, peak_value] = find_first_peak(x_out, integral_values);
end

function [peak_x, peak_value] = find_first_peak(x, y)
    % 寻找函数的第一个极大值点
    % 极大值点定义：该点的值大于其左右相邻点的值
    
    n = length(y);
    peak_x = [];
    peak_value = [];
    
    % 从第二个点开始检查到倒数第二个点
    for i = 2:n-1
        if y(i) > y(i-1) && y(i) > y(i+1)
            peak_x = x(i);
            peak_value = y(i);
            break; % 找到第一个极大值点后立即退出
        end
    end
    
    % 如果没有找到极大值点，检查最后一个点（可能是单调递增的情况）
    if isempty(peak_x) && n > 1 && y(n) > y(n-1)
        peak_x = x(n);
        peak_value = y(n);
    end
    
    % 如果仍然没有找到，检查第一个点（可能是单调递减的情况）
    if isempty(peak_x) && n > 1 && y(1) > y(2)
        peak_x = x(1);
        peak_value = y(1);
    end
    
    % 如果还是没有找到极大值点，返回空值
    if isempty(peak_x)
        warning('未找到极大值点。函数可能是单调的或常数函数。');
        peak_x = NaN;
        peak_value = NaN;
    end
end