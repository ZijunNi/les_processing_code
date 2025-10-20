function [R_avg, tau, T_L, T_L_u, T_L_v, T_L_w, R_u, R_v, R_w] = calculate_velocity_timescale(t, velocity, varargin)
%CALCULATE_VELOCITY_TIMESCALE 计算速度自相关时间尺度。
%   [R_AVG, TAU, T_L, T_L_U, T_L_V, T_L_W] =
%   CALCULATE_VELOCITY_TIMESCALE(T, VEL) 会基于速度序列 VEL 的三个分量
%   （列依次对应 u、v、w）计算空间平均自相关函数与积分时间尺度。
%
%   [...] = CALCULATE_VELOCITY_TIMESCALE(T, VEL, 'MaxLag', LAG) 指定计算自相关时考虑的最大延迟（采样点数）。

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 't', @(v) isnumeric(v) && isvector(v));
    addRequired(p, 'velocity', @(v) isnumeric(v) && size(v, 2) >= 1);
    addParameter(p, 'MaxLag', [], @(v) isempty(v) || (isscalar(v) && v >= 0));
    parse(p, t, velocity, varargin{:});

    t = p.Results.t(:);
    velocity = p.Results.velocity;
    maxLag = p.Results.MaxLag;

    if size(velocity, 1) ~= numel(t)
        error('calculate_velocity_timescale:SizeMismatch', ...
              'Length of time vector must match number of rows in velocity array.');
    end

    if isempty(maxLag)
        maxLag = min(length(t) - 1, floor(length(t) / 2));
    else
        maxLag = min(maxLag, length(t) - 1);
    end

    if maxLag < 1
        error('calculate_velocity_timescale:MaxLagTooSmall', ...
              'MaxLag must be at least 1 for meaningful correlation.');
    end

    % 保证拥有三个速度分量（不足时用零补足）
    if size(velocity, 2) < 3
        velocity = [velocity, zeros(size(velocity, 1), 3 - size(velocity, 2))];
    end

    [R_u, tau] = local_autocorr(t, velocity(:, 1), maxLag);
    [R_v, ~  ] = local_autocorr(t, velocity(:, 2), maxLag);
    [R_w, ~  ] = local_autocorr(t, velocity(:, 3), maxLag);

    R_avg = (R_u + R_v + R_w) / 3;

    [T_L_u, zero_cross_idx_u] = integrate_until_zero(tau, R_u);
    [T_L_v, zero_cross_idx_v] = integrate_until_zero(tau, R_v);
    [T_L_w, zero_cross_idx_w] = integrate_until_zero(tau, R_w);
    [T_L  , zero_cross_idx  ] = integrate_until_zero(tau, R_avg);

    % 若自相关未过零点，则给出简单的诊断提示
    if zero_cross_idx_u == numel(tau)
        warning('calculate_velocity_timescale:NoZeroCrossU', ...
                'u autocorrelation did not cross zero; integrated over entire range.');
    end
    if zero_cross_idx_v == numel(tau)
        warning('calculate_velocity_timescale:NoZeroCrossV', ...
                'v autocorrelation did not cross zero; integrated over entire range.');
    end
    if zero_cross_idx_w == numel(tau)
        warning('calculate_velocity_timescale:NoZeroCrossW', ...
                'w autocorrelation did not cross zero; integrated over entire range.');
    end
    if zero_cross_idx == numel(tau)
        warning('calculate_velocity_timescale:NoZeroCrossAvg', ...
                'Average autocorrelation did not cross zero; integrated over entire range.');
    end

    if nargout < 7
        clear R_u R_v R_w;
    end
end

function [rho, lags] = local_autocorr(t, component, maxLag)
    component = component(:);
    component = component - mean(component);
    [rho, lags] = my_xcorr(component, t, maxLag);
    rho = rho(:);
    lags = lags(:);
end

function [T_int, zero_idx] = integrate_until_zero(tau, rho)
    zero_idx = find(rho < 0, 1);
    if isempty(zero_idx)
        zero_idx = numel(rho);
    end
    T_int = trapz(tau(1:zero_idx), rho(1:zero_idx));
end
