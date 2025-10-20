function [A, B, r_squared, x_fit, y_fit] = fit_exp_model(x, y, fixedA)
%FIT_EXP_MODEL 拟合形式为 y = A + exp(B * x) 的指数模型。
%   [A, B, R2, X_FIT, Y_FIT] = FIT_EXP_MODEL(X, Y) 对数据 (X, Y) 执行非线性
%   最小二乘拟合，目标模型为 y = A + exp(B * x)。
%
%   ... = FIT_EXP_MODEL(X, Y, FIXEDA) 将偏移量 A 固定为 FIXEDA，仅拟合指数系数 B。
%
%   返回拟合得到的参数 A、B，自决定系数 R^2，以及在 X 取值范围上采样的拟合曲线。

    narginchk(2, 3);

    if nargin < 3 || isempty(fixedA)
        fixedA = 1;
    end

    if numel(x) ~= numel(y)
        error('fit_exp_model:InputSizeMismatch', ...
              'x and y must have the same number of elements.');
    end

    % 移除 NaN 并按 x 排序以提高拟合稳定性
    validIdx = ~(isnan(x) | isnan(y));
    x = x(validIdx);
    y = y(validIdx);

    if isempty(x)
        error('fit_exp_model:EmptyInput', ...
              'No valid data points remain after removing NaNs.');
    end

    [x, sortIdx] = sort(x);
    y = y(sortIdx);

    % 为可视化/输出准备一组稠密的 x 网格
    x_fit = linspace(min(x), max(x), 1000);

    options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8);

    if ~isnan(fixedA)
        % 固定 A，只拟合指数系数 B
        A = fixedA;
        model = @(B, xq) fixedA + exp(B .* xq);

        yAdj = y - fixedA;
        posMask = yAdj > 0;
        if nnz(posMask) > 1
            p = polyfit(x(posMask), log(yAdj(posMask)), 1);
            B0 = p(1);
        else
            B0 = -1;
        end

        [B, ~, residual] = lsqcurvefit(model, B0, x, y, [], [], options);
        y_fit = model(B, x_fit);
    else
        % 同时拟合 A 与 B
        model = @(p, xq) p(1) + exp(p(2) .* xq);

        A0 = min(y) - 0.1;
        yAdj = y - A0;
        posMask = yAdj > 0;
        if nnz(posMask) > 1
            p = polyfit(x(posMask), log(yAdj(posMask)), 1);
            B0 = p(1);
        else
            B0 = -1;
        end

        [params, ~, residual] = lsqcurvefit(model, [A0, B0], x, y, [], [], options);
        A = params(1);
        B = params(2);
        y_fit = model(params, x_fit);
    end

    yMean = mean(y);
    ssTotal = sum((y - yMean).^2);
    ssResidual = sum(residual.^2);
    r_squared = 1 - ssResidual / ssTotal;
end
