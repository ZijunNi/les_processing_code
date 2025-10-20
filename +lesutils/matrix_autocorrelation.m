function autocorr = matrix_autocorrelation(matrix_sequence, max_lag)
%MATRIX_AUTOCORRELATION 计算矩阵序列的自相关函数。
%   A = MATRIX_AUTOCORRELATION(SEQUENCE) 会对矩阵时间序列 SEQUENCE 进行标准化
%   自相关计算。SEQUENCE 可以是尺寸为 (nRows x nCols x nTime) 的三维数组，
%   也可以是所有元素尺寸一致的元胞数组。
%
%   A = MATRIX_AUTOCORRELATION(SEQUENCE, MAX_LAG) 将计算限制在给定的最大时间延迟 MAX_LAG 内。

    if nargin < 2 || isempty(max_lag)
        max_lag = [];
    end

    if iscell(matrix_sequence)
        T = numel(matrix_sequence);
        if T == 0
            error('matrix_autocorrelation:EmptyInput', ...
                  'Input matrix sequence is empty.');
        end
        [n, m] = size(matrix_sequence{1});
        data3d = zeros(n, m, T);
        for tIdx = 1:T
            current = matrix_sequence{tIdx};
            if ~isequal(size(current), [n, m])
                error('matrix_autocorrelation:InconsistentSizes', ...
                      'All matrices in the sequence must have identical size.');
            end
            data3d(:, :, tIdx) = current;
        end
    else
        if ndims(matrix_sequence) ~= 3
            error('matrix_autocorrelation:InvalidDimensions', ...
                  'Numeric input must be an n-by-m-by-T array.');
        end
        data3d = matrix_sequence;
        [n, m, T] = size(data3d);
    end

    if isempty(max_lag)
        max_lag = T - 1;
    else
        max_lag = min(max_lag, T - 1);
    end

    if max_lag < 0
        error('matrix_autocorrelation:InvalidLag', ...
              'max_lag must be non-negative and less than the sequence length.');
    end

    % 将矩阵序列向量化以便计算相关系数
    sequence2d = reshape(data3d, n * m, T).';

    meanVec = mean(sequence2d, 1);
    centered = sequence2d - meanVec;
    denom = sum(centered(:).^2);

    autocorr = zeros(1, max_lag + 1);
    for tau = 0:max_lag
        if tau == 0
            numerator = denom;
        else
            numerator = sum(sum(centered(1:end-tau, :) .* centered(1+tau:end, :)));
        end
        autocorr(tau + 1) = numerator / denom;
    end
end
