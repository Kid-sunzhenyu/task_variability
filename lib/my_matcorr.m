function C = my_matcorr( s_series, t_series )
%MY_MATCORR compute diagnoal correlation between two matrices with the same size
% Input: s_series is T x N, t_series is T x N; Size must be the same
% Output: C is vector = diag(CoRRMat)

% edited by Jianxun Ren, 20171110
    s_series = bsxfun(@minus, s_series, mean(s_series, 1));
    s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));

    t_series = bsxfun(@minus, t_series, mean(t_series, 1));
    t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));
    
    C = sum(bsxfun(@times, s_series, t_series));
end

