function corr_mat = my_corr(s_series, t_series)

% Input: s_series is T x N1, t_series is T x N2
% Output: corr_mat is N1 x N2
%
% Input, Output are identical to build-in CORR,
% but my_corr is faster 10 more times than it
%
% created by B.T.Thomas Yeo
% Combine my_corr and self_corr togather.
% edited by Jianxun Ren, 20171110

if nargin < 2
    t_series = s_series;
end

s_series = bsxfun(@minus, s_series, mean(s_series, 1));
s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));

t_series = bsxfun(@minus, t_series, mean(t_series, 1));
t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));

corr_mat = s_series' * t_series;