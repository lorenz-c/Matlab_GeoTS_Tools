function xx = eqm(x, y, yy, val)

% empricial Quantile-Quantile Correction
%

% Compute CDF for models during "historic" period
[f, y1] = ecdf(y);
y1 = y1(2:end);
f = f(2:end);

% Get the 
yy_rank = interp1(y1, f, yy);

xx = quantile(x, yy_rank);

end