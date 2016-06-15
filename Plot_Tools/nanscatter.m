function h = nanscatter(inpt, nan_method, plot_method)

if nargin < 3, plot_method = 'single'; end
if nargin < 2, nan_method  = 'rowwise'; end

% Remove all NaNs from the data
if strcmp(nan_method, 'rowwise')
    mask = zeros(size(inpt));
    mask(isnan(inpt)) = NaN;
    mask = sum(mask, 2);
    inpt(isnan(mask), :) = [];
end



if strcmp(plot_method, 'single')
    h = figure
    for i = 1:size(inpt, 2) - 1
        plot(inpt(:, 1), inpt(:, i+1), '.', 'markersize', 10)
        hold on
    end
end
