function [f, g] = plotglbl(inpt_map, theta, lambda, method)


if nargin < 2, theta  = 89.75:-0.5:-89.75;  end
if nargin < 3, lambda = -179.75:0.5:179.75; end
if nargin < 4, method = 'im'; end


% Load some coastlines
load coast

% Get the size of the screen
scrsz = get(0,'ScreenSize');

% Open a new figure
f     = figure('OuterPosition',[1 scrsz(4)/4 scrsz(3)/3 scrsz(4)/2]);

% Plot the figure
if strcmp(method, 'im')
    imagesc(lambda, theta, inpt_map);
elseif strcmp(method, 'pc')
    pcolor(lambda, theta, inpt_map);
    shading interp
end

hold on

% Plot some coastlines
plot(long, lat, 'k', 'linewidth', 1.5)

% Fix ratio 
pbaspect([2 1 1])
axis xy
g = colorbar('eastoutside', 'fontsize', 14);
colormap(linspecer)

