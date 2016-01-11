function [f, varargout] = tsplot(nrplots, period, indxs, datatype, yy, ttle, varargin)
% Simple function for plotting time-series. It can display a single or
% multiple time-series, which would be arranged according to nrplots (e.g.
% nrplots = [3 2] would create a figure "matrix" with three columns and two
% rows). Period control the time-steps which will be drawn. With indxs, the
% user can choose specific time-series in every prodivded data-set. With
% datatype, one can select between daily, monthly, yearly, or mean monthly
% input data. Finally, varargin contains all the datasets which should be
% plotted.
%--------------------------------------------------------------------------
% Input:        nrplots   [m x n]   Scalar or matrix controling the number
%                                   and the arrangement of plots
%               period    [1 x 2]   Start- and end-date of the period which
%                                   should be plotted
%               indxs     [1 x l]   Indices of the desired time-series
%               datatype  'string'  One of the following:
%                                   - 'daily'
%                                   - 'monthly' (default)
%                                   - 'yearly'
%                                   - 'monthly_mn'
%               varargin            One or more matrices containing the
%                                   different time-series. 
% Output:       h                   Figure identifier
%--------------------------------------------------------------------------
% Author: Christof Lorenz
% Date:   July 2011
%--------------------------------------------------------------------------
% Uses: linspecer.m (Mathworks)
%--------------------------------------------------------------------------
% Updates: 
%--------------------------------------------------------------------------


if isempty(datatype), datatype = 'monthly'; end
if isempty(yy),       yy       = 0;         end

if strcmp(datatype, 'daily')
    tmevec = 4;  
elseif strcmp(datatype, 'monthly')
    tmevec = 3;
elseif strcmp(datatype, 'annual') | strcmp(datatype, 'monthly_mn')
    tmevec = 1;
end

for i = 1:length(varargin)
    for j = 1:length(indxs)
%         if strcmp(datatype, 'daily') | strcmp(datatype, 'monthly')
%             varargin{i} = findtstps_ts(varargin{i}, [period(1) 1], [period(2) 12]);
%         end
        ccol(i,j)   = find(varargin{i}(1,:) == indxs(j));
    end
end

if strcmp(datatype, 'monthly_mn'), mnths = mnthnms('vshort'); end
   
clr = linspecer(length(varargin));

scrsz = get(0,'ScreenSize');

f = figure('OuterPosition',[1 scrsz(4)/4 scrsz(3)/3 scrsz(4)/2]);

if yy == 1
    for i = 1:nrplots(1)*nrplots(2)
        h(i) = subplot(nrplots(1), nrplots(2), i);
        
        [AX,H1,H2] = plotyy(varargin{1}(2:end, tmevec), ...
                            varargin{1}(2:end, ccol(1, i)), ...
                            varargin{2}(2:end, tmevec), ...
                            varargin{2}(2:end, ccol(2, i)));
           
        if strcmp(datatype, 'daily') | strcmp(datatype, 'monthly')
           datetick(AX(1),'x','yyyy');
           set(AX(1), 'xtick', 0:1e10:10e10);
           datetick(AX(2),'x','yyyy');
           
           set(AX(1),'Xlim', [datenum(period(1), 1, 15) ...
                              datenum(period(2), 12, 15)]);
           set(AX(2),'Xlim', [datenum(period(1), 1, 15) ...
                              datenum(period(2), 12, 15)]);
            
            
        elseif strcmp(datatype, 'monthly_mn')
            set(AX(1),'xlim',[1 12]);
            set(AX(1), 'xticklabel', mnths); 
            set(AX(2),'XTick',[]);
        elseif strcmp(datatype, 'annual')
            set(AX(2),'XTick',[]);
            set(AX(1),'Xlim', [period(1) period(2)]);
            set(AX(2),'Xlim', [period(1) period(2)]);
        end
    end    
else
    for i = 1:nrplots(1)*nrplots(2)
        h(i) = subplot(nrplots(1), nrplots(2), i);

        for j = 1:length(varargin)
            plot(varargin{j}(2:end, tmevec), varargin{j}(2:end, ccol(j, i)),  'Color', clr(j,:), 'Linewidth', 1.7);
            hold on
        end
    
        if strcmp(datatype, 'daily') | strcmp(datatype, 'monthly') 
            datetick('x');
            set(gca,'Xlim', [datenum(period(1), 1, 15) ...
                             datenum(period(2), 12, 15)]);
        elseif strcmp(datatype, 'monthly_mn')
            set(gca, 'xlim',  [1 12]);
            set(gca, 'xticklabel', mnths);
        end
        
        if ~isempty(ttle)
            title(ttle{i}, 'fontsize', 14);
        end
        
    end
    
end


set(h, 'Fontsize', 12);
varargout{1} = h;
