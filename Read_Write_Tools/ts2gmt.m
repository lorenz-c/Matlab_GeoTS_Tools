function  [] = ts2gmt(data, fname, variable)
% Conversion from Matlab time-series to ascii-files which can be read and
% plotted by GMT.
%--------------------------------------------------------------------------
% Input:        data   [m x n] Matrix which contains one (or several)
%                              time-series. The first three (four) columns
%                              must contain (day), month, year and matlab
%                              timestamp. Data will be read from the fourth
%                              (fifth) column. 
%--------------------------------------------------------------------------
% Author: Christof Lorenz
% Date:   June 2013
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------
% Updates: 
%--------------------------------------------------------------------------

% Get the time-vector
tme = data.Data.time;

tme_out = datestr(tme, 'yyyymmdd');
tme_out = str2num(tme_out);

mat_out = NaN(length(tme_out), size(data.Data.(variable), 2) + 1);
mat_out(:, 1) = tme_out; 
mat_out(:, 2:end) = data.Data.(variable);

frmt = repmat('%g ', 1, size(data.Data.(variable), 2));
frmt = ['%i ', frmt, '\n'];

fid = fopen(fname, 'w');
fprintf(fid, frmt, mat_out');
fclose(fid);










% if nargin < 3, dtatype = 'monthly'; end
% if data(1,1) == 0, data = data(2:end, :); end
% 
% if strcmp(dtatype, 'monthly')        % Monthly Data
%     
%     tvec = data(:, 2)*10000 + data(:, 1)*100 + 15;
%     out  = [tvec data(:, 4:end)];
% 
%     frmt = repmat('%g ', 1, size(data, 2) - 3);
%     frmt = ['%i ', frmt, ' \n'];
% 
%     fid  = fopen(fname, 'w');
%     fprintf(fid, frmt, out');
%     fclose(fid);
%     
% elseif strcmp(dtatype, 'daily')      % Daily Data
%     
%     tvec = data(:, 1)*10000 + data(:, 2)*100 + data(:, 3);
%     out  = [tvec data(:, 5:end)];
% 
%     frmt = repmat('%g ', 1, size(data, 2) - 4);
%     frmt = ['%i ', frmt, ' \n'];
% 
%     fid  = fopen(fname, 'w');
%     fprintf(fid, frmt, out');
%     fclose(fid);
%     
% elseif strcmp(dtatype, 'mean_monthly')
%     data(:, 1) = 20000000+100*data(:, 1)+15;
%     
%     frmt = repmat('%g ', 1, size(data, 2)-1);
%     frmt = ['%8i ', frmt, ' \n'];
% 
%     fid  = fopen(fname, 'w');
%     fprintf(fid, frmt, data');
%     fclose(fid);
% 
% elseif strcmp(dtatype, 'annual')
%     
%     frmt = repmat('%g ', 1, size(data, 2)-1);
%     frmt = ['%i ', frmt, ' \n'];
%     
%     fid  = fopen(fname, 'w');
%     fprintf(fid, frmt, data');
%     fclose(fid)
%     
% end








