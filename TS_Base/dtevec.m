function [DateTime, TimeStamp] = dtevec(sdte, edte, tres)
%--------------------------------------------------------------------------
% Create a date-time matrix and the appropriate numeric dates for a given
% time period, which is defined by its start- and end-date. The function
% allows different "styles" of input dates and the resolution of the output
% values is defined by the tres-parameter.
%--------------------------------------------------------------------------
% INPUT:
% - sdte      Start-date of the time-period. The lenght of this vector
%             depends on the chosen temporal resolution. For e.g. daily
%             data, sdte could be [2000 01 01]. For monthly data, the day
%             is neglected, i.e. sdte = [2000 01 01] and sdte = [2000 01]
%             gives the same result. 
%             If not a "full" date is provided, the start month/day/hour 
%             is set to 01/01/00, depending on the chosen temporal
%             resolution.
% - edte      End-date. See description of sdte
% - tres      Resolution of the output date-time. Can be set to 'hourly',
%             'daily', 'monthly' (days are set to 15), 'monthly_alt' (days 
%             are set to 1), 'seasonal'
%--------------------------------------------------------------------------
% OUTPUT:
% - DateTime  Date-time matrix with YYYY MM DD HH MM SS
% - TimeStamp Numeric date (from Matlab's datenum-command)
%--------------------------------------------------------------------------
% EXAMPLE:
% >> [DateTime, TimeStamp] = doy2date([2000 05], 2004, 'daily');
% --> The function creates the dates from 2000-05-01 to 2004-12-31.
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         October 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------
if nargin < 3, tres = 'monthly'; end

if strcmp(tres, 'hourly') | strcmp(tres, '1 hour')
    
    if length(sdte) >= 4
        % yyyy/mm/dd/hh
        sdte_num = datenum([sdte(1:4) 0 0]);
    elseif length(sdte) == 3
        % yyyy/mm/dd
        sdte_num = datenum([sdte 0 0 0]);
    elseif length(sdte) == 2
        % yyyy/mm
        sdte_num = datenum([sdte 1 0 0 0]); 
    elseif length(sdte) == 1
        % yyyy
        sdte_num = datenum([sdte 1 1 0 0 0]); 
    end
    
    if length(edte) >= 4
        % yyyy/mm/dd/hh
        edte_num = datenum([edte(1:4) 0 0]);
    elseif length(sdte) == 3
        % yyyy/mm/dd
        edte_num = datenum([edte 23 0 0]);
    elseif length(sdte) == 2
        % yyyy/mm
        edte_num = datenum([edte eomday(edte) 23 0 0]); 
    elseif length(sdte) == 1
        % yyyy
        edte_num = datenum([edte 12 31 23 0 0]); 
    end
    
    TimeStamp = (sdte_num:1/24:edte_num)';
    DateTime  = datevec(TimeStamp);
    
    
elseif strcmp(tres, 'daily') | strcmp(tres, '1 day')
    
    if length(sdte) >= 3
        % yyyy/mm/dd
        sdte_num = datenum([sdte(1:3) 0 0 0]);
    elseif length(sdte) == 2
        % yyyy/mm
        sdte_num = datenum([sdte 1 0 0 0]); 
    elseif length(sdte) == 1
        % yyyy
        sdte_num = datenum([sdte 1 1 0 0 0]); 
    end
    
    if length(edte) >= 3
        % yyyy/mm/dd
        edte_num = datenum([edte(1:3) 0 0 0]);
    elseif length(edte) == 2
        % yyyy/mm
        edte_num = datenum([edte eomday(edte(1), edte(2)) 0 0 0]); 
    elseif length(edte) == 1
        % yyyy
        edte_num = datenum([edte 12 31 0 0 0]); 
    end
    
    TimeStamp = (sdte_num:1:edte_num)';
    DateTime  = datevec(TimeStamp);
    
elseif strcmp(tres, 'monthly_alt')
    
    if length(sdte) >= 2
        % yyyy/mm
        sdte_num = datenum([sdte(1:2) 1 0 0 0]); 
    elseif length(sdte) == 1
        % yyyy
        sdte_num = datenum([sdte 1 1 0 0 ]); 
    end
    
    if length(edte) >= 2
        % yyyy/mm
        edte_num = datenum([edte(1:2) 1 0 0 0]); 
    elseif length(sdte) == 1
        % yyyy
        edte_num = datenum([edte 12 1 0 0 0]); 
    end
    
    TimeStamp = (sdte_num:1:edte_num)';
    DateTime  = datevec(TimeStamp);
    
    TimeStamp = TimeStamp(DateTime(:, 3) == 1, :);
    DateTime  = DateTime(DateTime(:, 3) == 1, :);
    
elseif strcmp(tres, 'monthly') | strcmp(tres, '1 month')
    
    if length(sdte) >= 2
        % yyyy/mm
        sdte_num = datenum([sdte(1:2) 15 0 0 0]); 
    elseif length(sdte) == 1
        % yyyy
        sdte_num = datenum([sdte 1 15 0 0 0]); 
    end
    
    if length(edte) >= 2
        % yyyy/mm
        edte_num = datenum([edte(1:2) 15 0 0 0]); 
    elseif length(sdte) == 1
        % yyyy
        edte_num = datenum([edte 12 15 0 0 0]); 
    end
    
    TimeStamp = (sdte_num:1:edte_num)';
    DateTime  = datevec(TimeStamp);
    
    TimeStamp = TimeStamp(DateTime(:, 3) == 15, :);
    DateTime  = DateTime(DateTime(:, 3) == 15, :);
    
elseif strcmp(tres, 'seasonal') 
    
    if length(sdte) >= 2
        % yyyy/mm
        if ~ismember(sdte(2), [1 4 7 10])
            error('For seasonal time-scales, the first month has to be 1, 4, 7, or 10!')
        else
            sdte_num = datenum([sdte(1:2) 15 0 0 0]); 
        end 
    elseif length(sdte) == 1
        % yyyy
        sdte_num = datenum([sdte 1 15 0 0 0]); 
    end
    
    if length(edte) >= 2
        % yyyy/mm
        if ~ismember(edte(2), [1 4 7 10])
            error('For seasonal time-scales, the last month has to be 1, 4, 7, or 10!')
        else
            edte_num = datenum([edte(1:2) 15 0 0 0]); 
        end 
    elseif length(edte) == 1
        % yyyy
        edte_num = datenum([edte 10 15 0 0 0]); 
    end
    
    TimeStamp = (sdte_num:1:edte_num)';
    DateTime  = datevec(TimeStamp);
    
    TimeStamp = TimeStamp(DateTime(:, 3) == 15 & ismember(DateTime(:, 2), [1 4 7 10]));
                      
    DateTime  = datevec(TimeStamp);
        
elseif strcmp(tres, 'annual') | strcmp(tres, '1 year') | strcmp(tres, 'yearly')
    
    yrs       = (sdte(1):edte(1))';
    DateTime  = [yrs ones(length(yrs), 2) zeros(length(yrs), 3)];
    TimeStamp = datenum(DateTime);
            
end


    
    


