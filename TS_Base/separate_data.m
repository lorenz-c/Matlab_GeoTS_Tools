function [ts_out, tme_vec, num_vec, indx_vec] = separate_data(ts_in, tres_in, indx_row)


if isstruct(ts_in)
    nts = size(ts_in.Time, 1);
    
    tme_vec = ts_in.Time;
    
    if isfield(ts_in, 'Timestamp')
        num_vec = ts_in.Timestamp;
    else
        num_vec = [];
    end
    
    if indx_row == true
        indx_vec = ts_in.Index;
    end
   
else
    
    if strcmp(tres_in, 'seconds')       % ss mm hh dd mm yy  
        tme_vec = ts_in(:, 1:6);
        num_vec = ts_in(:, 7);
        ts_out  = ts_in(:, 8:end);
    elseif strcmp(tres_in, 'hourly')    % hh dd mm yy  
        tme_vec = ts_in(:, 1:4);
        num_vec = ts_in(:, 5);
        ts_out  = ts_in(:, 6:end);
    elseif strcmp(tres_in, 'daily')     % dd mm yy  
        tme_vec = ts_in(:, 1:3);
        num_vec = ts_in(:, 4);
        ts_out  = ts_in(:, 5:end);
    elseif strcmp(tres_in, 'monthly')   % mm yy  
        tme_vec = ts_in(:, 1:2);
        num_vec = ts_in(:, 3);
        ts_out  = ts_in(:, 4:end);
    elseif strcmp(tres_in, 'annual')    % yy  
        tme_vec = ts_in(:, 1);
        num_vec = [];
        ts_out  = ts_in(:, 2:end);
    elseif strcmp(tres_in, 'none')
        tme_vec = (1:nts)';
        num_vec = [];
    end

    if indx_row == true
        indx_vec = ts_out(1, :);
        ts_out   = ts_out(2:end, :);
    
        if strcmp(tres_in, 'none')
            tme_vec = tme_vec(1:end-1);
        else
            tme_vec = tme_vec(2:end, :);
        end
    
        if ~isempty('num_vec')
            num_vec = num_vec(2:end, :);
        end 
    else
        indx_vec = [];
    end
    
end

