function valid_options = set_valid_options

% -------------------------------------------------------------------------
%               Mean 2D long term error fields (full signal)
% -------------------------------------------------------------------------
valid_options.lt_2d.rmse      = 1;
valid_options.lt_2d.nse       = 1;
valid_options.lt_2d.mae       = 1;
valid_options.lt_2d.corr      = 1;
valid_options.lt_2d.pbias     = 1;

% -------------------------------------------------------------------------
%               Mean 2D long term error fields (anomalies)
% -------------------------------------------------------------------------
valid_options.lt_2d_anom.rmse   = 1;
valid_options.lt_2d_anom.nse    = 1;
valid_options.lt_2d_anom.mae    = 1;
valid_options.lt_2d_anom.corr   = 1;
valid_options.lt_2d_anom.pbias  = 1;


% -------------------------------------------------------------------------
%                 Mean 1D long term errors from 2D fields
% % -------------------------------------------------------------------------
valid_options.lt_1d.rmse      = 0;
valid_options.lt_1d.nse       = 0;
valid_options.lt_1d.mae       = 0;
valid_options.lt_1d.corr      = 0;
valid_options.lt_1d.pbias     = 1;

% % -------------------------------------------------------------------------
% %               Mean 2D seasonal error fields (full signal)
% % -------------------------------------------------------------------------
% valid_options.ssnl_2d.rmse      = 0;
% valid_options.ssnl_2d.nse       = 0;
% valid_options.ssnl_2d.mae       = 0;
% valid_options.ssnl_2d.corr      = 0;
% valid_options.ssnl_2d.pbias     = 0;



% -------------------------------------------------------------------------
% %                 Mean 1D seasonal errors (full signal)
% % -------------------------------------------------------------------------
% valid_options.rmse_1d_ssnl      = 0;
% valid_options.nse_1d_ssnl       = 0;
% valid_options.mae_1d_ssnl       = 0;
% valid_options.corr_1d_ssnl      = 0;
% valid_options.pbias_1d_ssnl     = 0;

% -------------------------------------------------------------------------
%               Mean 1D long term error fields (anomalies)
% % -------------------------------------------------------------------------
% valid_options.anom_rmse_1d_lt      = 0;
% valid_options.anom_nse_1d_lt       = 0;
% valid_options.anom_mae_1d_lt       = 0;
% valid_options.anom_corr_1d_lt      = 0;
% valid_options.anom_pbias_1d_lt     = 0;

% -------------------------------------------------------------------------
%                     Spatial averages (full signal)
% -------------------------------------------------------------------------
valid_options.spataverage.rmse      = 1;
valid_options.spataverage.nse       = 1;
valid_options.spataverage.mae       = 1;
valid_options.spataverage.corr      = 1;
valid_options.spataverage.pbias     = 1;

% -------------------------------------------------------------------------
%                     Spatial averages (anomalies)
% -------------------------------------------------------------------------
valid_options.spataverage_anom.rmse  = 1;
valid_options.spataverage_anom.nse   = 1;
valid_options.spataverage_anom.mae   = 1;
valid_options.spataverage_anom.corr  = 1;
valid_options.spataverage_anom.pbias = 1;

