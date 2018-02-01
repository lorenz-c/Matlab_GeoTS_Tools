
mdl_p_crr_qm  = NaN(length(mdl_p), 1);
mdl_p_crr_qdm = NaN(length(mdl_p), 1);
mdl_p_crr_dqm = NaN(length(mdl_p), 1);

for i = 2:length(tme_p)-1
    
    mnth_act = tme_p(i, 2);
    yr_act   = tme_p(i, 1);
    
    if mnth_act == 1
        mnth_bw = 12;
        yr_bw   = yr_act - 1;
    else
        mnth_bw = mnth_act - 1;
        yr_bw   = yr_act;
    end
    
    if mnth_act == 12
        mnth_fw = 1;
        yr_fw   = yr_act + 1;
    else
        mnth_fw = mnth_act + 1;
        yr_fw   = yr_act;
        
    end
    
    
    mnth_bw  = mnth_act - 1;
    
    ids_h_bw  = find(tme_h(:, 2) == mnth_bw);
    ids_h_act = find(tme_h(:, 2) == mnth_act);
    ids_h_fw  = find(tme_h(:, 2) == mnth_fw);
    
    ids_p_bw  = find(tme_p(:, 2) == mnth_bw  & tme_p(:, 1) == yr_bw);
    ids_p_act = find(tme_p(:, 2) == mnth_act & tme_p(:, 1) == yr_act);
    ids_p_fw  = find(tme_p(:, 2) == mnth_fw  & tme_p(:, 1) == yr_fw);
    
    mdl_h_mnth = mdl_h([ids_h_bw; ids_h_act; ids_h_fw]);
    obs_h_mnth = obs_h([ids_h_bw; ids_h_act; ids_h_fw]);
    
    mdl_p_mnth = mdl_p([ids_p_bw; ids_p_act; ids_p_fw]);
    obs_p_mnth = obs_p([ids_p_bw; ids_p_act; ids_p_fw]);
    
    mdl_p_crr_qm([ids_p_bw; ids_p_act; ids_p_fw]) = bias_correction(mdl_h_mnth, obs_h_mnth, mdl_p_mnth, 'qm', 'ecdf', 0.05);
    mdl_p_crr_dqm([ids_p_bw; ids_p_act; ids_p_fw]) = bias_correction(mdl_h_mnth, obs_h_mnth, mdl_p_mnth, 'dqm', 'ecdf', 0.05);
    mdl_p_crr_qdm([ids_p_bw; ids_p_act; ids_p_fw]) = bias_correction(mdl_h_mnth, obs_h_mnth, mdl_p_mnth, 'qdm', 'ecdf', 0.05);
end


    
    
    