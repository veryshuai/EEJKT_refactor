function iter_out = simulateForeignMatches(pt_ndx,macro_state_f,mm,policy)

[iter_in, iter_out] = simulateForeignInnerInitialize(mm, pt_ndx);

iter_in.year = 1;
%tlag = 1;
season = 1;
iter_in.N_match = 0;
iter_in.firm_yr_sales_lag = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),4);
% firm_yr_sales_lag will contain: [firmID,sales,#shipments,firm age]
tic
my_flag = 0;
for t = 2:1:mm.periods

    iter_in.t = t;

    % update year, type and pmat_cum_t
    iter_in.year = floor((iter_in.t-1)/mm.pd_per_yr);
    iter_in.mic_type = find(policy.firm_type_prod_succ_macro(:,3) == mm.pt_type(pt_ndx,2) & policy.firm_type_prod_succ_macro(:,4) == mm.pt_type(pt_ndx,1),1,'first');
    type     = find(policy.firm_type_prod_succ_macro(:,2) == macro_state_f(iter_in.t) & policy.firm_type_prod_succ_macro(:,3) == mm.pt_type(pt_ndx,2) & policy.firm_type_prod_succ_macro(:,4) == mm.pt_type(pt_ndx,1),1,'first');

    % Load relevant transition probabilites for current micro type & macro state.
    % pmat_cum_t holds cumulative transition probs across (#success, #meeting) pairs.
    % It's created in inten_sim using the relevant Q matrix.
    pmat_cum_t = policy.pmat_cum_f{type}; % type-specific cum trans. probs., micro state

    % reset season when previous period completes a year
    if abs(floor((iter_in.t-1)/mm.pd_per_yr) - (iter_in.t-1)/mm.pd_per_yr) <= 0.0001
        season = 1;
    end

    %% gross additions to clients, before drops and deaths between t-1 and t

    %  Need to treat firms that have maxed out learning separately from others.
    %  To get to the no-learning state, firms must first land on or above n = mm.n_size-3.
    %  where mm.n_size is the maximum number of meetings that firms learn from.
    %  The -3 keeps firms transiting to the no-learning states from first having to
    %  land first on mm.n_size. They can land on mm.n_size,mm.n_size-1,mm.n_size-2, or mm.n_size-3
    stay     = ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1);  % will flag no-learning firms that continue, t-1 to t
    no_learn = iter_in.cum_meets(:,iter_in.t-1) >= mm.n_size-3; % for picking off seasoned exporters
    learn    = iter_in.cum_meets(:,iter_in.t-1) <  mm.n_size-3; % for picking off learning exporters

    % First deal with firms that have not maxed out learning using pmat_cum, based on Q matrix.

    % Find N_learn randomly selected micro states next period, given
    % macro state (common to all firms), initial micro states, and pmat_cum.
    N_learn = sum(learn,1) ;
    if  N_learn > 0
        trans_rands = pmat_cum_t(iter_in.micro_state(learn,iter_in.t-1),:)> rand(N_learn,1)*ones(1,size(policy.pmat_cum_f{1},2));
        iter_in.micro_state(learn,iter_in.t) = int16(size(policy.pmat_cum_f{1},2) + 1 - sum(trans_rands,2)); % drawn new micro states
        iter_in.cum_meets(learn,iter_in.t)   = policy.pmat_to_meets_succs(iter_in.micro_state(learn ,iter_in.t),2) - 1; % trials, new state, matrix for all firms (iter_in.t+1)
        iter_in.cum_succ(learn,iter_in.t)    = policy.pmat_to_meets_succs(iter_in.micro_state(learn ,iter_in.t),3) - 1; % successes, new state, matrix for all firms starting (t+1)
    end

    stay(learn)  = iter_in.micro_state(learn,iter_in.t-1) - iter_in.micro_state(learn,1) ~= 0; % wasn't in initial state last period

    % Now deal with firms that have maxed-out learning.

    % Update cumulative meetings and successes for firms with >= mm.n_size meets, if
    % any. These firms are presumed to know their thetas, hence nn1 = floor(succ_prob*nn2)
    % These calculations don't involve the Q matrix or the associated pmat trans. probs.
    N_no_learn  = sum(no_learn);
    if N_no_learn >0
        % no exog. death, and shipments prev. year
        stay(no_learn) = (rand(N_no_learn,1) > 1-exp(-mm.firm_death_haz));

        % CHECK THIS LINE: arguments of policy.lambda_f right?
        iter_in.cum_meets(no_learn,iter_in.t) = (iter_in.cum_meets(no_learn,iter_in.t-1) + poissinv(rand(N_no_learn,1),mm.theta2(mm.pt_type(pt_ndx,2)) ...
            * reshape(policy.lambda_f(floor(mm.theta2(mm.pt_type(pt_ndx,2))*(mm.n_size+1)),(mm.n_size+1),1,min((mm.net_size+1),iter_in.cum_succ(no_learn,iter_in.t-1)+1),...
            mm.pt_type(pt_ndx,1),macro_state_f(iter_in.t-1)),N_no_learn,1))).*stay(no_learn); % resets to 0 if stay==0 or new firm this period

        iter_in.cum_succ(no_learn,iter_in.t) = iter_in.cum_succ(no_learn,iter_in.t-1).*stay(no_learn)...  %.*(1-iter_in.new_firm(no_learn)) ...
            + random('bino',iter_in.cum_meets(no_learn,iter_in.t)-iter_in.cum_meets(no_learn,iter_in.t-1).*stay(no_learn),mm.theta2(mm.pt_type(pt_ndx,2)));
        % Note: cum_succ includes single shipment matches (with bad z's) that are dropped after 1 period.
        % cum_succ and cum_meets will always be 0 after a period t reset due to exit (stay==0).
    end

    %(succ,trial,common succ rate (defunct), network size, prod of firm, macro shock)

    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    % Hereafter can treat learning and no-learning firms together

    iter_in.add_cli_cnt(:,iter_in.t) = max(iter_in.cum_succ(:,iter_in.t) - iter_in.cum_succ(:,iter_in.t-1),0); %
    % max() resets count to 0 for neg. differences to deal with new exporters

    % identify first period in which a new exporter is active. First
    % condition is for learning firms; second (stay==0) is for no-learning

    iter_in.new_firm(:,iter_in.t)    = max((iter_in.cum_meets(:,iter_in.t) - iter_in.cum_meets(:,iter_in.t-1) < 0),(stay == 0));
    iter_in.exit_firm(:,iter_in.t-1) = iter_in.cum_meets(:,iter_in.t) - iter_in.cum_meets(:,iter_in.t-1) < 0 ; % last period before exit

    % NOTE: iter_in.new_firm = 1 each period that stay = 0.

    incumb        = ones(size(iter_in.new_firm(:,iter_in.t),1),1)- iter_in.new_firm(:,iter_in.t);
    % Careful: incumb means no reset after previous period.  It is not the same
    % as cont_expr, which tracks whether an exporter had at least one
    % shipment during the previous year.


    %% Deal with endogenous and exogenous drops

    % identify z values at which exporters keep current matches from t-1 to t
    iter_in.keep_cli = policy.c_val_f(:,mm.pt_type(pt_ndx,1),macro_state_f(iter_in.t-1))' > 0; % = 1 if want to keep type for t
    drop_Zcut = size(mm.Z,1) - sum(iter_in.keep_cli); % cutoff: matches dropped at z value <= drop_Zcut

    % count endogenous drops for all exporter hotel rooms (exporters): z too low to continue
    drop_cnt = sum(iter_in.lag_cli_zst.*(1-iter_in.keep_cli),2);

    % draw the number of exogenous deaths of remaining matches between t-1 and t
    ddum = find(iter_in.cur_cli_cnt(:,iter_in.t-1)-drop_cnt > 0);
    if sum(ddum)>0
        iter_in.exog_deaths(ddum,iter_in.t-1) =...
            random('bino',iter_in.cur_cli_cnt(ddum,iter_in.t-1)-drop_cnt(ddum),1-exp(-mm.delta));
    end

    %% update current count for new matches, drops, and exogenous deaths
    iter_in.cur_cli_cnt(:,iter_in.t) = iter_in.add_cli_cnt(:,iter_in.t) + iter_in.cur_cli_cnt(:,iter_in.t-1) ...
        - drop_cnt - iter_in.exog_deaths(:,iter_in.t-1) ;

    %% break down by buyer types (z)

    for i=1:mm.sim_firm_num_by_prod_succ_type(pt_ndx)
        % break down new clients that occur between t-1 and t into e.o.p. z-types
        iter_in.new_cli_zst(i,:) = new_vec_C(iter_in.add_cli_cnt(i,iter_in.t),size(mm.Z,1),cumsum(mm.erg_pz)); % distribute gross additions
        if sum(iter_in.new_cli_zst(i,:)) ~= iter_in.add_cli_cnt(i,iter_in.t) % cheap patch--better to clean up C++ code
            'C routine for random draws failed in matchdat_gen_f. Trying Matlab version'
            iter_in.new_cli_zst(i,:) = new_vec(iter_in.add_cli_cnt(i,iter_in.t),size(mm.Z,1),cumsum(mm.erg_pz));
        end
        if iter_in.exog_deaths(i,iter_in.t-1) > 0
            % break down exogenous deaths that occur between t-1 and t down by b.o.p. z state:
            iter_in.die_cli_zst(i,:) = die_vec(iter_in.lag_cli_zst(i,:).*iter_in.keep_cli,iter_in.exog_deaths(i,iter_in.t-1),size(mm.Z,1));
        end
        %trans_count_test = trans_count;
        iter_in.trans_count(2:size(mm.Z,1)+1,1,i) = (iter_in.lag_cli_zst(i,:).*(1-iter_in.keep_cli))' + iter_in.die_cli_zst(i,:)';
        % For each firm (i) of a particular type, column 1 of trans_count(:,:,i)
        % now contains counts of all exiting matches (endog. and exog.), by buyer type (row).

        % Update surviving client counts by z type using transition matrix for z.
        % Do this for those that don't die for endogenous or exogenous reasons.
        iter_in.surviv_zst(i,:) = iter_in.lag_cli_zst(i,:).*iter_in.keep_cli - iter_in.die_cli_zst(i,:);
        N_sur = sum(iter_in.surviv_zst(i,:),2); % number of survivors from t-1 by b.o.p. type, firm i

        if N_sur > 0
            sur_typ = find(iter_in.surviv_zst(i,:)); % addresses for z-states populated by at least one survivor, firm i
            for jj = sur_typ  % loop over initial (b.o.p.) states of surviving matches, exporter i
                draw = rand(iter_in.surviv_zst(i,jj),1);
                % identify destination z states for each surviving client (could be multiple survivors per initial type):
                trans_z = ones(iter_in.surviv_zst(i,jj),1)*policy.pmat_cum_z(jj,:) > draw;
                % count # clients in each destination z state for each beginning z state.
                % Record e.o.p. counts in cols 2:N_Z+1 of trans_count. Row indices are b.o.p. states, plus 1:
                iter_in.trans_count(jj+1,2:size(mm.Z,1)+1,i) = sum(trans_z(:,1:size(mm.Z,1)) - [zeros(size(draw,1),1),trans_z(:,1:size(mm.Z,1)-1)],1);
                % cumulate over b.o.p. z types to get row vector of surviving client e.o.p. types. Rows (i) index exporter hotel rooms:
                iter_in.trans_zst(i,:) = iter_in.trans_zst(i,:) +  iter_in.trans_count(jj+1,2:size(mm.Z,1)+1,i);
            end
        end
        if sum(iter_in.new_cli_zst(i,:),2)>0
            iter_in.trans_count(1,2:size(mm.Z,1)+1,i) = iter_in.new_cli_zst(i,:); % load new clients for exporter i in first row of trans_count
        end



    end
    cur_cli_zst = iter_in.new_cli_zst + iter_in.trans_zst;


    %% Deal with dormant firms: if a firm has no clients for more than one year,
    %  kick it out of its slot and start a new firm

    if iter_in.t > mm.pd_per_yr
        dmt = sum(iter_in.cur_cli_cnt(:,iter_in.t-mm.pd_per_yr+1:iter_in.t),2)==0; % identify dormant firms
        iter_in.micro_state(dmt,iter_in.t)  = 1;  % reset initial micro state to 1 (entrant)
        iter_in.new_firm(dmt,iter_in.t)     = 1;  % will mark new firms that haven'iter_in.t made a match
        iter_in.exit_firm(dmt,iter_in.t)    = 1;  % will mark last period of exiting firm (same thing)
        iter_in.cum_meets(dmt,iter_in.t)    = 0;  % cumulative number of meetings
        iter_in.cum_succ(dmt,iter_in.t)     = 0;  % cumulative number of successes
    end

    % Calculate time (in periods) in export market
    flr       = max(iter_in.flrlag,iter_in.t*iter_in.new_firm(:,iter_in.t)); % floor resets to current year for new exporters.
    age       = iter_in.t*ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1) - flr;     % age in periods. age=0 all year for firms with no shipment previous year
    iter_in.flrlag    = flr ; % carry floor forward for continuing matches
    iter_in.cumage    = cat(2,iter_in.cumage,age);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Construct period-specific variables

    %  Load season to season transitions into mat_tran, which describes matches
    %  of all mm.sim_firm_num_by_prod_succ_type(pt_ndx) of a particular type for a particular transition (t-1 to t).

    mat_tran_all_zeros = ~any(iter_in.trans_count(:));
    if mat_tran_all_zeros
        mat_tran = zeros(0,4);ship_cur = zeros(0,1); age_vec = zeros(0,1);
    else

        mkt =1; % =1 for foreign market
        [mat_tran,ship_cur,age_vec] = match_sales(mkt,mm,iter_in.trans_count,age,pt_ndx,macro_state_f(iter_in.t));

        % [mat_tran,ship_cur,age_vec] =...
        %     match_sales(mm.scale_f,mm.eta,iter_in.trans_count(:,:,:),age,mm.X_f(macro_state_f(t)),t,mm.poisCDF_shipments,...
        %     mm.max_ships,size(mm.Z,1),mm.Z,mm.Phi(mm.pt_type(pt_ndx,1)));

        % mat_tran:  [initial state, exporter id, ending state, match revenue]
        % ship_cur:   match's number of shipments within the current period
        % age_vec:    firm age (# periods)
    end

    if size(mat_tran,1)>0
        try
            % check that ending Z, if positive, is always greater than Zcut.
            assert(min(  (mat_tran(:,3)>0).*((mat_tran(:,3)>0) - drop_Zcut.*ones(size(mat_tran,1),1) )>0,[],1)>=0)
        catch
            warning('ending Z positive and less than Zcut')
            [mat_tran(:,3),drop_Zcut.*ones(size(mat_tran,1),1)];
        end
        if season > 1
            try
                % check that beginning Z, if positive, is always greater than last period's Zcut.
                assert(min((mat_tran(:,1)>0).*((mat_tran(:,1)>0) - iter_in.seas_Zcut(season-1).*ones(size(mat_tran,1),1))>0,[],1)>=0)
            catch
                warning('beginning Z positive and less than last period Zcut')
                [mat_tran(:,1),iter_in.seas_Zcut(season-1).*ones(size(mat_tran,1),1)]
            end
        end
    end

    if season == 1
        iter_in.N_match = size(mat_tran,1);
    end

    %if t==2
    %    seas_Zcut_lag = iter_in.seas_Zcut;
    %    seas_tran_lag = iter_in.seas_tran;
    %end % needed to get started (will be discarded with the burn-in)

    % Concatenate time index, season index and year index with match variables
    % and collect results for all seasons in year in seas_tran
    iter_in.seas_tran{1,season} = [[iter_in.t,season,iter_in.year].*ones(size(mat_tran,1),1),mat_tran,ship_cur,age_vec];
    % iter_in.seas_tran: [t, season, year, initial state, exporter id, ending state, match revenue,
    %             # shipments, firm age (# periods)]
    iter_in.seas_Zcut(season)   = drop_Zcut;

    iter_in.N_match = iter_in.N_match + size(mat_tran,1); % cumulate match count within current year

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% construct annualized variables
    if season == mm.pd_per_yr

        [iter_in,iter_out] = simulateForeignInnerAnnualize(iter_in,iter_out,mm);

    end   % season == mm.pd_per_yr if statement
    

    season = season + 1;

    %% load lagged client state matrix and re-initialize objects

    iter_in.lag_cli_zst  = cur_cli_zst;
    iter_in.new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iter_in.die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iter_in.trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iter_in.trans_count  = zeros(size(mm.Z,1)+1,size(mm.Z,1)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx));

end     % end of time loop
end
