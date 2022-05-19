function iter_in = simulateForeignMatchesInnerSim(iter_in,mm,policy)

    %  Need to treat firms that have maxed out learning separately from others.
    %  To get to the no-learning state, firms must first land on or above n = mm.n_size-3.
    %  where mm.n_size is the maximum number of meetings that firms learn from.
    %  The -3 keeps firms transiting to the no-learning states from first having to
    %  land first on mm.n_size. They can land on mm.n_size,mm.n_size-1,mm.n_size-2, or mm.n_size-3
    no_learn = iter_in.cum_meets(:,iter_in.t-1) >= mm.n_size-3; % for picking off seasoned exporters
    learn    = iter_in.cum_meets(:,iter_in.t-1) <  mm.n_size-3; % for picking off learning exporters
    stay     = ones(mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx),1);  % will flag no-learning firms that continue, t-1 to t

    [iter_in, stay] = simulateForeignMatchesInnerSimLearners(learn, iter_in, policy, stay);

    [stay, iter_in] = simulateForeignMatchesInnerSimNoLearners(no_learn, stay, mm, iter_in, policy);

    iter_in.add_cli_cnt(:,iter_in.t) = max(iter_in.cum_succ(:,iter_in.t) - iter_in.cum_succ(:,iter_in.t-1),0); %max() resets count to 0 for neg. differences to deal with new exporter
    iter_in.new_firm(:,iter_in.t)    = max((iter_in.cum_meets(:,iter_in.t) - iter_in.cum_meets(:,iter_in.t-1) < 0),(stay == 0)); %first period for new exporter
    iter_in.exit_firm(:,iter_in.t-1) = iter_in.cum_meets(:,iter_in.t) - iter_in.cum_meets(:,iter_in.t-1) < 0 ; % last period before exit
    
    % DAVID QUESTION: Does exit firm pick up no learning firms that exit?
    % If so, why is stay==0 condition needed?

    [iter_in, drop_Zcut, drop_cnt] = simulateForeignMatchesInnerSimDrops(iter_in, policy, mm);

    iter_in.cur_cli_cnt(:,iter_in.t) = iter_in.add_cli_cnt(:,iter_in.t) + iter_in.cur_cli_cnt(:,iter_in.t-1) ...
        - drop_cnt - iter_in.exog_deaths(:,iter_in.t-1) ;

    %% break down by buyer types (z)

    for i=1:mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx)
        % break down new clients that occur between t-1 and t into e.o.p. z-types
        iter_in.new_cli_zst(i,:) = new_vec_C(iter_in.add_cli_cnt(i,iter_in.t),size(mm.Z,1),cumsum(mm.erg_pz)); % distribute gross additions
        if sum(iter_in.new_cli_zst(i,:)) ~= iter_in.add_cli_cnt(i,iter_in.t) % cheap patch--better to clean up C++ code
            'C routine for random draws failed in matchdat_gen_f. Trying Matlab version'
            iter_in.new_cli_zst(i,:) = new_vec(iter_in.add_cli_cnt(i,iter_in.t),size(mm.Z,1),cumsum(mm.erg_pz));
        end
        if iter_in.exog_deaths(i,iter_in.t-1) > 0
            % break down exogenous deaths that occur between t-1 and t down by b.o.p. z state:
            iter_in.die_cli_zst(i,:) = createDieVec(iter_in.lag_cli_zst(i,:).*iter_in.keep_cli,iter_in.exog_deaths(i,iter_in.t-1),size(mm.Z,1));
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
    iter_in.cur_cli_zst = iter_in.new_cli_zst + iter_in.trans_zst;

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
    age       = iter_in.t*ones(mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx),1) - flr;     % age in periods. age=0 all year for firms with no shipment previous year
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
        [mat_tran,ship_cur,age_vec] = match_sales(mkt,mm,iter_in.trans_count,age,iter_in.pt_ndx,iter_in.macro_state_f(iter_in.t));

        % [mat_tran,ship_cur,age_vec] =...
        %     match_sales(mm.scale_f,mm.eta,iter_in.trans_count(:,:,:),age,mm.X_f(macro_state_f(t)),t,mm.poisCDF_shipments,...
        %     mm.max_ships,size(mm.Z,1),mm.Z,mm.Phi(mm.pt_type(pt_ndx,1)));

        % mat_tran:  [initial state, exporter id, ending state, match revenue]
        % ship_cur:   match's number of shipments within the current period
        % age_vec:    firm age (# periods)
    end

    if iter_in.season == 1
        iter_in.N_match = size(mat_tran,1);
    end

    % Concatenate time index, season index and year index with match variables
    % and collect results for all seasons in year in seas_tran
    iter_in.seas_tran{1,iter_in.season} = [[iter_in.t,iter_in.season,iter_in.year].*ones(size(mat_tran,1),1),mat_tran,ship_cur,age_vec];
    % iter_in.seas_tran: [t, season, year, initial state, exporter id, ending state, match revenue,
    %             # shipments, firm age (# periods)]
    iter_in.seas_Zcut(iter_in.season)   = drop_Zcut;

    iter_in.N_match = iter_in.N_match + size(mat_tran,1); % cumulate match count within current year
end