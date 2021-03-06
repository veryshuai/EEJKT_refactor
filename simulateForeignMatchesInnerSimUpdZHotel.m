function iter_in = simulateForeignMatchesInnerSimUpdZHotel(mm, iter_in, policy)
    for i=1:mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx)
        % break down new clients that occur between t-1 and t into e.o.p. z-types
               
        if iter_in.add_cli_cnt(i,iter_in.t)> 0
           iter_in.new_cli_zst(i,:) = new_vec_C(iter_in.add_cli_cnt(i,iter_in.t),size(mm.Z,1),cumsum(mm.erg_pz)); % distribute gross additions
        else
            iter_in.new_cli_zst(i,:)= zeros(1,size(mm.Z,1));
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
end