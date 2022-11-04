%%              SimulateForeignMatchesInnerSimUpdZHotel

 % This script maps Z states onto each active match for a given firm type
 % (pt_ndx) and period (t).

 % iter_in.trans_count(:,:,i):  For firm # i, each row, corresponds to a
 % b.o.p. z state and figures in columns indicate number of its matches starting 
 % from that state and transiting to a particular e.o.p z state. Col 1 contains 
 % counts of exiting matches, other columns correspond to z values. Row
 % 1 corresponds to new matches

 % iter_in.trans_zst(i,:): For each firm (i), gives # surviving matches 
 % in each e.o.p. z state (cols). Columns correspond to z values

%% 

function iter_in = simulateForeignMatchesInnerSimUpdZHotel(mm, iter_in, policy)

    for i=1:mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx)
        % break down new clients that occur between t-1 and t into e.o.p. z-types
               
        if iter_in.add_cli_cnt(i,iter_in.t)> 0
         % distribute gross additions across entry states
           iter_in.new_cli_zst(i,:) = new_vec_C(iter_in.add_cli_cnt(i,iter_in.t),size(mm.Z,1),cumsum(mm.erg_pz)); 

          % detect and redo rare cases where random draws failed (not needed on UNIX systems)
           flag = find(sum(iter_in.new_cli_zst(i,:),2) - iter_in.add_cli_cnt(i,iter_in.t) ~= 0);
             if flag == 1
             iter_in.new_cli_zst(i,:) = new_vec_C(iter_in.add_cli_cnt(i,iter_in.t),size(mm.Z,1),cumsum(mm.erg_pz));
             end          
                           
        else
            iter_in.new_cli_zst(i,:)= zeros(1,size(mm.Z,1));
        end        
                       
       % break down exogenous deaths that occur between t-1 and t down by b.o.p. z state:
        if iter_in.exog_deaths(i,iter_in.t-1) > 0
            iter_in.die_cli_zst(i,:) = createDieVec(iter_in.lag_cli_zst(i,:).*iter_in.keep_cli,iter_in.exog_deaths(i,iter_in.t-1),size(mm.Z,1));
        end
        
        % get rid of all clients when firm dies      
        if iter_in.new_firm(i,iter_in.t)*(1-iter_in.new_firm(i,iter_in.t-1)) == 1 
             iter_in.die_cli_zst(i,:) = iter_in.lag_cli_zst(i,:);
        end
                 
        iter_in.trans_count(2:size(mm.Z,1)+1,1,i) =  max([(iter_in.lag_cli_zst(i,:).*(1- iter_in.keep_cli))',  iter_in.die_cli_zst(i,:)']')';            
        % Update column 1 of trans_count(:,:,i) so that it contains counts of
        % all exiting matches (endog. and exog.), by buyer type (row).

        iter_in.surviv_zst(i,:) =  iter_in.lag_cli_zst(i,:) -  iter_in.trans_count(2:size(mm.Z,1)+1,1,i)';       
        % Update surviving client counts by z type using transition matrix for z.
        % Do this for those that don't die for endogenous or exogenous reasons.
        
        N_sur = sum(iter_in.surviv_zst(i,:),2); 
        % number of survivors from t-1 by b.o.p. z, firm i
        
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