%%                 SimulateHomeMatchesInnerSimUpdZHotel  

% This script maps Z states onto each active match for each type-pt_ndx
% buyer in each period t.

function[iterH_in] = simulateHomeMatchesInnerSimUpdZHotel(iterH_in, mm, policy)

t = iterH_in.t;
pt_ndx = iterH_in.pt_ndx;

% break down by buyer types (z)

    for i=1:mm.sim_firm_num_by_prod_succ_type(pt_ndx)
        % break down new clients that occur between t-1 and t into z-types

        % distribute gross additions        
        if iterH_in.add_cli_cnt(i,t) > 0
           iterH_in.new_cli_zst(i,:) = new_vec_C(iterH_in.add_cli_cnt(i,t),size(mm.Z,1),cumsum(mm.erg_pz)); 
        else
           iterH_in.new_cli_zst(i,:) = zeros(1,size(mm.Z,1));
        end

        % break down exogenous match deaths that occur between t-1 and t down by z state, and
        % record number of endogenous plus exogenous exits by z state in 1st column of trans_count:
        if iterH_in.exog_deaths(i,t-1) > 0
            iterH_in.die_cli_zst(i,:) = createDieVec(iterH_in.lag_cli_zst(i,:).*iterH_in.keep.cli,iterH_in.exog_deaths(i,t-1),size(mm.Z,1));
        end
        
         % get rid of all clients when the firm slot turns over
         if iterH_in.new_firm(i,t)*(1-iterH_in.new_firm(i,t-1)) == 1
             iterH_in.die_cli_zst(i,:) = iterH_in.lag_cli_zst(i,:);
             %.*iterH_in.keep_cli;     
         end
        
      iterH_in.trans_count(2:size(mm.Z,1)+1,1,i) = max([(iterH_in.lag_cli_zst(i,:).*(1-iterH_in.keep.cli))', iterH_in.die_cli_zst(i,:)']')';
        % For each firm (i) of a particular type, column 1 of trans_mat now
        % contains counts of all exiting matches, by buyer type (row).

        % Update surviving client counts by z type using transition matrix for
        % z. Do this for those that don't die for endogenous reasons, minus those
        % that die for exogenous reasons:

        iterH_in.surviv_zst(i,:) = iterH_in.lag_cli_zst(i,:) - iterH_in.trans_count(2:size(mm.Z,1)+1,1,i)';
 %        iterH_in.surviv_zst(i,:) = (iterH_in.lag_cli_zst(i,:) - iterH_in.trans_count(2:size(mm.Z,1)+1,1,i)').*iterH_in.keep_cli;   
        
        N_sur = sum(iterH_in.surviv_zst(i,:),2); % number of survivors from t-1, firm i

        if N_sur > 0
            sur_typ = find(iterH_in.surviv_zst(i,:)); % addresses for z-states populated by at least one survivor
            for jj = sur_typ  % loop over initial states of surviving matches, exporter i
                draw = rand(iterH_in.surviv_zst(i,jj),1);
                % identify destination z states for each surviving client (could be multiple survivors per initial type):
                trans_z = ones(iterH_in.surviv_zst(i,jj),1)*policy.pmat_cum_z(jj,:) > draw;
                % count # clients in each destination z state for each beginning z state.
                % Record counts in cols 2:size(mm.Z,1)+1 of trans_mat. Rows are initial states, plus 1:
                iterH_in.trans_count(jj+1,2:size(mm.Z,1)+1,i) = sum(trans_z(:,1:size(mm.Z,1)) - [zeros(size(draw,1),1),trans_z(:,1:size(mm.Z,1)-1)],1);
                % cumulate over b.o.p. z types to get vector of surviving client e.o.p. types. Rows (i) are exporters:
                iterH_in.trans_zst(i,:) = iterH_in.trans_zst(i,:) +  iterH_in.trans_count(jj+1,2:size(mm.Z,1)+1,i);
            end
        end
        if sum(iterH_in.new_cli_zst(i,:),2)>0
            iterH_in.trans_count(1,2:size(mm.Z,1)+1,i) = iterH_in.new_cli_zst(i,:);
        end
    end
    
    iterH_in.cur_cli_zst = iterH_in.new_cli_zst + iterH_in.trans_zst;
    
