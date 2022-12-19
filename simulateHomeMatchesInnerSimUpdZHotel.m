%%                 SimulateHomeMatchesInnerSimUpdZHotel  

 % This script maps Z states onto each active match for each type-pt_ndx
 % buyer in each period t.

 % iterH_in.trans_count(:,:,i):  For firm i, each row, corresponds to a
 % b.o.p. z state and #'s in columns indicate number of its matches starting 
 % from that state and transiting to a particular e.o.p z state. Col 1 contains 
 % counts of exiting matches, other columns correspond to z values. Row
 % 1 corresponds to new matches

 % iterH_in.trans_zst(i,:): For each firm (i), gives # surviving matches 
 % in each e.o.p. z state (cols). Columns correspond to z values
 
 % JT: Essentially the same structure as simulatForeignMatchesInnerSimUpdZHotel. 
 % We could easily consolidate these functions.

%%

function[iterH_in] = simulateHomeMatchesInnerSimUpdZHotel(iterH_in, mm, policy)

t = iterH_in.t;
pt_ndx = iterH_in.pt_ndx;

%fprintf('\rIn SimulateHomeMatchesInnerSimUpdZHotel, t =%4.0f, pt_ndx = %4.0f\n', [t, pt_ndx] )
% if t >= 469
%     'pause in simulateHomeMatchesInnerSimUpdZHotel'
% end

% break down by buyer types (z)

    for i=1:mm.sim_firm_num_by_prod_succ_type(pt_ndx)

        % distribute new clients between t-1 and t across initial z-states       
        if iterH_in.add_cli_cnt(i,t) > 0
           iterH_in.new_cli_zst(i,:) = new_vec_C(iterH_in.add_cli_cnt(i,t),size(mm.Z,1),cumsum(mm.erg_pz)); 
           
        % PATCH: detect and redo rare cases where random draws failed (not needed on UNIX systems)
             flag = find(sum(iterH_in.new_cli_zst(i,:),2) - iterH_in.add_cli_cnt(i,t) ~= 0);
             if flag == 1
             iterH_in.new_cli_zst(i,:) = new_vec_C(iterH_in.add_cli_cnt(i,t),size(mm.Z,1),cumsum(mm.erg_pz));
             end         
             
        else
           iterH_in.new_cli_zst(i,:) = zeros(1,size(mm.Z,1));
        end

        % Among matches that die endogenously, break down exogenous match 
        % deaths that occur between t-1 and t down by z state:
        if iterH_in.exog_deaths(i,t-1) > 0
            iterH_in.die_cli_zst(i,:) = createDieVec(iterH_in.lag_cli_zst(i,:).*iterH_in.keep_cli_lag,iterH_in.exog_deaths(i,t-1),size(mm.Z,1));
        end
        
         % Get rid of all clients when the firm slot turns over
         if iterH_in.new_firm(i,t)*(1-iterH_in.new_firm(i,t-1)) == 1
             iterH_in.die_cli_zst(i,:) = iterH_in.lag_cli_zst(i,:);   
         end
 
         %SLOW!
%         iterH_in.trans_count(2:size(mm.Z,1)+1,1,i)...
%         = max([(iterH_in.lag_cli_zst(i,:).*(1-iterH_in.keep_cli_lag))', iterH_in.die_cli_zst(i,:)']')';
%     
        temp = max([(iterH_in.lag_cli_zst(i,:).*(1-iterH_in.keep_cli_lag))', iterH_in.die_cli_zst(i,:)']')'; 
        for rndx = 1:size(mm.Z,1)
            iterH_in.trans_count(rndx+1,1,i) = temp(rndx,1);
        end
        
        % Update column 1 of trans_count(:,:,i) so that it contains counts of
        % all exiting matches (endog. and exog.), by buyer type (row).

        iterH_in.surviv_zst(i,:) = iterH_in.lag_cli_zst(i,:) - iterH_in.trans_count(2:size(mm.Z,1)+1,1,i)';
        % Update surviving client counts by z type using transition matrix for z.
        % Do this for those that don't die for endogenous or exogenous reasons.
        
        N_sur = sum(iterH_in.surviv_zst(i,:),2); 
        % Number of survivors from t-1 by b.o.p. z, firm i

        % Draw new z states for surviving matches; update trans_zst and trans_count
        if N_sur > 0
            sur_typ = find(iterH_in.surviv_zst(i,:)); % addresses for z-states populated by at least one survivor
            for jj = sur_typ  % loop over initial states of surviving matches, exporter i
                draw = rand(iterH_in.surviv_zst(i,jj),1);
                % identify destination z states for each surviving client (could be multiple survivors per initial type):
                trans_z = ones(iterH_in.surviv_zst(i,jj),1)*policy.pmat_cum_z(jj,:) > draw;
                % count # clients in each destination z state for each beginning z state.
                
      %SLOW!: % Record counts in cols 2:size(mm.Z,1)+1 of trans_count. Rows are initial states, plus 1:
      %         iterH_in.trans_count(jj+1,2:size(mm.Z,1)+1,i) = sum(trans_z(:,1:size(mm.Z,1)) - [zeros(size(draw,1),1),trans_z(:,1:size(mm.Z,1)-1)],1);
       
               temp = sum(trans_z(:,1:size(mm.Z,1)) - [zeros(size(draw,1),1),trans_z(:,1:size(mm.Z,1)-1)],1);
               for rndx=1:size(mm.Z,1)
                 iterH_in.trans_count(jj+1,rndx+1,i) = temp(rndx);
               end

                
      %SLOW!:  % cumulate over b.o.p. z types to get vector of surviving client e.o.p. types. Rows (i) are exporters:
%              iterH_in.trans_zst(i,:) = iterH_in.trans_zst(i,:) +  iterH_in.trans_count(jj+1,2:size(mm.Z,1)+1,i);

%               temp = iterH_in.trans_zst(i,:) +  iterH_in.trans_count(jj+1,2:size(mm.Z,1)+1,i);
               for cndx=1:size(mm.Z,1)
 %                 iterH_in.trans_zst(i,cndx) = temp(cndx);
                  iterH_in.trans_zst(i,cndx) = iterH_in.trans_zst(i,cndx) +  iterH_in.trans_count(jj+1,cndx+1,i);
               end
               
            end
        end
        
        % load new client z counts into first row of iterH_in.trans_count
        if sum(iterH_in.new_cli_zst(i,:),2)>0
            iterH_in.trans_count(1,2:size(mm.Z,1)+1,i) = iterH_in.new_cli_zst(i,:);
        end
    end
    
    iterH_in.cur_cli_zst = iterH_in.new_cli_zst + iterH_in.trans_zst;
    
