function [time_gap2,mktexit_data] = time_gaps(iter_in,mm)

          pd_per_yr = mm.pd_per_yr;
          t         = iter_in.t;
          exit_firm = iter_in.exit_firm;
          cum_meets = iter_in.cum_meets;
          cum_succ  = iter_in.cum_succ;
   
          yr_lag3   = t-3*pd_per_yr+1; % reach back as far as 3 years looking for last match before first current year match  
          firm_flip = exit_firm(:,yr_lag3:t); % isn't needed but need to be careful removing it from matrices
          new_meet  = cum_meets(:,yr_lag3:t) - cum_meets(:,yr_lag3-1:t-1); % deal with firm ID changes later
          
%% create variables for firm exit regression
          
          exit         = exit_firm(:,t-pd_per_yr+1:t) == 1;  % flags firm_IDs that flip ownership during year
          no_exit      = sum(exit_firm(:,t-pd_per_yr+1:t)==0,2) == pd_per_yr; % flags firm IDs that don't flip during year
          cum_meets_xt = cum_meets(:,t-pd_per_yr:t-1).*exit; % cum meetings of exiting firms entering current yr
          cum_succ_xt  = cum_succ(:,t-pd_per_yr:t-1).*exit;  % cum successes of exiting firms entering current yr
          cum_meets_ct = cum_meets(:,t-pd_per_yr).*no_exit;  % cum meetings of incumbents entering current yr
          cum_succ_ct  = cum_succ(:,t-pd_per_yr).*no_exit;   % cum successes of incumbents entering current yr
          ff_exit      = find(cum_succ_xt > 0);
          ff_cont      = find(cum_succ_ct > 0);  
          Nexit        = size(ff_exit,1)*size(ff_exit,2);
          Ncont        = size(ff_cont,1)*size(ff_cont,2);
          
          exit_regdat  = [ones(Nexit,1),reshape(cum_meets_xt(ff_exit),Nexit,1),reshape(cum_succ_xt(ff_exit),Nexit,1)];
          cont_regdat  = [zeros(Ncont ,1),reshape(cum_meets_ct(ff_cont),Ncont,1),reshape(cum_succ_ct(ff_cont),Ncont,1)];
          mktexit_data = [exit_regdat;cont_regdat];

          try
              assert(sum(mktexit_data(:,2)-mktexit_data(:,3)<0)==0)
          catch
              warning('successes exceed meetings in transformed data')
              mktexit_data
          end         

%%   Create variables for match hazard rate regression

           [rr,cc]   = find(new_meet>0);          % new meetings over past 3 years
           if size(rr,1)*size(rr,2)>0
               
           cell_add      = find(new_meet>0);      % firm-months with new meetings (vector of row addresses, stacked column by column)          
           cum_meet_int = cum_meets(:,yr_lag3:t); % cum meetings, submatrix for 3 yr interval 
           cum_succ_int = cum_succ(:,yr_lag3:t);  % cum successes, submatrix for 3 yr interval
                             
           temp = sortrows([rr cc new_meet(cell_add) firm_flip(cell_add) cum_meet_int(cell_add) cum_succ_int(cell_add)],1);
           %  temp contains observations on all firm-period pairs in which new meetings take place. 
           % (1) firm_ID, (2) period w/in interval, (3)# new meetings, (4) firm replacement, (5) cum. meetings, (6) cum succeses 
           nnn = size(temp,1);
           
           exit_flag    = exit_firm(sortrows(rr),yr_lag3:t);
           same_firm_ID = temp(2:nnn,1)-temp(1:nnn-1,1)==0;
                            
          tdiff = [temp(2:nnn,1:2), (temp(2:nnn,2)-temp(1:nnn-1,2)), temp(2:nnn,3:4), temp(1:nnn-1,5:6)]...
                   .*same_firm_ID;  % zero out diffs across different firms
          % (1) firm_ID, (2) period w/in interval, (3) time gap between meetings, 
          % (4) # new meetings (5) firm replacement dummy, (6) cum. meetings as of previous meeting, 
          % (7) cum succeseses as of previous meeting

%%        Drop rows that correspond to same firm ID but different firms      

          t_start = temp(1:nnn-1,2);
          t_end   = temp(2:nnn,2);
          
          [rr2,~]   = find(tdiff(:,1)>0); 
          tdiff2    = tdiff(rr2,:);
          t_gap_lag = t_start(rr2);
          t_span    = [t_gap_lag, max([t_end(rr2) t_gap_lag],[],2)];
          
          same_firm = zeros(length(rr2),1);
          for jj = 1:length(rr2)
            same_firm(jj) = sum(exit_flag(jj,t_span(jj,1):t_span(jj,2)),2)==0;
          end
          
          time_gap = tdiff2(same_firm==1,:);
          % (1) firm (2) time w/in interval (3) time gap (4) # new meetings, 
          % (5) firm replacement, (6) lagged cum. meetings, (7) lagged cum succeses

          % only keep matches that happen in the current year
          time_gap = time_gap(time_gap(:,2)>2*pd_per_yr,:)  ;        

           
%%         Deal with cases of multiple new meetings within a single period

           % Each of the k matches is treated as spread evenly
           % within the period. Hence the time gap before the first one goes
           % back to the previous period in which a match occurred, and the time gaps
           % for the others are intra-period.
           
           multi_meet = time_gap(:,4)>1;
           extra_meet0 = time_gap(multi_meet,:);
           extra_meet0(:,4) = extra_meet0(:,4)-1; % first instance of within-period match already accounted for
           extra_meet1 = zeros(sum(extra_meet0(:,4)),size(extra_meet0,2));
           cntr = 1;
           for ii=1:size(extra_meet0,1)
               extra_meet1(cntr:cntr+extra_meet0(ii,4)-1,:) = ones(extra_meet0(ii,4),1).*extra_meet0(ii,:);
               extra_meet1(cntr:cntr+extra_meet0(ii,4)-1,3) = 1./(1+extra_meet1(cntr:cntr+extra_meet0(ii,4)-1,4));
               cntr = cntr+extra_meet0(ii,4);
           end
           extra_meet1(:,4) = ones(size(extra_meet1,1),1); % now each row is a single meeting
  
           % stack the multi-meeting observations with the others:
           time_gap2 = [time_gap; extra_meet1];
           time_gap2 = sortrows(time_gap2, [1 2]); % could do without the sort--it's just for eye-balling
           time_gap2(:,5) = time_gap2(:,2)+ (t-3*pd_per_yr).*ones(size(time_gap2,1),1); % replace firm replacement indicator with time  

%            if size(time_gap,1) > 20
%                'pause here'
%            end
           
           else
               time_gap2 = double.empty(0,7);
           end
          
end
