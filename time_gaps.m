function [time_gap2,mktexit_data] = time_gaps(iter_in,mm)

          pd_per_yr = mm.pd_per_yr;
          t         = iter_in.t;

%         exit_firm = iter_in.exit_firm;
%  Reach back 3 years to find last meeting by current firm
          cum_meets   = iter_in.cum_meets(:,t-3*pd_per_yr:t);
          cum_succ    = iter_in.cum_succ(:,t-3*pd_per_yr:t);
          cur_cli_cnt = iter_in.cur_cli_cnt(:,t-3*pd_per_yr:t);
         
          % identify periods where firm ID switches to new entrant
          new_meets  = cum_meets(:,2:end) - cum_meets(:,1:end-1);
          firm_flip  = new_meets(:,1:end)<0; 
          
%% create variables for firm exit regression
          yr_lag       = 2*pd_per_yr; % end of year period, last year          
          exit         = logical((min(new_meets(:,yr_lag+1:end),[],2) < 0).*(cur_cli_cnt(:,yr_lag)>0));  % flags firm_IDs that flip ownership during year
          no_exit      = logical((min(new_meets(:,yr_lag+1:end),[],2) == 0).*(cur_cli_cnt(:,yr_lag)>0)); % flags firm IDs that don't flip during previous yr

          cum_meets_xt = cum_meets(exit,yr_lag);   % cum meetings of exiting firms entering current yr
          cum_succ_xt  = cum_succ(exit,yr_lag);    % cum successes of exiting firms entering current yr
          cum_meets_ct = cum_meets(no_exit,yr_lag);% cum meetings of non-exiting entering current yr
          cum_succ_ct  = cum_succ(no_exit,yr_lag); % cum successes of non-exiting entering current yr
          Nexit        = size(cum_meets_xt,1);
          Ncont        = size( cum_meets_ct,1);
          
          exit_regdat  = [ones(Nexit,1),cum_meets_xt,cum_succ_xt];
          cont_regdat  = [zeros(Ncont,1),cum_meets_ct,cum_succ_ct];
          mktexit_data = [exit_regdat;cont_regdat];
       

%%   Create variables for match hazard rate regression
           Nr        = size(new_meets,1);
           [rr,cc]   = find(new_meets>0);  % new meetings over past 3 years
  
           if isempty(rr)==0             
           cell_add      = find(new_meets>0);  % firm-months with new meetings (vector of row addresses, stacked column by column)          
           cum_meet_int = cum_meets(:,2:end); % cum meetings, submatrix for 3 yr interval 
           cum_succ_int = cum_succ(:,2:end);  % cum successes, submatrix for 3 yr interval
  
           try
               
          % sort meeting by firm ID, carrying along relevant variables 
          if Nr>1
            temp = sortrows([rr cc new_meets(cell_add) cum_meet_int(cell_add) cum_succ_int(cell_add) firm_flip(rr,:)],1);
          elseif Nr==1 % patch dealing with annoying Matlab transposing convention
            rr = rr';
            cc = cc';
            cell_add = cell_add';
            temp = sortrows([rr cc new_meets(cell_add)' cum_meet_int(cell_add)' cum_succ_int(cell_add)' firm_flip(rr,:)],1);
          end
          rr            = temp(:,1);
          cc            = temp(:,2);
          new_meets_add = temp(:,3);
          cum_meet_add  = temp(:,4);
          cum_succ_add  = temp(:,5);
          firm_flip     = temp(:,6:end);

          % measure gap length and identify gaps that span a flipping period
           nnn = length(rr);
           gap = zeros(nnn-1,1);
%          flip_in_gap = false(nnn,1);
           same_firm = false(nnn,1);
           for j=2:nnn
               gap(j,1) = cc(j) - cc(j-1);
               same_firm(j) = ...
               logical(rr(j)-rr(j-1)==0 && sum(firm_flip(j,cc(j-1):cc(j)))==0) ;
           end
           
           %  temp2 contains observations on all firm-period pairs in which new meetings take place. 
           temp2 = [rr, cc, gap, new_meets_add, same_firm, cum_meet_add, cum_succ_add];
           % (1) firm_ID, (2) period w/in interval, (3) gap (4) # new meetings,(5)same firm
           % (6) cum. meetings, (7) cum succeses
          
          catch
          fprintf('\r\n problem in time gaps line 44-65, firm type  = %.3f\n',iter_in.pt_ndx);
           end
          time_gap = temp2(same_firm,:); 

          
          % only keep matches that happen in the current year
          time_gap = time_gap(time_gap(:,2)>2*pd_per_yr,:)  ;        
           % (1) firm_ID, (2) period w/in interval, (3) gap (4) # new meetings,(5)same firm
           % (6) cum. meetings, (7) cum succeses   
           
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
           % (1) firm_ID, (2) period w/in interval, (3) gap 
           % (4) # new meetings,(5) t (6) cum. meetings, (7) cum succeses 
           
           else
               time_gap2 = double.empty(0,7);
           end
          
end
