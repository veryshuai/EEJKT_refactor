function [pmat_type,pmat_cum,Q_size] = Q2pmat_v2(mm,nn2,typemat,lambda,Q_size,Q_index,f_dummy)

% Called from inten_sim_v1, this function uses intensity matricies to 
% calculates the discrete time period-to-period state transition probabilities 
% implied by the policy functions.

%% construct intensity matrix for each firm type 

% lambda_f (succ, trial, common succ rate (defunct), network size, prod of firm, F macro shock) 
% lambda_h (common succ rate (defunct), known succ rate, network size, prod of firm, H macro shock)

qmat_type = cell(size(typemat,1),1);
pmat_type = cell(size(typemat,1),1);
pmat_cum  = cell(size(typemat,1),1);

%   mo_increm       = mm.mo_increm;  % months per unit time interval
  frac_of_year    = 1/mm.pd_per_yr;  % Portion of year to which transition probabilites will apply
  firm_death_haz  = mm.d;            % per-period (not per year) death hazard
  
%% Create a matrix for taking out rounding errors in matrix exponentiation  
%  Don't want impossible events to occur with non-zero probability.
  nonzero_tran = zeros(Q_size,Q_size); 
    for i = 1:Q_size
        for j = 1:Q_size
               cond1 = (Q_index(j,2) - Q_index(i,2)<0)*(Q_index(j,1)>1);
               cond2 = (Q_index(j,3) - Q_index(i,3)<0)*(Q_index(j,1)>1);
               imposs_tran = max(cond1,cond2);             
               % when cum matches or cum successes fall, they must go all the way to zero
            if f_dummy == 1
               cond3 = ((Q_index(j,3) - Q_index(i,3)) > (Q_index(j,2) - Q_index(i,2)))*(Q_index(j,1)>1);
               % # new successes cannot be greater than # new meetings for continuing firms
               imposs_tran = max(imposs_tran,cond3);
            end
            nonzero_tran(i,j) = 1 - imposs_tran;
        end
    end
%% Now ready to loop over types 
for typ_indx = 1:size(typemat,1)  
    %set type
    ms             = typemat(typ_indx,2);
    succ_prob      = mm.theta2(typemat(typ_indx,3));
    prod_lvl       = typemat(typ_indx,4);
    
    counter = 0;
  
 if f_dummy==1  % for foreign market transitions
    q_index_list = zeros(nn2*(nn2+1)/2,3);    
    for i=1:1:nn2  % number of meetings, plus 1
        pos_i = sum(1:i-1); %advance to the correct number of meetings
        for ss=1:1:i % number of successes, plus 1
             r_ind = pos_i+ss;
             search_inten = lambda(ss,i,1,ss,prod_lvl,ms);
             
             if i ~= nn2 && i ~= 1 % deal separately with first and last nn1.  It is impossible to learn more from last trial    
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,r_ind,-(search_inten + firm_death_haz)]; %diag
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,pos_i+i+ss,(1-succ_prob)*search_inten]; %trial and failure
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,pos_i+i+ss+1,succ_prob*search_inten]; %trial and success
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,1,firm_death_haz]; %firm death
            elseif i == nn2
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,r_ind,-(firm_death_haz)]; %no learning (or searching for now)
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,1,firm_death_haz]; %firm death hazard
            elseif i == 1
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,r_ind,-(search_inten)]; %no firm death hazard, no clients to lose
                 counter = counter + 1;
                q_index_list(counter,:) = [r_ind,pos_i+i+ss,(1-succ_prob)*search_inten]; %trial and failure
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,pos_i+i+ss+1,succ_prob*search_inten]; %trial and success
            end
        end
    end   
    qmat_type{typ_indx} = sparse(q_index_list(:,1),q_index_list(:,2),q_index_list(:,3));
 
 else % for domestic market transitions (f_dummy==2)
      
        q_index_list = zeros(nn2,3);
        for ss=1:1:Q_size % number of successes, plus 1
             if ss <= nn2
                  search_inten = lambda(1,typemat(typ_indx,3),ss,prod_lvl,ms);
             else
                  search_inten = lambda(1,typemat(typ_indx,3),nn2,prod_lvl,ms);
             end
             
% deal separately with first and last nn1. It is impossible to learn more from last trial
             if ss ~= Q_size && ss ~= 1 
                counter = counter + 1;
                q_index_list(counter,:) = [ss,ss,-(succ_prob*search_inten + firm_death_haz)]; %diag
                counter = counter + 1;
                q_index_list(counter,:) = [ss,ss+1,succ_prob*search_inten]; %trial and success
                counter = counter + 1;
                q_index_list(counter,:) = [ss,1,firm_death_haz]; %firm death
            elseif ss == Q_size
                counter = counter + 1;
                q_index_list(counter,:) = [ss,ss,-(firm_death_haz)]; 
                counter = counter + 1;
                q_index_list(counter,:) = [ss,1,firm_death_haz]; 
            elseif ss == 1
                counter = counter + 1;
                q_index_list(counter,:) = [1,1,-succ_prob*search_inten]; %no firm death hazard, no success
                 counter = counter + 1;
                q_index_list(counter,:) = [1,2,succ_prob*search_inten]; %trial and success
            end
        end
       qmat_type{typ_indx} = sparse(q_index_list(:,1),q_index_list(:,2),q_index_list(:,3));
 end   
%% construct transition probabilities for discrete time intervals

%  pmat_type{typ_indx} = nonzero_tran.*expm(frac_of_year.*qmat_type{typ_indx}); 
   pmat_type{typ_indx} = nonzero_tran.*expm(qmat_type{typ_indx}); 
   pmat_type{typ_indx} = pmat_type{typ_indx}./sum(pmat_type{typ_indx},2);
 % Division by sum() is needed because of rounding error in expm()
   pmat_cum{typ_indx}  = cumsum(pmat_type{typ_indx}')';
end

end




