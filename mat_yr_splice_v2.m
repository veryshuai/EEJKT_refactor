% Called from SimulateForeignMatchesInnerAnnualize

function [mat_cont_2yr,mat_yr_sales,mat_yr_sales_lag,year_lag] =...
    mat_yr_splice_v2(mat_yr_sales,mat_yr_sales_lag,mm,year)

% This function splices the current year's records on matches for a given
% firm type with last year's records for the same firm type. Splicing is
% done by firm ID and by matching last year's eoy Z with this year's boy Z.
% Once the two years are spliced match age variables are
% created. Note that the count of one-year olds does not include singletons
% that sent sample shipments but did not establish a successful match.

%  mat_yr_sales: [(1) firm ID, (2) match-specific sales, (3) shipments,   
%      (4) boy Z, (5) eoy Z, (6) match age in yrs, (7) firm age in periods]
   
    mat_yr_sales = sortrows(mat_yr_sales,[1,4,6]);

%% find matches to splice, recognizing firm_ID slots that flip occupants

% Find matches in current year that correspond to boy incumbents, and 
% drop matches that correspond to post-flip periods.
   
    incumb = mat_yr_sales(:,1)==(floor(mat_yr_sales(:,1))).*(mat_yr_sales(:,4)>0);
    contin = mat_yr_sales_lag(:,5)>0;
    tmp_tran = mat_yr_sales(incumb,:);
    tmp_tran_lag = mat_yr_sales_lag(contin,:);
   

%% calculate match ages and deal with firm turnover
try
  % update match ages for continuing firms
    floorID      = floor(tmp_tran_lag(:,1));
    tmp_tran_lag = sortrows([tmp_tran_lag,floorID],[8,5,6]);
    tmp_tran_lag =  tmp_tran_lag(:,1:7);
    tmp_tran     = sortrows(tmp_tran,[1,4,6]);
    
 fprintf('\rNow at mat_yr_splice_v2, line 37. Evaluating year %2.0f\n', year)   

    cont_find = tmp_tran(:,7) - tmp_tran_lag(:,7) >  0;


    n_firms = length(cont_find);
    for ss=1:n_firms
        firm_ss = tmp_tran_lag(:,1) == ss;
        tmp_tran(firm_ss,6) = cont_find(ss).*tmp_tran_lag(firm_ss,6) + mm.pd_per_yr*ones(sum(firm_ss),1);
    end
catch
        'problem in mat_yr_splice_v2, lines 42-46'
end
 %%  
 try
%  mat_lastyr_lag = mat_yr_sales_lag(logical(include_lag-ff_lag),:);
last_yr_exit = logical(ones(size(mat_yr_sales_lag,1),1)-contin);
 mat_lastyr_lag = mat_yr_sales_lag(last_yr_exit,:);

 if sum(contin,1)>0

 % update match age for continuing matches
    mat_yr_sales(contin,6) = tmp_tran(:,6);     
     
% Check consistency of match records, then spice:
try
    assert(sum((tmp_tran_lag(:,5) - tmp_tran(:,4)).^2)==0) % match e.o.y. Z in t-1 sames as b.o.y. Z in t
catch
    'problem at line 63 of mat_yr_splice_v2'
    find_prob = tmp_tran_lag(:,5)-tmp_tran(:,4)~=0;
    prob_obs = [tmp_tran_lag(find_prob,:),tmp_tran(find_prob,:)];
end
    mat_cont_2yr = [tmp_tran_lag, tmp_tran];
     
%  Stack continuing lagged matches with non-continuing lagged matches
    
    mat_yr_sales_lag  = [tmp_tran_lag;mat_lastyr_lag];
      
else % if no matches for this firm type in the previous year
    mat_cont_2yr = zeros(0,2*size(mat_yr_sales,2));
    mat_yr_sales_lag  = double.empty(0,size(mat_yr_sales,2));
end

year_lag = year;

 catch
     'problem in last block of mat_yr_splice_v2'
 end

end
