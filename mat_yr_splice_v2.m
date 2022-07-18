% Called from SimulateForeignMatchesInnerAnnualize

function [mat_cont_2yr,mat_yr_sales,mat_yr_sales_adj,year_lag] =...
    mat_yr_splice_v2(mat_yr_sales,mat_yr_sales_lag,mm,year_lag,year)

% This function splices the current year's records on matches for a given
% firm type with last year's records for the same firm type. Splicing is
% done by firm ID and by matching last year's eoy Z with this year's boy Z.
% Once the two years are spliced match age variables are
% created. Note that the count of one-year olds does not include singletons
% that sent sample shipments but did not establish a successful match.

%  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in yrs, firm age in periods] 
   
    mat_yr_sales = sortrows(mat_yr_sales,[1,4,6]);

%% find matches to splice

    ff_cont      = mat_yr_sales(:,4)>0;        % current year, boy Z > 0
    ff_lag       = mat_yr_sales_lag(:,5)>0;    % lagged year, eoy Z > 0
    ff_lag_noco  = mat_yr_sales_lag(:,5)==0;   % lagged match not continuing: eoy Z = 0
    tmp_tran     = mat_yr_sales(ff_cont,:);    % continuing matches beginning of this year (boy Z > 0)
    tmp_tran_lag = mat_yr_sales_lag(ff_lag,:); % continuing matches end of last year (eoy Z > 0)
%     mat_cont_lag = sortrows(tmp_tran_lag,[1,5]); % sort lagged matches by firm ID and eoy Z
%     mat_noco_lag = mat_yr_sales_lag(ff_lag_noco,:); % lagged matches with eoy Z == 0
          

%% calculate match ages and deal with firm turnover
    
%  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, cum. match age, firm age] 
   
% if year ==5
%     'pause here'
% end

  % update match ages for continuing firms
    tmp_tran_lag = sortrows(tmp_tran_lag,[1,5,6]);
    tmp_tran     = sortrows(tmp_tran,[1,4,6]);
    cont_find = tmp_tran(:,7) - tmp_tran_lag(:,7) >  0;
    n_firms = length(cont_find);
    for ss=1:n_firms
        firm_ss = tmp_tran_lag(:,1) == ss;
        tmp_tran(firm_ss,6) = cont_find(ss).*tmp_tran_lag(firm_ss,6) + mm.pd_per_yr*ones(sum(firm_ss),1);
    end
        
 %%  
 

 if sum(ff_cont,1)>0

 % update match age for continuing matches
    mat_yr_sales(ff_cont,6) = tmp_tran(:,6);     
     
% Check consistency of match records, then spice:
    assert(sum((tmp_tran_lag(:,5) - tmp_tran(:,4)).^2)==0) % match e.o.y. Z in t-1 sames as b.o.y. Z in t
    mat_cont_2yr = [tmp_tran_lag, tmp_tran];
     
%  Stack non-continuing lagged matches with continuing lagged matches, modifying 
%  the latter so that eop Z = 0 if the match generates zero sales next period.
%     ff_ghost = logical((mat_yr_sales(ff_cont,3)==0).*(mat_cont_lag(:,5)>0)); 
%     mat_cont_lag(ff_ghost,5) = 0; % set eoy Z=0 for continuers that generate no further sales
%     mat_yr_sales_adj  = cat(1,mat_cont_lag,mat_noco_lag);
    
    mat_yr_sales_adj  = mat_yr_sales_lag;
      
else % if no matches for this firm type in the previous year
    mat_cont_2yr = zeros(0,2*size(mat_yr_sales,2));
%   mat_yr_sales(:,6) = 1;
%   mat_yr_sales(:,7) = mm.pd_per_yr/2;
%   mat_yr_sales_adj  = mat_yr_sales_lag;
    mat_yr_sales_adj  = double.empty(0,size(mat_yr_sales,2));
end

year_lag = year;

end