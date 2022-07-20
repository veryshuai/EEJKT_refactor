% Called from SimulateForeignMatchesInnerAnnualize

function [mat_cont_2yr,mat_yr_sales,mat_yr_sales_lag,year_lag] =...
    mat_yr_splice_v2(mat_yr_sales,mat_yr_sales_lag,mm,year)

% This function splices the current year's records on matches for a given
% firm type with last year's records for the same firm type. Splicing is
% done by firm ID and by matching last year's eoy Z with this year's boy Z.
% Once the two years are spliced match age variables are
% created. Note that the count of one-year olds does not include singletons
% that sent sample shipments but did not establish a successful match.

%  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in yrs, firm age in periods] 
   
    mat_yr_sales = sortrows(mat_yr_sales,[1,4,6]);

%% find matches to splice, recognizing firm_ID slots that flip occupants

% find matches in current year that correspond to boy incumbents
    active_ID = sortrows(unique(floor(mat_yr_sales(:,1))));
    ID_occur = sort(floor(unique(mat_yr_sales(:,1))));
    flipper = zeros(length(active_ID),1);
    include = ones(size(mat_yr_sales,1),1);
    for j = 1:length(active_ID)
        flipper(j) = sum(ID_occur==j)-1;
        if flipper(j)>0
        include = include - (mat_yr_sales(:,1)==(active_ID(j)+0.5));
        end
    end
            
 % find last year's matches corresponding to eoy active firms    
    active_ID = sortrows(unique(floor(mat_yr_sales_lag(:,1))));
    ID_occur = sort(floor(unique(mat_yr_sales_lag(:,1))));
    flipper_lag = zeros(length(active_ID),1);
    include_lag = ones(size(mat_yr_sales_lag,1),1);
    for j = 1:length(active_ID)
        flipper_lag(j) = sum(ID_occur==j)-1;
        if flipper_lag(j)>0
        include_lag = include_lag - (mat_yr_sales_lag(:,1)==(active_ID(j)));
        end
    end  
    
    % further select on whether last year's matches were active eoy and
    % this year's matches were active boy.
    
    ff_cont      = include.*mat_yr_sales(:,4)>0;         % current year, boy Z > 0
    ff_lag       = include_lag.*mat_yr_sales_lag(:,5)>0; % lagged year, eoy Z > 0
    tmp_tran     = mat_yr_sales(ff_cont,:);    % continuing matches beginning of this year (boy Z > 0)
    tmp_tran_lag = mat_yr_sales_lag(ff_lag,:); % continuing matches end of last year (eoy Z > 0)
             

%% calculate match ages and deal with firm turnover
    
%  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, cum. match age, firm age] 
   
if year == 5
    'pause here'
end

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
 mat_lastyr_lag = mat_yr_sales_lag(logical(include_lag-ff_lag),:);

 if sum(ff_cont,1)>0

 % update match age for continuing matches
    mat_yr_sales(ff_cont,6) = tmp_tran(:,6);     
     
% Check consistency of match records, then spice:
    assert(sum((tmp_tran_lag(:,5) - tmp_tran(:,4)).^2)==0) % match e.o.y. Z in t-1 sames as b.o.y. Z in t
    mat_cont_2yr = [tmp_tran_lag, tmp_tran];
     
%  Stack continuing lagged matches with non-continuing lagged matches
    
    mat_yr_sales_lag  = [tmp_tran_lag;mat_lastyr_lag];
      
else % if no matches for this firm type in the previous year
    mat_cont_2yr = zeros(0,2*size(mat_yr_sales,2));
    mat_yr_sales_lag  = double.empty(0,size(mat_yr_sales,2));
end

year_lag = year;

end