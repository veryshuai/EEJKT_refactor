% Called from SimulateForeignMatchesInnerAnnualize

function [mat_cont_2yr,mat_yr_sales,mat_yr_sales_lag,year_lag] =...
    mat_yr_splice_v2(mat_yr_sales,mat_yr_sales_lag,mm,year)

% if year==16 || year==17
%     
%     'pause here'
% 
% end

% This function splices the current year's records on matches for a given
% firm type with last year's records for the same firm type. Splicing is
% done by firm ID and by matching last year's eoy Z with this year's boy Z.
% Once the two years are spliced match age variables are
% created. Note that the count of one-year olds does not include singletons
% that sent sample shipments but did not establish a successful match.

%  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in yrs, firm age in periods] 
   
    mat_yr_sales = sortrows(mat_yr_sales,[1,4,6]);

%% find matches to splice, recognizing firm_ID slots that flip occupants

% Find matches in current year that correspond to boy incumbents, and 
% drop matches that correspond to post-flip periods.

%   active_ID = sortrows(unique(floor(mat_yr_sales(:,1))));
    ID_occur = sortrows(unique(mat_yr_sales(:,1)));
    flipper  = zeros(length(ID_occur),1);
    include  = ones(size(mat_yr_sales,1),1);
    for j = 1:length(ID_occur)
        flipper(j) = ID_occur(j) ~= floor(ID_occur(j)); % not boy incumbents
        if flipper(j)>0
        include = include - (mat_yr_sales(:,1)==(floor(mat_yr_sales(:,1))+0.5));
        end
    end
            
 % Find last year's matches corresponding to eoy active firms and drop matches
 % that correspond to pre-flip periods
 
%   active_ID_lag = sortrows(unique(floor(mat_yr_sales_lag(:,1))));
    ID_occur_lag  = sortrows((unique(mat_yr_sales_lag(:,1))));
    flipper_lag   = zeros(length(ID_occur_lag),1);
    include_lag   = ones(size(mat_yr_sales_lag,1),1);
    for j = 1:length(ID_occur_lag)
        flipper_lag(j) = ID_occur_lag(j)~=floor(ID_occur_lag(j));
        if flipper_lag(j)>0
        include_lag = include_lag - (mat_yr_sales_lag(:,1)==floor(ID_occur_lag(j)));
        end
    end  
    
    % further select on whether last year's matches were active eoy and
    % this year's matches were active boy.
    
    ff_cont      = logical(include.*(mat_yr_sales(:,4)>0));         % current year, boy Z > 0
    ff_lag       = logical(include_lag.*(mat_yr_sales_lag(:,5)>0)); % lagged year and firm age, eoy Z > 0
 try
    tmp_tran     = mat_yr_sales(ff_cont,:);    % continuing matches beginning of this year (boy Z > 0)
    tmp_tran_lag = mat_yr_sales_lag(ff_lag,:); % continuing matches end of last year (eoy Z > 0)
 catch
        'problem in mat_yr_splice_v2 line 53-54'
  end            

%% calculate match ages and deal with firm turnover
    
%  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, cum. match age, firm age] 
   
% if year == 5
%     'pause here'
% end

  % update match ages for continuing firms
    tmp_tran_lag = sortrows(tmp_tran_lag,[1,5,6]);
    tmp_tran     = sortrows(tmp_tran,[1,4,6]);
    
    
    'now in mat_yr_splice_v2'
    year
try
    cont_find = tmp_tran(:,7) - tmp_tran_lag(:,7) >  0;
catch
    'problem in mat_yr_splice_v2 72'
end
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
try
    assert(sum((tmp_tran_lag(:,5) - tmp_tran(:,4)).^2)==0) % match e.o.y. Z in t-1 sames as b.o.y. Z in t
catch
    'problem in line 80 of mat_yr_splice_v2'
end
    mat_cont_2yr = [tmp_tran_lag, tmp_tran];
     
%  Stack continuing lagged matches with non-continuing lagged matches
    
    mat_yr_sales_lag  = [tmp_tran_lag;mat_lastyr_lag];
      
else % if no matches for this firm type in the previous year
    mat_cont_2yr = zeros(0,2*size(mat_yr_sales,2));
    mat_yr_sales_lag  = double.empty(0,size(mat_yr_sales,2));
end

year_lag = year;

end