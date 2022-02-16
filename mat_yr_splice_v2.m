% Called from matchdat_gen_f.m

function [mat_cont_2yr,mat_yr_sales,mat_yr_sales_adj,year_lag] =...
    mat_yr_splice_v2(mat_yr_sales,mat_yr_sales_lag,mm,year_lag,year)

% This function splices the current year's records on matches for a given
% firm type with last year's records for the same firm type. Splicing is
% done by firm ID and by matching last year's eoy Z with this year's boy Z.
% Once the two years are spliced match age and firm age variables are
% created. Note that the count of one-year olds does not include singletons
% that sent sample shipments but did not establish a successful match.

%  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in yrs, firm age] 

% if year == year_lag+1 % condition for updating year and forming mat_cont_2yr
    
    mat_yr_sales = sortrows(mat_yr_sales,[1,4]);
    mat_yr_sales(:,6) = mat_yr_sales(:,3)>0; % if shipments, initialize match age to one year; will be replaced with cum age for cont. matches below
%   mat_yr_sales(:,7) = 0; % initialize firm age to 0; will be replaced with cum firm age for cont. firms below
%% find matches to splice

    ff_cont      = mat_yr_sales(:,4)>0;        % current year, boy Z > 0
    ff_lag       = mat_yr_sales_lag(:,5)>0;    % lagged year, eoy Z > 0
    ff_lag_noco  = mat_yr_sales_lag(:,5)==0;   % lagged match not continuing: eoy Z = 0
    tmp_tran     = mat_yr_sales(ff_cont,:);    % continuing matches beginning of this year (boy Z > 0)
    tmp_tran_lag = mat_yr_sales_lag(ff_lag,:); % continuing matches end of last year (eoy Z > 0)
    mat_cont_lag = sortrows(tmp_tran_lag,[1,5]); % sort lagged matches by firm ID and eoy Z
    mat_noco_lag = mat_yr_sales_lag(ff_lag_noco,:); % lagged matches with eoy Z == 0
      
%% cumulate match age
      
 %     increment match age by one year if match was active at beginning of period (boy Z>0, or tmp_tran(:,4)>0) 
 %     and shipments occur (tmp_tran(:,3)>0). If shipments occur but match  wasn't active at beginning of period, 
 %     set age to 1. If no shipments occur during period, set age=0, regardless of whether match is active.
     
%   keep incrementing match age in years so long as shipments occur or it remains active: boy Z > 0 and eoy Z > 0 
    add_match_yr =  tmp_tran(:,3)>0; % match began the current year active and generated shipments
    tmp_tran(:,6) = add_match_yr.*(tmp_tran_lag(:,6) + ones(sum(ff_cont,1),1));
    mat_yr_sales(ff_cont,:)  = tmp_tran;
    


%% convert firm age to years
    
%  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, cum. match age, firm age] 

   % if other matches correspond to continuing firm, retreive lagged age
   ff_ndx = unique(mat_yr_sales(:,1));
   ff_ndx_lag = unique(mat_yr_sales_lag(:,1));
   finder = ff_ndx*ones(1,size(ff_ndx_lag,1)) - ones(size(ff_ndx,1),1)*ff_ndx_lag' == 0; 
   try
    if size(ff_ndx_lag,2)>0
    lag_age = zeros(size(ff_ndx_lag,1),1);
     for nn = 1:size(ff_ndx_lag,1)
       fndx_lag = mat_yr_sales_lag(:,1) == ff_ndx_lag(nn);
     % find age of continuing firms, which had eop Z > 0. Others get a 0
       lag_age(nn)= mean(mat_yr_sales_lag(fndx_lag,7))*max(mat_yr_sales_lag(fndx_lag,5)>0);
     end
    end
   trans_age = finder*lag_age;
   if size(ff_ndx,2)>0
     for nn = 1:size(ff_ndx,1)
       fndx = mat_yr_sales(:,1) == ff_ndx(nn);
       pos_ship = max(mat_yr_sales(fndx,3)>0);
       mat_yr_sales(fndx,7) = pos_ship*(1 + trans_age(nn));
     end
   end
   catch
       'pause here'
   end
 %%  
 if sum(ff_cont,1)>0
 %%      
    mat_cont_2yr = [mat_cont_lag, mat_yr_sales(ff_cont,:)];
    
%%  
%  Stack non-continuing lagged matches with continuing lagged matches, modifying 
%  the latter so that eop Z = 0 if the match generates zero sales next period.
    ff_ghost = logical((mat_yr_sales(ff_cont,3)==0).*(mat_cont_lag(:,5)>0)); 
    mat_cont_lag(ff_ghost,5) = 0; % set eoy Z=0 for continuers that generate no further sales
    mat_yr_sales_adj  = cat(1,mat_cont_lag,mat_noco_lag);
%%    
%   year_lag = year;
      
else % if no matches for this firm type in the previous year
    mat_cont_2yr = zeros(0,2*size(mat_yr_sales,2));
    mat_yr_sales(:,6) = 1;
    mat_yr_sales(:,7) = 1;
    mat_yr_sales_adj  = mat_yr_sales_lag;
end

year_lag = year;

end