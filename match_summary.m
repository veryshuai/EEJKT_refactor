function [exit_by_age,brooks] =  match_summary(match_recs,mm);


% This function generates failure rates by initial sales and match age.
% More precisely, by type group and match age, where types are classified
% according to their initial sales.
 
% agg_mat_yr_sales includes all generated matches, so it can be used to analyze first-year
% matches that don't continue to the next year, and it can be used to analyze cohorts at the
% firm level.

% match_recs: [t,type,firm ID, match sales, shipments, boy Z, eoy Z, match age, firm age]    

% If shipments occur but a match wasn't active at beginning of period, we set its age to 1. 
% If no shipments occur during period, we set age=0, regardless of whether a match is active. 
% (see mat_yr_splice, line 17.) Match age increments by one per year if it remains active  
% and shipments occur.

% Firm age is 1 plus lagged firm age for continuing firms, given that shipments occur for 
% at least one match. When no shipments occur, firm age resets to 0.

%%  Preliminary data preparation

N_theta2 = size(mm.theta2',1);    % number of success rates (thetas)
N_Phi    = size(mm.Phi,1);        % number of exporter productivities
pt_type = [kron((1:N_Phi)',ones(N_theta2,1)),kron(ones(N_Phi,1),(1:N_theta2)')];
% prod_ndx  = pt_type(pt_ndx,1);
% theta_ndx = pt_type(pt_ndx,2);


 %   match_recs: [t,type,firm ID, match sales, shipments, boy Z, eoy Z, match age, firm age] 
     ff_ship = match_recs(:,5) > 0;
     match_recs = match_recs(ff_ship,:);
     match_recs(:,1) = floor(match_recs(:,1)/mm.pd_per_yr); % create integer year variable
       
     all_matches = match_recs(match_recs(:,1)>5,:); % throw out 5-yr burn-in period
   % all_matches = agg_mat_yr_sales(agg_mat_yr_sales(:,4)>0,:); % throw out matches with no sales
   % all_matches: [year,type,firm ID, sales, shipments, boy Z, eoy Z, match age, firm age]   
     
     all_matches = all_matches(all_matches(:,9)>0,:);  % exclude the few firms with age 0
     yr1 = all_matches(:,8)==1; % pick matches in first year 
  
     new_matches = all_matches(yr1,:); % matches in their 1st yr.
   % new_matches: [year, type, firm_ID, sales, shipments, boy Z, eoy Z, match age w/in yr, firm age] 
      
%% Sort types into 4 groups by initial sales using the new_matches data   
  
     type_yr = new_matches(:,2) + 0.001*new_matches(:,1); % create type-year id, new matches
     typ_yr_list = unique(type_yr);
     active_typ_yr = size(typ_yr_list,1);
     xx_new = zeros(active_typ_yr,1); % total sales by type and year
     nx_new = zeros(active_typ_yr,1); % number of new matches by type and year
     ag_new = zeros(active_typ_yr,1); % age
     tp_new = zeros(active_typ_yr,1); % type
     yr_new = zeros(active_typ_yr,1); % calendar year
     th_new = zeros(active_typ_yr,1); % theta index 
     
   % means and sums by type-yr 
     for jj=1:active_typ_yr
        ff = type_yr == ones(size(type_yr,1),1).*typ_yr_list(jj); % pick off type-year of interest
        xx_new(jj) = sum(new_matches(ff,4),1);       % total type-year specific sales
        nx_new(jj) = sum(ff);                        % number of firms in type-year category 
        ag_new(jj) = mean(floor(new_matches(ff,9))); % avg. age of firms making sales
        tp_new(jj) = mean(new_matches(ff,2));        % type index  
        th_new(jj) = mm.theta2(pt_type(tp_new(jj),2)); % theta index 
        yr_new(jj) = mean(new_matches(ff,1));        % calendar year 
     end
       
     % find mean initial (new match) sales by type
     new_sales = [tp_new,yr_new,nx_new,xx_new];
     type_list     = unique(new_sales(:,1)); %list of active types, pooling years
     active_typ    = size(type_list,1);
     th_typ        = zeros(active_typ,1);
     new_typ_sales = zeros(active_typ,1);
     new_type      = zeros(active_typ,1);
     new_cnt       = zeros(active_typ,1);
     dud_cnt       = zeros(active_typ,1);
     dud_sales     = zeros(active_typ,1);
     dud_size      = ones(0,1);
     avg_sales     = mean(new_sales(:,4));
     for jj=1:active_typ % get mean sales and counts, by type
        ff = new_sales(:,1) == ones(size(new_sales,1),1).*type_list(jj);
        th_typ(jj)        = mean(th_new(ff,1),1);
        new_type(jj)      = mean(new_sales(ff,1));
        new_cnt(jj)       = sum(new_sales(ff,3),1);
        dud_cnt(jj)       = sum(new_sales(ff,3).*(1-th_new(ff))./th_new(ff),1);
        new_typ_sales(jj) = mean(new_sales(ff,4));
        dud_sales(jj)     = 0.1*mean(new_sales(ff,4));
        dud_sales(jj)     = 0.5*avg_sales;
        dud_size          = cat(1,dud_size,dud_sales(jj)*ones(sum(ff,1),1));
     end
       
%%    group types according to mean initial sales quartile    
%    new_sales2 = sort([new_type,new_cnt,new_typ_sales],3); % sort types by mean sales of new matches
     adj_sales = th_typ.*new_typ_sales+(1-th_typ).*dud_sales;
     % adjust initial sales to include duds, making figures compatible with
     % C.J.'s Brooks table figures
     new_sales2 = sortrows([new_type,new_cnt+dud_cnt,adj_sales,dud_sales],3);
     cum_cnt = cumsum(new_sales2(:,2))./sum(new_sales2(:,2));
     ndx_q = cell(4,1);
     typ_q = cell(4,1);
     mean_q = zeros(4,1);
     match_cnt_q = zeros(4,1);
     upb = 0.25;
     lowb = 0;
     for ii = 1:4
        ndx_q{ii}  = cum_cnt>lowb & cum_cnt<=upb; % vector of row identifiers for types in sales quartile ii
        lowb = upb;
        upb = upb + 0.25;
        % quartile-wide avg. sales:
        mean_q(ii) = sum(new_sales2(ndx_q{ii},3).*new_sales2(ndx_q{ii},2))/sum(new_sales2(ndx_q{ii},2));
        match_cnt_q(ii) = sum(new_sales2(ndx_q{ii},2),1);
        typ_q{ii}  = new_sales2(ndx_q{ii},1); % vector of types in quartile i      
     end
 
%% Reorganize columns in all_matches (holdover from earlier version of code)
 
 % all_matches: [year, type, firm ID, sales, shipments, boy Z, eoy Z, match age, firm age] 
   succ_rate   = mm.theta2(pt_type(all_matches(:,2),2))';

 % for each new match, (1-succ_rate)/succ_rate duds. Successful matches plus duds: 
   effec_match = yr1.*(1-succ_rate)./succ_rate + 1;
 % recall yr1 is a dummy for matches with age = 1 
   
 % effective matches scales matches in their first year to impute duds.
 % Again the idea is to be compatible with C.J.'s Brooks tables.
 
   splicedat  = [all_matches(:,2),all_matches(:,4),all_matches(:,6:9),...
                 all_matches(:,3),all_matches(:,1), effec_match];    
  % splicedat: [type, sales, boy Z, eoy Z, match age, firm age integer, firm_ID, year, effective matches]

 % Create data sets for each initial size quartile; include new and
 % continuing matches.
     matdat_q = cell(4,1);
     for qq = 1:4
     aggdat_q = zeros(0,9); 
         for jj=1:size(typ_q{qq},1)     
         ndx_q    = splicedat(:,1)==typ_q{qq}(jj);
         aggdat_q = [aggdat_q;splicedat(ndx_q,:)];
         end
         matdat_q{qq} = aggdat_q;
     end
%% match survival and age 

   max_mat_age = max(splicedat(:,5));   
   exit_rate = zeros(max_mat_age,4);
      for qq=1:4
      ff_age = unique(matdat_q{qq}(matdat_q{qq}(:,5)>0,5)); % ignore matches yielding zero sales
      matdat_qs = matdat_q{qq}(:,4:5); % eoy Z and age, all matches in quantile qq
      matdat_qm = matdat_q{qq}(:,9);   % effective matches, all matches in quantile qq
      sales_qq  = matdat_q{qq}(:,2);   % sales, all  matches in quantile qq
      theta_qq  = mm.theta2(pt_type(matdat_q{qq}(:,1),2))';
      age_types = size(ff_age,1);
      for aa = 1:age_types
        gg = find(matdat_qs(:,2) == ff_age(aa));  
        if ff_age(aa)==1
%        exit_rate(aa,qq) = (sum(matdat_qs(gg,1)==0) + sum(1-theta_qq(gg)))/sum(matdat_qm(gg,1),1);
        exit_rate(aa,qq) = (sum(matdat_qs(gg,1)==0) + sum((1-theta_qq(gg)).*matdat_qm(gg,:)))/...
            sum(matdat_qm(gg,:),1);
        else
       exit_rate(aa,qq) = sum(matdat_qs(gg,1)==0)/sum(matdat_qm(gg,:),1);
        end
 %       exit_rate(aa,qq) = 1- sum(matdat_qs(gg,1)>0)/sum(matdat_qm(gg,1),1);
        
       end
     end
     exit_by_age = exit_rate(1:5,:);
   
     % NOTE: the age 0 firms have no shipments, although they are active.

%% Brooks tables 

%all_matches: [year, type, firm_ID, sales, shipments, boy Z, eoy Z, match age, firm age] 

% Pooling all matches through time, find total sales (xx) for each type/firm-ID/year/firm_age
% combination. Need this stage to avoid averaging over different firms with
% the same type, firm-ID, and age.
    firm_ndx = all_matches(:,2) + 1000*all_matches(:,3) + 0.01*all_matches(:,1) + 100000*floor(all_matches(:,9));
  % firm_ndx1 = type + 1000*firm_ID + 0.01*year + 100000*firm_age
    
    big_id_list = unique(firm_ndx);
    ntypes = size(big_id_list,1);
    xx = zeros(ntypes,1);
    nx = zeros(ntypes,1);
    ag = zeros(ntypes,1);
    tp = zeros(ntypes,1);
    id =zeros(ntypes,1);
    for jj=1:ntypes
        ff = firm_ndx == ones(size(firm_ndx,1),1).*big_id_list(jj);
        xx(jj) = sum(all_matches(ff,4),1); % sales
        nx(jj) = sum(ff);                  % # obs for firm_ndx1
        ag(jj) = ceil(mean(floor(all_matches(ff,9)))); % firm age, current yr
        tp(jj) = mean(all_matches(ff,2));  % firm type
        id(jj) = mean(all_matches(ff,3));  % firm ID
    end
    
 % now ignore type and firm_ID; just distinguish by firm age 
    age_list = unique(ag);
    n_age = size(age_list,1);
    nx_age = zeros(n_age,1); 
    ag_age = zeros(n_age,1);
    mx_age = zeros(n_age,1);
    xx_age = zeros(n_age,1);
    for jj=1:n_age
        ff = ag == ones(size(ag,1),1).*age_list(jj);
        nx_age(jj) = sum(ff); 
        ag_age(jj) = mean(ag(ff));
        mx_age(jj) = mean(xx(ff));
        xx_age(jj) = sum(xx(ff));
    end
    
    brooks = [nx_age,xx_age,mx_age];


end