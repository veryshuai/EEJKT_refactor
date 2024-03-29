% Called from SimulateForeignMatchesInnerAnnualize

function [mat_cont_2yr,mat_yr_sales,mat_yr_sales_lag,year_lag] =...
    mat_yr_splice_v2(iterX_in,mm,year)

% This function splices the current year's records on matches for a given
% firm type with last year's records for the same firm type. Splicing is
% done by firm ID and by matching last year's eoy Z with this year's boy Z.
% Once the two years are spliced match age variables are
% created. Note that the count of one-year olds does not include duds
% that sent sample shipments but did not establish a successful match.

%  mat_yr_sales: [(1) firm ID, (2) match-specific sales, (3) shipments,   
%      (4) boy Z, (5) eoy Z, (6) match age in periods, (7) firm age in periods]

    mat_yr_sales     = iterX_in.mat_yr_sales;
    mat_yr_sales_lag = iterX_in.mat_yr_sales_lag;
    Zcut_eoy_lag     = iterX_in.Zcut_eoy_lag;   
    
    mat_yr_sales     = sortrows(mat_yr_sales,[1,4,6]);

%% find matches to splice, recognizing firm_ID slots that flip occupants

 % Find matches in current year that correspond to boy incumbents, and 
 % drop matches that correspond to post-flip periods. Firms that flip in 
 % period 1 of current year require extra attention
    boy_noflip = iterX_in.new_firm(floor(mat_yr_sales(:,1)),iterX_in.t-11)==0;
    incumb = logical((mat_yr_sales(:,1)==(floor(mat_yr_sales(:,1))).*(mat_yr_sales(:,4)>Zcut_eoy_lag)).*boy_noflip); % boy Z > Zcut_eoy_lag 
    tmp_tran = mat_yr_sales(incumb,:);
 
  % find matches that were active at the end of last year 
    boy_noflip_lag = iterX_in.new_firm(floor(mat_yr_sales_lag(:,1)),iterX_in.t-11)==0;
    contin = logical((mat_yr_sales_lag(:,5)>Zcut_eoy_lag).*boy_noflip_lag);
    tmp_tran_lag  = mat_yr_sales_lag(contin,:);

%% calculate match ages and deal with firm turnover
try
  % find matches active at end of last period (may be redundant)

    tmp_tran_lag = sortrows(tmp_tran_lag,[1,5,6]);  
    tmp_tran     = sortrows(tmp_tran,[1,4,6]);

    % firm is older this year (no flip) 
    cont_find = (tmp_tran_lag(:,7)>0).*(tmp_tran(:,7) - tmp_tran_lag(:,7)) >  0;  
    % cumulate match age
    tmp_tran(cont_find ,6) = tmp_tran_lag(cont_find ,6) + tmp_tran(cont_find,6);

catch
        'problem in mat_yr_splice_v2, lines 42-46';
        fileID4 = fopen('results/EEJKT_error_log.txt','a');
        fprintf(fileID4,'\r\n  ');
        fprintf(fileID4,'\r\n problem in mat_yr_splice_v2, lines 42-46');
        fprintf(fileID4,'\r\n period = %.2f, firm type = %.2f, market = %.2f', [iterX_in.t iterX_in.pt_ndx iterX_in.mkt]);
        fprintf(fileID4,'\r\n params = ');
        fprintf(fileID4,'\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',mm.param_vec(1:6));
        fprintf(fileID4,'\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',mm.param_vec(7:12));
        fprintf(fileID4, '\r\n  ');   
        fclose(fileID4);
end
 %%  
 try
    last_yr_exit = logical(ones(size(mat_yr_sales_lag,1),1)-contin);
    % match was there last year but not this year
    mat_lastyr_lag = mat_yr_sales_lag(last_yr_exit,:);

 if sum(contin,1)>0

 % load match age for continuing matches into mat_yr_sales
    mat_yr_sales(incumb,6) = tmp_tran(:,6);     
     
% Check consistency of match records, then spice:
try
    assert(sum((tmp_tran_lag(:,5) - tmp_tran(:,4)).^2)==0) 
    % match e.o.y. Z in t-1 sames as b.o.y. Z in t
    no_match = find(tmp_tran_lag(:,5) - tmp_tran(:,4)~=0);
    mismatch = [tmp_tran_lag(no_match,:),tmp_tran(no_match,:)];
catch
    'problem at line 63 of mat_yr_splice_v2'
    find_problem = tmp_tran_lag(:,5)-tmp_tran(:,4)~=0;
    prob_obs = [tmp_tran_lag(find_problem,:),tmp_tran(find_problem,:)];
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
     fileID4 = fopen('results/EEJKT_error_log.txt','a');
     fprintf(fileID4,'\r\n  ');
     fprintf(fileID4,'\r\n problem in last block of mat_yr_splice_v2');
     fprintf(fileID4,'\r\n period = %.2f, firm type = %.2f, market = %.2f', [iterX_in.t iterX_in.pt_ndx iterX_in.mkt]);
     fprintf(fileID4,'\r\n params = ');
     fprintf(fileID4,'\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',mm.param_vec(1:6));
     fprintf(fileID4,'\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',mm.param_vec(7:12));
     fprintf(fileID4, '\r\n  ');   
     fclose(fileID4);
 end

end
