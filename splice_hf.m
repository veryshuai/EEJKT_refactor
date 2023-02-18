function sim_out_pt_ndx = splice_hf(sim_out_pt_ndx,policy,mm,pt_ndx)
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% This function takes panel of realizations on domestic and foreign sales
% for a particular type (theta_f and phi), and splices them together to 
% create synthetic firms. It also counts firms that service the export
% market only, the home market only, or both, and constructs export to 
% total sales ratios for each exporter.

% transF contains records for populated rows of the following exporter transition variables from iter_in: 
 %   transF{pt_inx,1:6}: [find_xcli, cur_cli_cnt, cum_succ, cumage, new_firm, cum_meets]
 
% transH contains records for populated rows of the following home mkt. transition variables from iterH_in: 
 %   transH{pt_inx,1:5}: [find_hcli, cur_cli_cnt, cum_succ, cumage, new_firm]

% sim_out.firm_f_yr_sales: [t,type,firm ID, total exports,total # foreign shipments,firm age in export mkt.]
% sim_out.firm_h_yr_sales: [t,type,firm ID, total dom. sales, total # dom. shipments,firm age in dom. mkt.]
% sim_out.dud matches:     [t,type,firm_ID, shipment=1, sales, bop Z = eop Z = match_age = 0, firm age in export mkt.]

% policy.firm_type_macro_succ_prod: [type, macro state index, theta index, productivity index]
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%% Initialize counts

  sim_out_pt_ndx.nfirm   = 0;   % count, all active firms 
  sim_out_pt_ndx.nhfirms = 0;   % count, home mkt. firms 
  sim_out_pt_ndx.nexptr  = 0;   % count, foreign mkt. firms  
  sim_out_pt_ndx.hf_nobs = 0;   % count, firm-yrs in both mkts.   
 
  sim_out_pt_ndx.expt_rate = double.empty(0,1); 

  %% Sum duds within a year for each firm_ID after sorting by firm and period

  dud_matches = [sim_out_pt_ndx.dud_matches(:,1:3),sim_out_pt_ndx.dud_matches(:,6),...
                 sim_out_pt_ndx.dud_matches(:,4),sim_out_pt_ndx.dud_matches(:,9)];
  dud_firmyr        = dud_matches(:,1) + 0.000001*dud_matches(:,3);
  dud_temp_matches  = sortrows([dud_firmyr,dud_matches],1);
% dud_matches:                  [t, type, firm_ID, sales, shipments (1), firm age in export mkt.]
% dud_temp_matches: [dud_firmyr, t, type, firm_ID, sales, shipments (1), firm age in export mkt.]

  f_firmyr     = sim_out_pt_ndx.firm_f_yr_sales(:,1) + 0.000001*sim_out_pt_ndx.firm_f_yr_sales(:,3);
  f_temp_firmyr = sortrows([f_firmyr,sim_out_pt_ndx.firm_f_yr_sales],1);
% f_temp_firmyr contains non-dud export sales, aggregated to the annual level: 
%     [f_firmyr, t, type, firm_ID, sales, shipments (1), firm age in export mkt.]

% Aggregate dud sales to the annual level
dud_yr_sales = double.empty(0,7);
if ~isempty(dud_temp_matches)
    
    dud_yr_sales(1,:) = dud_temp_matches(1,:);
    for i = 2:size(dud_temp_matches,1)
        if dud_temp_matches(i,1)-dud_temp_matches(i-1,1)==0
            dud_yr_sales(end,5:6) = dud_yr_sales(end,5:6) + dud_temp_matches(i,5:6);
        else
           dud_yr_sales = [dud_yr_sales; dud_temp_matches(i,:)];
        end
    end
%dud_yr_sales: [dud_firmyr, t, type, firm_ID, total dud sales, shipment count, firm age in export mkt.]

%%  Consolidate annual dud sales with annual non-dud export sales, firm-year by firm-year
   
if size(unique(f_firmyr),1) - size(f_firmyr,1) ~= 0
  fprintf('Warning: multiple records for same firm-yr in firm_f_yr_sales  \n');
end

% Add dud records to successes when both are present for the same firm-yr.
 hh = ismember(dud_yr_sales(:,1),f_firmyr);  
 ff = ismember(f_firmyr,dud_yr_sales(:,1));  
 assert(sum((f_temp_firmyr(ff,1)-dud_yr_sales(hh,1)).^2)==0) 
 f_temp_firmyr(ff,5:6) = f_temp_firmyr(ff,5:6) + dud_yr_sales(hh,5:6);
% Stack remaining dud records onto the end of firm_f_yr_sales
 firm_f_yr_salesD = [f_temp_firmyr(ff,:); dud_yr_sales(logical(1-hh),:)]; 
% The D indicates duds are included

else
 firm_f_yr_salesD = f_temp_firmyr;
end
% firm_f_yr_salesD: [f_firmyr, t, type, firm_ID, all sales (incl. duds), 
%                    all shipments (incl. duds), firm age in export mkt.]


%% Splice export sales records with domestic market sales records

h_firmyr = sim_out_pt_ndx.firm_h_yr_sales(:,1) ...
           + 0.000001*sim_out_pt_ndx.firm_h_yr_sales(:,3);
firm_h_yr_sales = sortrows([h_firmyr,sim_out_pt_ndx.firm_h_yr_sales],1);
% firm_f_yr_salesD: [h_firmyr, t, type, firm_ID, all sales (incl. duds), 
%                    all shipments (incl. duds), firm age in export mkt.]

some_shpmts_h = find(firm_h_yr_sales(:,6)>0);
some_shpmts_f = find(firm_f_yr_salesD(:,6)>0);

firm_h_yr_sales  = firm_h_yr_sales(some_shpmts_h,:);
firm_f_yr_salesD = firm_f_yr_salesD(some_shpmts_f,:);
Nh_firm_yrs      = size(firm_h_yr_sales,1);
Nf_firm_yrs      = size(firm_f_yr_salesD,1);


if size(unique(h_firmyr),1) - size(h_firmyr,1) ~= 0
   fprintf('Warning: multiple records for same firm-yr in firm_h_yr_sales  \n');
end

 ff = ismember(firm_f_yr_salesD(:,1),firm_h_yr_sales(:,1));
 hh = ismember(firm_h_yr_sales(:,1),firm_f_yr_salesD(:,1));

 f_bothmkt = sortrows(firm_f_yr_salesD(ff,:),1);
 h_bothmkt = sortrows(firm_h_yr_sales(hh,:),1);

 type_f  = f_bothmkt(:,3);
 theta_f = mm.pt_type(type_f,2); 
 prod_f  = mm.pt_type(type_f,1); 

 type_h  = h_bothmkt(:,3);
 theta_h = policy.firm_type_macro_succ_prod(type_h,3); 
 prod_h  = policy.firm_type_macro_succ_prod(type_h,4);

try
 assert(sum((f_bothmkt(:,1)-h_bothmkt(:,1)).^2)==0); % same firm-period?
 assert(sum((prod_f-prod_h).^2)==0);                 % same productivity?
catch
  fprintf('\r Warning: year, firm ID or productivity mismatch for home versus foreign sales');
end

%% 

if sum(ff)*sum(hh) == 0
   both_mkt    = double.empty(0,15);
   export_rate = double.empty(0,1);
   Nhf_firm_yrs = 0;
else
    
  try
% Combine export and domestic sales records by firm type and year in both_mkt      
  both_mkt = [theta_f, theta_h ,prod_f, f_bothmkt(:,2:end),h_bothmkt(:,2:end)];       
  % both_mkt: [(1) theta_f, (2) theta_h, (3) prod_h , (4-9) firm_h_yr_sales(:,2:end),
  %            (10-15) firm_f_yr_sales(:,2:end)]
  catch
   fprintf('\r Warning: dimension compatibility problem in splice_hf line 122');
  end
  
%% Count the number of distinct firm-yrs in the home and the foreign database

% number of firms in both markets, aggregating over firms and years
  Nhf_firm_yrs = size(both_mkt,1);
  
 % number of firms with sales in at least one market
  N_firm_yrs  = Nh_firm_yrs + Nf_firm_yrs - Nhf_firm_yrs;  
 % number of pure exportes 
  Nexprt_only = Nf_firm_yrs - Nhf_firm_yrs;
 
  sim_out_pt_ndx.nfirm   = N_firm_yrs;   % load count, all active firms 
  sim_out_pt_ndx.nhfirms = Nh_firm_yrs;  % load count, home mkt. firms 
  sim_out_pt_ndx.nexptr  = Nf_firm_yrs;  % load count, foreign mkt. firms  
  sim_out_pt_ndx.hf_nobs = Nhf_firm_yrs; % load count, firms in both mkts.   
  
  % load vector of export shares among active exporters, including pure exporters
  sim_out_pt_ndx.expt_rate = [both_mkt(:,7)./(both_mkt(:,7)+both_mkt(:,14)); ones(Nexprt_only,1)];

  % sim_out_pt_ndx.hffirms = Nhf_firm_yrs; 
  sim_out_pt_ndx.nexptr  = Nf_firm_yrs;
  sim_out_pt_ndx.nhfirms = Nh_firm_yrs;

end
  sim_out_pt_ndx.firm_h_yr_sales = firm_h_yr_sales;
% sim_out_pt_ndx.firm_h_yr_sales: [t,type,firm ID,total dom. sales, total # dom. shipments,firm age in dom. mkt.]
  sim_out_pt_ndx.firm_f_yr_sales = firm_f_yr_salesD;
% sim_out_pt_ndx.firm_f_yr_sales: [t,type,firm ID,total exports,total # foreign shipments,firm age in export mkt.]

% uncomment next 3 lines to view obs. used for regressions
% sales_hf_temp     = sales_splice(same_firm,:);
% eyeball           = sortrows(sales_hf_temp,[1 2 5 4 8]);
% eyeball_expt_rate = [eyeball, eyeball(:,7)./(eyeball(:,6)+eyeball(:,7))];

%% Load matrices into sim_out_pdx
  sim_out_pt_ndx.firm_f_yr_sales = firm_f_yr_salesD(:,2:end);
  sim_out_pt_ndx.firm_h_yr_sales = firm_h_yr_sales(:,2:end);


%% moments for regression of log foreign sales on log domestic sales

if size(both_mkt,1)>0
   sim_out_pt_ndx.y_hf      = log(both_mkt(:,7));
   sim_out_pt_ndx.x_hf      = [ones(sim_out_pt_ndx.hf_nobs,1),log(both_mkt(:,13))];
   sim_out_pt_ndx.hfmoms_xx = sim_out_pt_ndx.x_hf'*sim_out_pt_ndx.x_hf;
   sim_out_pt_ndx.hfmoms_xy = sim_out_pt_ndx.x_hf'*sim_out_pt_ndx.y_hf;
   sim_out_pt_ndx.hfysum    = sum(sim_out_pt_ndx.y_hf);    
else
   sim_out_pt_ndx.hf_nobs   = 0;
   sim_out_pt_ndx.x_hf      = zeros(0,2);
   sim_out_pt_ndx.y_hf      = zeros(0,1);
   sim_out_pt_ndx.hfmoms_xx = zeros(2,2);
   sim_out_pt_ndx.hfmoms_xy = zeros(2,1);
   sim_out_pt_ndx.hfysum    = 0;
   sim_out_pt_ndx.hf_nobs   = 0; 
   sim_out_pt_ndx.expt_rate = 0;
end


end

