function sim_out = splice_hf(sim_out,transF,transH,policy,mm,pt_ndx)
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% This function takes panel of realizations on domestic and foreign sales
% for a particular type (theta_f and phi), and splices them together to 
% create synthetic firms

% transF contains records for populated rows of the following exporter transition variables from iter_in: 
 %   transF{pt_inx,1:6}: [find_xcli, cur_cli_cnt, cum_succ, cumage, new_firm, cum_meets]
 
% transH contains records for populated rows of the following home mkt. transition variables from iterH_in: 
 %   transH{pt_inx,1:5}: [find_hcli, cur_cli_cnt, cum_succ, cumage, new_firm]

% sim_out.firm_f_yr_sales: [t,type,firm ID, total exports,total # foreign shipments,firm age in export mkt.]
% sim_out.firm_h_yr_sales: [t,type,firm ID, total dom. sales, total # dom. shipments,firm age in dom. mkt.]

% policy.firm_type_prod_succ_macro: [type, macro state index, theta index, productivity index]
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%% Get rid of obs. with 0 shipments, market by market

some_shpmts_h = find(sim_out.firm_h_yr_sales(:,5)>0);
some_shpmts_f = find(sim_out.firm_f_yr_sales(:,5)>0);

sim_out.firm_h_yr_sales = sim_out.firm_h_yr_sales(some_shpmts_h,:);
sim_out.firm_f_yr_sales = sim_out.firm_f_yr_sales(some_shpmts_f,:);

%% Count the number of distinct firms in the home and the foreign database

% temp_h = sortrows(sim_out.firm_h_yr_sales,[3 1]);
% temp_f = sortrows(sim_out.firm_f_yr_sales,[3 1]);
% sim_out.nfirm  = sum((temp_h(2:end,3)-temp_h(1:end-1,3)~=0)) + ...
%            sum((temp_h(2:end,6)-temp_h(1:end-1,6)~=mm.pd_per_yr).*(temp_h(2:end,3)-temp_h(1:end-1,3)==0)) ;
% sim_out.nexptr = sum((temp_f(2:end,3)-temp_f(1:end-1,3)~=0)) + ...
%            sum((temp_f(2:end,6)-temp_f(1:end-1,6)~=mm.pd_per_yr).*(temp_f(2:end,3)-temp_f(1:end-1,3)==0)) ;

% new_entryH = [zeros(rowsH,1),((transH{pt_ndx,5}(:,2:end) - transH{pt_ndx,5}(:,1:end-1))==-1)];
% firm_ndxH = cumsum(new_entryH,2);
% sim_out.nfirm = sum(firm_ndxH(:,end));

Nh_firm_yrs = sum(sum(transH{pt_ndx,5}(:,2:mm.pd_per_yr) - transH{pt_ndx,5}(:,1:mm.pd_per_yr-1)==-1),2);
Nf_firm_yrs = sum(sum(transF{pt_ndx,5}(:,2:mm.pd_per_yr) - transF{pt_ndx,5}(:,1:mm.pd_per_yr-1)==-1),2);
sim_out.nfirm = Nh_firm_yrs/(mm.periods/mm.pd_per_yr);  % home mkt. firms per year
sim_out.nexptr = Nf_firm_yrs/(mm.periods/mm.pd_per_yr); % foreign mkt. firms per year

if Nh_firm_yrs > 0
  export_rate = Nf_firm_yrs/Nh_firm_yrs;
  fprintf('pt_ndx =%2.0f, t=%3.0f\n',[pt_ndx, max(sim_out.firm_h_yr_sales(:,1))])  
  fprintf('%2.0f exporters and %3.0f firms, for export rate of %3.2f\n',[Nf_firm_yrs Nh_firm_yrs export_rate])  
end

rowsH = size(transH{pt_ndx,5},1);
rowsF = size(transF{pt_ndx,5},1);
if rowsH*rowsF>0

% firm_idH  = transH{pt_ndx,1}*ones(1,mm.periods) + 0.001*firm_ndxH;

% N_firm_yrs  = 0;
% N_Xfirm_yrs  = 0;
% yr_counter = 0;
% ub=12;
% while ub <= mm.periods - mm.pd_per_yr
%     lb=ub+1;
%     ub=lb+mm.pd_per_yr-1;
%     N_firm_yrs = N_firm_yrs + sum(sum(transH{pt_ndx,5}(:,lb:ub) - ...
%           (transH{pt_ndx,5}(:,lb-1:ub-1)==-1) )); 
%     N_Xfirm_yrs = N_Xfirm_yrs + length(unique(firm_idH(1:rowsF,lb:ub) ...
%             .*(transF{pt_ndx,2}(:,lb:ub)>0) ));
%     yr_counter = yr_counter  + 1;
% end
% sim_out.nfirm = N_firm_yrs/yr_counter;
% sim_out.nexptr = N_Xfirm_yrs/yr_counter;

%% Extract observations on home and foreign sales with same firm_ID and date

% create unique identifiers for each period/firm_ID pair
obs_id_h = sim_out.firm_h_yr_sales(:,3) + (1/mm.periods+1)*sim_out.firm_h_yr_sales(:,1); % micro type and period
obs_id_f = sim_out.firm_f_yr_sales(:,3) + (1/mm.periods+1)*sim_out.firm_f_yr_sales(:,1); % micro type and period

% find double occurances of firm ID/year pairs (due to firm switching)

type_h  = sim_out.firm_h_yr_sales(:,2);
type_f  = sim_out.firm_f_yr_sales(:,2);
theta_f = policy.firm_type_prod_succ_macro(type_f,3); % Each element of this vector is common to all firms
theta_h = sim_out.theta_h_firm(some_shpmts_h);      % This is a vector of random draws--one per firm
prod_h  = policy.firm_type_prod_succ_macro(type_h,4);
prod_f  = policy.firm_type_prod_succ_macro(type_f,4);

sim_out_h_dat = sortrows([obs_id_h,sim_out.firm_h_yr_sales,theta_h,prod_h],1);
sim_out_f_dat = sortrows([obs_id_f,sim_out.firm_f_yr_sales,theta_f,prod_f],1);
find_same_h   = find([1;sim_out_h_dat(2:end,1)-sim_out_h_dat(1:end-1,1)==0]);
find_same_f   = find([1;sim_out_f_dat(2:end,1)-sim_out_f_dat(1:end-1,1)==0]);

%% when duplicate firm_ID/year pairs occur, drop the one with larger firm age

keeper_h = ones(size(sim_out_h_dat(:,1),1),1);
for jj = find_same_h(2:end)
   drop = (sim_out_h_dat(jj,6)<=sim_out_h_dat(jj-1,6));
   keeper_h(jj) = drop;
   keeper_h(jj-1) = 1-drop;
end

keeper_f = ones(size(sim_out_f_dat(:,1),1),1);
for jj = find_same_f(2:end)
   drop = (sim_out_f_dat(jj,6)<=sim_out_f_dat(jj-1,6));
   keeper_f(jj) = drop;
   keeper_f(jj-1) = 1-drop;
end

 dup_h  = length(find_same_h);
 ntot_h = length(sim_out_h_dat(:,1)); 
%  if dup_h > 0 && ntot_h >0
%    fprintf('\rWarning: %2.0f of %5.0f records in firm_h_yr_sales have same firm_ID-period \n',[dup_h, ntot_h])  
%  end
%  dup_f  = length(find_same_f);
%  ntot_f = length(sim_out_f_dat(:,1)); 
%  if dup_f > 0 && ntot_f > 0
%    fprintf('Warning: %2.0f of %5.0f records in firm_f_yr_sales have same firm_ID-period \n',[dup_f, ntot_f])  
%  end

%% unbundle sim_out_f_dat, now that it's variables are compatibly sorted
obs_id_f                = sim_out_f_dat(logical(keeper_f),1);
sim_out.firm_f_yr_sales = sim_out_f_dat(logical(keeper_f),2:7);
theta_f                 = sim_out_f_dat(logical(keeper_f),8); 
prod_f                  = sim_out_f_dat(logical(keeper_f),9);

obs_id_h                = sim_out_h_dat(logical(keeper_h),1);
sim_out.firm_h_yr_sales = sim_out_h_dat(logical(keeper_h),2:7);
theta_h                 = sim_out_h_dat(logical(keeper_h),8);      
prod_h                  = sim_out_h_dat(logical(keeper_h),9);

assert(length(obs_id_h)==length(unique((obs_id_h))));
assert(length(obs_id_f)==length(unique((obs_id_f))));

%% Merge home and foreign sales records by firm_ID and period (t)

try % confirm domestic and foreign records are for the same firm productivity type
    minobs = min(size(type_h,1),size(type_f,1));
    assert(sum((prod_h(1:minobs)- prod_f(1:minobs)).^2)==0);
catch
     fprintf('Warning: productivity different in different markets \n');
end

aa = ismember(obs_id_h,obs_id_f);
temp1 = obs_id_h(aa);

% find foreign market obs. with identifiers that appear in home mkt.
bb = ismember(obs_id_f,temp1);
try
  assert(sum(aa)==sum(bb))
catch
  fprintf('\r Warning: home-foreign merge discrepancy in spice_hf');
end


if sum(bb) == 0
  temp3 = zeros(0,15);
%  sim_out.nexptr = 0;
else
    try
       
  temp3 = [theta_f(bb),theta_h(aa),prod_h(aa),sim_out.firm_h_yr_sales(aa,:),...
           sim_out.firm_f_yr_sales(bb,:)];
       
  % temp3: (1) theta_f, (2) theta_h, (3) prod_h , (4-9) firm_h_yr_sales, (10-15) firm_f_yr_sales]
  % where:
  % firm_f_yr_sales: [t,type,firm ID, annual for. sales, # expt. shpmts, firm age in export mkt.]
  % firm_h_yr_sales: [t,type,firm ID, annual dom. sales, # dom. shmts,   firm age in dom. mkt.]
 
    catch
        'problem in splice_hf line 60'
    end
    
both_mkts       = logical((temp3(:,7)>0).*(temp3(:,13)>0));
sim_out.nhfirms = sum(both_mkts>0);
sim_out.hf_nobs = sum(aa);

sim_out.nfirm = sim_out.nhfirms + sim_out.nexptr - sim_out.hf_nobs; % number of firms with sales in at least one market

end
% sim_out.firm_h_yr_sales: [t,type,firm ID,total dom. sales, total # dom. shipments,firm age in dom. mkt.]
% sim_out.firm_f_yr_sales: [t,type,firm ID,total exports,total # foreign shipments,firm age in export mkt.]

try % comfirm home and foreign records are for same period and firm ID
assert(sum((temp3(:,4) - temp3(:,10)).^2)==0); % periods match?
assert(sum((temp3(:,6) - temp3(:,12)).^2)==0); % firm ID's match?
sales_splice = [temp3(:,1:4),temp3(:,6),temp3(:,7),temp3(:,13),temp3(:,9),temp3(:,15)];

% sales splice contains records for firms with matched foreign and domestic sales: 
%  [1) theta_f, 2) theta_h, 3) prod_h, 4) t, 5) firm ID, 
%   6) total dom. sales, 7) total exports, 
%   8) age in dom. mkt, 9) age in export mkt.]
catch
   warning('time periods or firm IDs do not match for home versus foreign sales');
   sales_splice = zeros(0,9);
end

% exiting firms are replaced by new ones with the same ID, so need to
% exclude cases where the "firm" age in the export market exceeds the firm
% age in the domestic market to avoid inherited foreign client bases

same_firm   = sales_splice(:,8) >= sales_splice(:,9);
sales_hf    = sales_splice(same_firm,1:7);

% uncomment next 3 lines to view obs. used for regressions
% sales_hf_temp     = sales_splice(same_firm,:);
% eyeball           = sortrows(sales_hf_temp,[1 2 5 4 8]);
% eyeball_expt_rate = [eyeball, eyeball(:,7)./(eyeball(:,6)+eyeball(:,7))];

% export_rate = [sales_hf(:,7)./(sales_hf(:,6)+sales_hf(:,7));ones(sum(ex_only,1),1)];
% the vector of ones in export_rate accounts for firms that exclusively

sim_out.expt_rate = sales_hf(:,7)./(sales_hf(:,6)+sales_hf(:,7));

%% moments for regression of log foreign sales on log domestic sales

   sim_out.hf_nobs   = size(sales_hf,1);
   sim_out.y_hf      = log(sales_hf(:,7));
   sim_out.x_hf      = [ones(sim_out.hf_nobs,1),log(sales_hf(:,6))];
   sim_out.hfmoms_xx = sim_out.x_hf'*sim_out.x_hf;
   sim_out.hfmoms_xy = sim_out.x_hf'*sim_out.y_hf;
   sim_out.hfysum    = sum(sim_out.y_hf);    
else
   sim_out.hf_nobs   = 0;
   sim_out.x_hf      = zeros(0,2);
   sim_out.y_hf      = zeros(0,1);
   sim_out.hfmoms_xx = zeros(2,2);
   sim_out.hfmoms_xy = zeros(2,1);
   sim_out.hfysum    = 0;
   sim_out.hf_nobs   = 0; 
   sim_out.expt_rate = 0;
end


end

