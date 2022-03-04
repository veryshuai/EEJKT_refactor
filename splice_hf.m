function iter_out = splice_hf(sim_out_h,sim_out_f,policy,mm)

% This function takes panel of realizations on domestic and foreign sales
% for a particular type (theta_f and phi), and splices them together to 
% create synthetic firms

% sim_out_f.firm_f_yr_sales: [t,type,firm ID, total exports,total # foreign shipments,firm age in export mkt.]
% sim_out_h.firm_h_yr_sales: [t,type,firm ID, total dom. sales, total # dom. shipments,firm age in dom. mkt.]
% policy.firm_type_prod_succ_macro: [type, macro state index, theta index, productivity index]

% some firms with active export relationships have zero exports. Weed them out:
exporter = sim_out_f.firm_f_yr_sales(:,4)>0;
sim_out_f.firm_f_yr_sales = sim_out_f.firm_f_yr_sales(exporter,:); % stacks all firm-yr combs. with exports, given micro type
% same deal for domestic market shipments
dom_shipper = sim_out_h.firm_h_yr_sales(:,4)>0;  % stacks all firm-yr combs. with dom .sales, given micro type
sim_out_h.firm_h_yr_sales = sim_out_h.firm_h_yr_sales(dom_shipper,:);

iter_out.nexptr = sum(exporter>0);
iter_out.nhfirms = sum(dom_shipper>0);

type_h = sim_out_h.firm_h_yr_sales(:,2);
type_f = sim_out_f.firm_f_yr_sales(:,2);
prod_h  = policy.firm_type_prod_succ_macro(type_h,4);
prod_f  = policy.firm_type_prod_succ_macro(type_f,4);

try % confirm domestic and foreign records are for the same firm productivity type
    minobs = min(size(type_h,1),size(type_f,1));
    assert(sum((prod_h(1:minobs)- prod_f(1:minobs)).^2)==0);
catch
     fprintf('Warning: productivity different in different markets');
end

%% Extract observations on home and foreign sales with same firm_ID and date

% create unique identifiers for each period/firm_ID pair
obs_id_h = sim_out_h.firm_h_yr_sales(:,3) + (1/(mm.periods+1))*sim_out_h.firm_h_yr_sales(:,1); % micro type and period
obs_id_f = sim_out_f.firm_f_yr_sales(:,3) + (1/(mm.periods+1))*sim_out_f.firm_f_yr_sales(:,1); % micro type and period
% find home market obs. with identifiers that appear in export mkt. obs.
aa = ismember(obs_id_h,obs_id_f);
temp1 = obs_id_h(aa);
% find foreign market obs. with identifiers that appear in home mkt.
bb = ismember(obs_id_f,temp1);
try
  assert(sum(aa)==sum(bb))
catch
  fprintf('\r Warning: home-foreign merge discrepancy in spice_hf');
end
iter_out.hf_nobs  = sum(aa);
iter_out.nfirm = iter_out.nhfirms + iter_out.nexptr - iter_out.hf_nobs; % number of firms with sales in at least one market

% dom_only = ones(nobs_h,1) - aa > 0; % firms with dom. sales only
 ex_only  = ones(iter_out.nexptr,1) - bb > 0; % firms with exports only

theta_f = policy.firm_type_prod_succ_macro(type_f,3); % Each element of this vector is common to all firms
theta_h = sim_out_h.theta_h_firm;      % This is a vector of random draws--one per firm
if sum(bb) == 0
  temp3 = zeros(0,15);
else
  temp3 = [theta_f(bb),theta_h(aa),prod_h(aa),sim_out_h.firm_h_yr_sales(aa,:),sim_out_f.firm_f_yr_sales(bb,:)];
end
% sim_out_h.firm_h_yr_sales: [t,type,firm ID,total dom. sales, total # dom. shipments,firm age in dom. mkt.]
% sim_out_f.firm_f_yr_sales: [t,type,firm ID,total exports,total # foreign shipments,firm age in export mkt.]

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

% export_rate = [sales_hf(:,7)./(sales_hf(:,6)+sales_hf(:,7));ones(sum(ex_only,1),1)];
% the vector of ones in export_rate accounts for firms that exclusively

iter_out.expt_rate = sales_hf(:,7)./(sales_hf(:,6)+sales_hf(:,7));
% this alternative version excludes firms active solely in the export market

% serve foreign markets

%% moments for regression of log foreign sales on log domestic sales

   iter_out.hf_nobs = size(sales_hf,1);
   if iter_out.hf_nobs > 0
    iter_out.y_hf = log(sales_hf(:,7));
    iter_out.x_hf = [ones(iter_out.hf_nobs,1),log(sales_hf(:,6))];
    iter_out.hfmoms_xx = iter_out.x_hf'*iter_out.x_hf;
    iter_out.hfmoms_xy = iter_out.x_hf'*iter_out.y_hf;
    iter_out.hfysum = sum(iter_out.y_hf);    
   else
        iter_out.x_hf = zeros(0,2);
        iter_out.y_hf = zeros(0,1);
        iter_out.hfmoms_xx = zeros(2,2);
        iter_out.hfmoms_xy = zeros(2,1);
        iter_out.hfysum    = 0;
        iter_out.hf_nobs = 0;
    end

end

