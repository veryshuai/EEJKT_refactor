[exit_r, exit_c] = find(iter_in.exit_firm==1);
[dcum_meet_r, dcum_meet_c] = find(iter_in.cum_meets(:,2:end)-iter_in.cum_meets(:,1:end-1)<0);

% check whether exits always followed by new firm
for i = 1:size(exit_r)
  for j = 1:size(exit_c)
    if iter_in.exit_firm(exit_r(i),exit_c(j))==1
    assert(iter_in.new_firm(exit_r(i),exit_c(j)+1)==1)
    end
  end
end

% reformat sales matrix: sort by firm, then time
NNfirms = size(iter_in.exit_firm,1);
SalesObs =size(iter_out.firm_f_yr_sales,1);
SalesMat = double.empty(0,6);
TT = SalesObs/NNfirms;
for j =1:NNfirms
  for i =1:TT
    nn = (i-1)*NNfirms + j;
    SalesMat = [SalesMat; iter_out.firm_f_yr_sales(nn,:)];
  end
end
    
    