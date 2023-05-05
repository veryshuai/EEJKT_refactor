function [varF] = cum_cov(CDF_in,N)

% quantiles: cutoffs for categories
% CDF_in: if a row, mass below each cutoff; 
%          if a matrix, each row corresponds to an initial state, and gives 
%          transition probs. conditioned on that state

ncuts = size(CDF_in,2); 
nrows = size(CDF_in,1);

varF = zeros(ncuts,ncuts,nrows);

for k = 1:nrows
    for j = 1:ncuts
        for i=1:j
            varF(i,j,k) = CDF_in(k,i)*(1-CDF_in(k,j))/N(k);
            varF(j,i,k) = varF(i,j,k);
        end
    end
end

end