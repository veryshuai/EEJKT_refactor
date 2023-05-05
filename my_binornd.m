function [res] = my_binord(N,p)

% Alternative function for generating matrix of binominal draws


[row_cnt, col_cnt] = size(N);
res = zeros(row_cnt, col_cnt);
for ii=1:row_cnt
   for jj=1:col_cnt
      res(ii, jj) = sum(rand(1,N(ii,jj))<p(ii,jj));
   end
end

end

