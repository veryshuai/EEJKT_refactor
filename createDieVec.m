function [die_vec] = createDieVec(cli_vec,die_cnt,n_Z)
die_vec = zeros(1,n_Z);
for ii=1:1:die_cnt  % loop over number of deaths
 cdf_curr = cumsum(cli_vec'./sum(cli_vec))' ; % construct cum distr. of types
 die_typ = rand > cdf_curr; % randomly sample death from client distribution
 drop = [1,die_typ(:,1:n_Z-1)] - die_typ(:,1:n_Z); %identify type that died
 cli_vec = cli_vec - drop;  % adjust client vector for death (sampling w/out replacement)
 die_vec = die_vec + drop;  % increment the death count vector
end
end
