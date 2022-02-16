% This function distributes new clients across z-types according to 
% ergodic (cdf_Z) distribution 

function [dtype] = new_vec(n_cli,n_Z,cdf_Z)
    newcli_typ = rand(n_cli,1) > cdf_Z'; 
    dtype = sum([ones(n_cli,1),newcli_typ(:,1:n_Z-1)] - newcli_typ(:,1:n_Z),1);
end