function [Q,st_vec] = makeq(L,D,n)
% Takes jump hazard (L), jump size (D), and desired number  
% of states (2n+1). we require this to be the same as in the 
% calculation of the maximum likelihood parameters. Returns a 2n+1 
% by 2n+1 intensity matrix Q  as well as a 2n+1 by 1 state vector 
% st_vec centered at zero.

Q = zeros(2*n+1);

Q(1,1) = -L;
Q(1,2) = L;
Q(2*n+1,2*n) = L;
Q(2*n+1,2*n+1) = -L;

for k=2:2*n
   Q(k,k) = -L;
   Q(k,k-1) = L*(.5*(1+((k-n-1)/(n))));
   Q(k,k+1) = L*(.5*(1-((k-n-1)/(n))));
end

st_vec = ((1:2*n+1)'-(n+1))*D;

end