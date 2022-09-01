function erg = make_erg(L,D,X)
%This function takes the arrival rate, jump size, and shock vector
%and returns a size by 1 vector of ergodic unconditional probabilities over 
%the states.  The results about the way the discrete process converges to 
%an OU process are from Shimer (2005).

%last updated 2011-8-10 by David Jinkins

n = size(X,1);

gamma = L/n;
sig = L^.5*D;

var = sig^2/(2*gamma);

erg = normpdf(X,0,var^.5); %get normal pdf value at each shock value
erg = erg./sum(erg); %make them sum to 1

end

