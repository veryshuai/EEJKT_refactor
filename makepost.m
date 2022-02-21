function a = makepost(N,th_g,af,bf,alp)
%This function creates posterior success probabilties

a = zeros(N+1,N+1,size(th_g,2)); %rows are number of trials, columns are number of successes (1st row is zero trials, 1st column is zero successes)

% pdf of theta f distribution
h = @(x,y,z) gamma(y+z)/(gamma(y)*gamma(z))*x.^(y-1).*(1-x).^(z-1);

% Get the posterior success expectations for each success failure mixture 
th_f = linspace(.001,.999,1000);
for m = 1:size(th_g,2) 
    for j=1:N+1
        for k=1:j
            pri = binopdf(k-1,j-1,(th_g(m))*alp+th_f*(1-alp)).*h(th_f,af,bf);
            pri = pri/sum(pri);
            a(j,k,m) = th_g(m) * alp + (th_f*(1-alp))* pri';
        end
    end
end

end % end function
