function postSuccessProbs_new = makeForeignSuccessPosteriors(mm)

postSuccessProbs = zeros(mm.n_size+1,mm.n_size+1); 
for j=1:mm.n_size+1 %trials
    for k = 1:j %successes
        % postSuccessProbs_new(j,k) = (mm.af + (k-1)) / (mm.af + mm.bf + (j-1));
        postSuccessProbs(j,k) = (mm.af + mm.optimism + (k-1)) / (mm.af + mm.optimism + mm.bf + (j-1));
    end
end

% %FOR ARBITRARY PRIOR DISTRIBUTION
% %commonPrior = @(x,y,z) gamma(y+z)/(gamma(y)*gamma(z))*x.^(y-1).*(1-x).^(z-1);
% commonPrior = @(x,y,z) gamma(y+mm.optimism+z)/(gamma(y+mm.optimism)*gamma(z))*x.^(y+mm.optimism-1).*(1-x).^(z-1);
% 
% possibleSuccProb = linspace(.0001,.9999,10000);
% for j=1:mm.n_size+1 %trials
%     for k=1:j %successes
%         post_improp = binopdf(k-1,j-1,possibleSuccProb).*commonPrior(possibleSuccProb,mm.af,mm.bf);
%         post = post_improp/sum(post_improp);
%         postSuccessProbs(j,k) = possibleSuccProb * post'; %expectation conditional on j trials, k succ
%     end
% end


end

