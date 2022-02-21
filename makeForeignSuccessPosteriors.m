function postSuccessProbs = makeSuccessPosteriors(mm)

postSuccessProbs = zeros(mm.n_size+1,mm.n_size+1); 
commonPrior = @(x,y,z) gamma(y+z)/(gamma(y)*gamma(z))*x.^(y-1).*(1-x).^(z-1);

possibleSuccProb = linspace(.001,.999,1000);
for j=1:mm.n_size+1 %trials
    for k=1:j %successes
        post_improp = binopdf(k-1,j-1,possibleSuccProb).*commonPrior(possibleSuccProb,mm.af,mm.bf);
        post = post_improp/sum(post_improp);
        postSuccessProbs(j,k) = possibleSuccProb * post';
    end
end

end

