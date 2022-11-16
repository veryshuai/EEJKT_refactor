function D = distance(X)

try

    format long;

% X = [ -6.69905  -6.88534   0.22734   1.15397   0.40905  12.99814...
%        0.11426   1.88804   1.27705  -0.35441   7.44522  -9.96677];  % fit metric: 11.8512300863

    rng(80085,'twister');
    seed_crand(80085);
    
    mm = setModelParameters(X);
    % choose the firm type to use for spot checking
    mm.check_type = 69;
    
    policy = generatePolicyAndValueFunctions(mm);
    simMoms = simulateMomentsMain(policy,mm);
    [D,~] = calculateDistanceAndPrint(simMoms,mm,X);

catch

 fprintf('\r\n Failed to evaluate fit metric \n')
 D = 1e12;

end
