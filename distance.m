function D = distance(X)

try

    format long;
 
    rng(80085,'twister');
    seed_crand(80085);
    
    mm = setModelParameters(X);
    % choose the firm type to use for spot checking
    mm.check_type = 99;
    
    policy = generatePolicyAndValueFunctions(mm);
    simMoms = simulateMomentsMain(policy,mm);
    [D,~] = calculateDistanceAndPrint(simMoms,mm,X);

catch

 fprintf('\r\n Failed to evaluate fit metric \n')
 D = 1e12;

end
