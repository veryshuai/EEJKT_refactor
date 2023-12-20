function D = distance_nolearning(X)

try

    format long;
 
    rng(80085,'twister');
    seed_crand(80085);
    
    mm = setModelParameters(X);
    % choose the firm type to use for spot checking
    mm.check_type = 99;
    
    % create a different set of policies for each theta type
    % (This is slightly inefficient because we solve for identical
    % home policy functions every time)
    policy = generatePolicyAndValueFunctions(mm);
    simMoms = simulateMomentsMain_nolearning(policy,mm);
    [D,~] = calculateDistanceAndPrint(simMoms,mm,X);

catch

 fprintf('\r\n Failed to evaluate fit metric \n')
 D = 1e12;

end
