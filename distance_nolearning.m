function D = distance(X)

try

    format long;
 
    rng(80085,'twister');
    seed_crand(80085);
    
    mm = setModelParameters(X);
    % choose the firm type to use for spot checking
    mm.check_type = 99;
    
    policy_cell = cell(mm.dim2,1);
    for th_ind = 1:mm.dim2
        mm.th_ind_temp = th_ind;
        policy_cell{th_ind} = generatePolicyAndValueFunctions(mm);
    end
    simMoms = simulateMomentsMain_nolearning(policy_cell,mm);
    [D,~] = calculateDistanceAndPrint(simMoms,mm,X);

catch

 fprintf('\r\n Failed to evaluate fit metric \n')
 D = 1e12;

end
