function policy = generatePolicyAndValueFunctions(mm) 

    tic
    policy = solvePolicyMain(mm);    
    % lambda_f (succ, trial, common succ rate (defunct), network size, prod of firm, F macro shock) 
    % lambda_h (common succ rate (defunct), known succ rate, network size, prod of firm, H macro shock)
    % c_val (demand shock, prod of firm, macro shock)
        
    policy = makeExporterTransitionProbabilities(mm,policy);

    toc
end