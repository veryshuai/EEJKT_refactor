function D = distance(X)

%try

    format long;
    % X = [-3.60529  -3.87941  0.23156  2.46487  1.88176  15.42634 0.38246  11.72211  1.38618  -1.21819  13.00238  -6.13506];
    % X = [-3.60529  -3.87941  0.23156  2.46487  1.88176  15.42634 0.38246  11.72211  1.88618  -1.21819  13.00238  -6.13506];
    % X = [-3.60529  0  0.23156  2.46487  1.88176  15.42634 0.38246  11  0.5  -3  13.00238  0];
    
    % X = [ -3.82908  -2.83652   0.15152   2.08200   1.59488   5.39672...
    %        0.29474   0.58545   1.74938  -1.44383  15.28438  -7.80303]; % bad X for debugging
    
     X = [-6.45092  -5.31392   0.18694   1.69881   0.66313  12.37231...
            0.19393   1.55780   1.23497  -0.27172   8.55913  -7.31628]; % fit approx. 12.22

   
    rng(80085,'twister');
    seed_crand(80085);
    
    mm = setModelParameters(X);
    policy = generatePolicyAndValueFunctions(mm);
    simMoms = simulateMomentsMain(policy,mm);
    [D,~] = calculateDistanceAndPrint(simMoms,mm,X);

%catch

%  fprintf('\r\n Failed to evaluate fit metric \n')
%  D = 1e12;

end
