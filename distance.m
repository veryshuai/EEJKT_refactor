function D = distance(X)

try

    format long;
% 
% X = [ -6.69905  -6.88534   0.22734   1.15397   0.40905  12.99814...
%        0.11426   1.88804   1.27705  -0.35441   7.44522  -9.96677];  % fit metric: 11.8512300863

% X = [-5.75356966903725 -5.25223507169887 0.164825015774296 1.17936463019703 ...
%     0.378914704447539 14.3015666827662 0.106596797278845 1.58037006160446 ...
%     1.30315909881952 -0.376856695185435 5.35823230259007 -7.36635544725547]; % fit: 13.723717

X =  [-8.53630  -7.02531   0.23560   0.87633   0.40764  12.76228...
      0.08514   1.51229   1.24776  -0.28887   8.37010  -9.03099]; % fit = 11.6206

    rng(80085,'twister');
    seed_crand(80085);
    
    mm = setModelParameters(X);
    % choose the firm type to use for spot checking
    mm.check_type = 90;
    
    policy = generatePolicyAndValueFunctions(mm);
    simMoms = simulateMomentsMain(policy,mm);
    [D,~] = calculateDistanceAndPrint(simMoms,mm,X);

catch

 fprintf('\r\n Failed to evaluate fit metric \n')
 D = 1e12;

end
