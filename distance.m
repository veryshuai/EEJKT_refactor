function D = distance(X)

try

    format long;
 
%      X = [-6.45092  -5.31392   0.18694   1.69881   0.66313  12.37231...
%            0.19393   1.55780   1.23497  -0.27172   8.55913  -7.31628]; % fit approx. 12.22

% X = [ -6.69905  -6.88534   0.22734   1.15397   0.40905  12.99814...
%        0.11426   1.88804   1.27705  -0.35441   7.44522  -9.96677];  % fit metric: 11.8512300863

% X = [-6.79757952990269,-6.69783952264302,0.182682463632116,1.65396826998034,...
%     0.674671051170727,12.9981423348683,0.114262951381515,1.60040440594609,...
%     1.27705129904746,-0.291910617912570,7.44522115406195,-9.96676755945645]; % problem vector

X = [-5.75356966903725	-5.25223507169887	0.164825015774296	1.17936463019703 ...
    0.378914704447539	14.3015666827662	0.106596797278845	1.58037006160446 ...
    1.30315909881952	-0.376856695185435	5.35823230259007	-7.36635544725547];  % problem vector
 
   
    rng(80085,'twister');
    seed_crand(80085);
    
    mm = setModelParameters(X);
    % choose the firm type to use for spot checking
    mm.check_type = 91;
    
    policy = generatePolicyAndValueFunctions(mm);
    simMoms = simulateMomentsMain(policy,mm);
    [D,~] = calculateDistanceAndPrint(simMoms,mm,X);

catch

 fprintf('\r\n Failed to evaluate fit metric \n')
 D = 1e12;

end
