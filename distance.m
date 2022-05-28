function [D,W,error] = distance(X)

format long;
X = [-3.60529  -3.87941  0.23156  2.46487  1.88176  15.42634 0.38246  11.72211  1.38618  -1.21819  13.00238  -6.13506];
% X = [-3.60529  -3.87941  0.23156  2.46487  1.88176  15.42634 0.38246  11.72211  1.38618  -1.21819  13.00238  -5.13506];

rng(80085,'twister');
seed_crand(80085);

mm = setModelParameters(X);
policy = generatePolicyAndValueFunctions(mm);
simMoms = simulateMomentsMain(policy,mm);
calculateDistanceAndPrint(simMoms,mm,X);

end 
