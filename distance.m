function D = distance(X)

try
format long;

% X = [-6.35647  -5.27227   0.21626   1.56615   1.18474  13.08123 0.25626   3.64325   1.46860  -0.43414   9.10446  -9.30786]; 

   fprintf('\r\n begin evaluation of X = ');
   fprintf('\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',X(1:6));
   fprintf('\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',X(7:12));
   fprintf( '\r\n  ');   

rng(80085,'twister');
seed_crand(80085);

mm = setModelParameters(X);
policy = generatePolicyAndValueFunctions(mm);
simMoms = simulateMomentsMain(policy,mm);
D = calculateDistanceAndPrint(simMoms,mm,X);
catch
    fprintf('\r\n Failed to evaluate fit metric \n'); 
     D = 1e12;
end

end  
