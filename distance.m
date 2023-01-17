function D = distance(X)

%try

    format long;
 
 %X = [-10.0234360850160,-3.86822062058644,0.216000078376388,0.702389472395804,...
 %      0.374214942055985,11.6388730026429,0.164637004707622,0.907628494202198,...
 %      0.915592893141143,-0.254901443835102,8.77853609204470]; 
X = [-12.0234360850160,-8.86822062058644,0.216000078376388,0.702389472395804,...
       0.374214942055985,11.6388730026429,0.164637004707622,0.907628494202198,...
       5,-0.254901443835102,6.77853609204470,1]; 

    rng(80085,'twister');
    seed_crand(80085);
    
    mm = setModelParameters(X);
    % choose the firm type to use for spot checking
    mm.check_type = 108;
    
    policy = generatePolicyAndValueFunctions(mm);
    simMoms = simulateMomentsMain(policy,mm);
    [D,~] = calculateDistanceAndPrint(simMoms,mm,X);

%catch

 fprintf('\r\n Failed to evaluate fit metric \n')
 D = 1e12;

%end
