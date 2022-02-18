function [lambda_f,lambda_h,pi_h,pi_f,c_val_h,c_val_f,my_flag,value_h,value_f] = solve_v1(mm)

    policy = struct();
        
    [policy.pi_f,policy.c_val_f] = solveExpectedMatchProfit(mm.scale_f,mm.X_f,mm.Q_f,mm.Q_f_d,mm.F_f,mm); %pi[prod,macro]
    [policy.pi_h,policy.c_val_h] = solveExpectedMatchProfit(mm.scale_h,mm.X_h,mm.Q_h,mm.Q_h_d,mm.F_h,mm); %c_val[prod,macro,demand shk]

    
    %% CALCULATE POSTERIOR MATCH PROBABILITIES
    
    postSuccessProb_f = makeForeignSuccessPosteriors(mm); %[trials,successes]
    
    %% SOLVE FOR SEARCH EFFORTS AT HOME AND IN FOREIGN MARKET
    
     [valh,lh,flag_h] = solvePolicyHome(policy,mm);  
%     [~,lf,flag_f] = val_loop_f(Q0_f,Q0_f_d,a_f,pi_f,mm); 
 
% parallelize the foreign and domestic policy function calculations

    flag_test = zeros(2,1);
    ll = cell(2);
%  for nnn=1:2
    parfor nnn=1:2
        if nnn==1
    [val{nnn},ll{nnn},flag_test(nnn)] = solvePolicyHome(policy,mm); 
    % args of ll(1): [home macro & prod state, general type (not used),theta_h, network effect]
        elseif nnn==2
    [val{nnn},ll{nnn},flag_test(nnn)] = val_loop_f(Q0_f,Q0_f_d,postSuccessProb_f,pi_f,mm); 
     % args of ll(2): [foreign macro & prod state,trials,successes,general type(not used), network effect]
        end
    end
    lh = ll{1};
    lf = ll{2};
    flag_h = flag_test(1);
    flag_f = flag_test(2);
    val_h = val{1};
    val_f = val{2};
     
    
    my_flag = max(flag_h,flag_f); %flag = 1 means value function optimization did not end normally, flag = 0 means the loop ended normally.

    %For use in "no learning" counterfactual
    %punishment = 0;
    %lf = zeros(size(Q0_f,1),n_size+1,n_size+1,dim0,net_size+1);
    
    %% RESHAPE 
    % For now, put matrix result into arrays to match old estimation routine format 
    if my_flag == 0
    
        % reshape policy functions to split up macro and productivity shocks.
        % In practice, the reshape function below breaks the first dimension of the policy
        % function (combined productivity and macro shocks, which are organized as [(p1,m1);(p1,m2);...(pN,mN-1),(pN,mN)]
        % into two dimensions, the first for macro shocks, the second for
        % productivity.  Other dimensions continue to indicate general
        % product appeal, specific product appeal, and network size for home, and
        % trials, successes, general product appeal, and network size for
        % foreign.
        lh_new = reshape(lh,n_x,n_phi,dim0,dim1,net_size+1); %splitting up the exogenous shock vector to macro and productivity shocks
        lf_new = reshape(lf,n_x,n_phi,n_size+1,n_size+1,dim0,net_size+1); %same
        val_h_new = reshape(val_h,n_x,n_phi,dim0,dim1,net_size+1);
        val_f_new = reshape(val_f,n_x,n_phi,n_size+1,n_size+1,dim0,net_size+1);
 
       % lh = [home macro & prod state, general type (not used),theta_h, network effect] reshaped with dimensions: 
       % # macro shocks, # seller productivities, # gen types, # theta_h's, # network effects 
       
       % lf = [foreign macro & prod state,trials,successes,general type(not used), network effect] reshaped with dimensions:  
       % # macro shocks, # seller productivities, # network effects, # network effects, general type, # network effects 
       
       % read policy function into cell arrays
       % Each element of the cell arrays have two dimensions.  The rows are
       % productivity states, and the columns are macro states.
       lambda_h = cell(dim0,dim1,net_size); %(general product appeal,specific product appeal,net effect)
       value_h = cell(dim0,dim1,net_size); %(general product appeal,specific product appeal,net effect)
       lambda_f = cell(n_size+1,n_size+1,dim0,net_size+1); %(successes,trials,general product appeal,net effect)
       value_f = cell(n_size+1,n_size+1,dim0,net_size+1);

       for j = 1:n_size+1
           for k = 1:n_size+1
               for m = 1:dim0
                   for n = 1:net_size+1
                   lambda_f{k,j,m,n} = lf_new(:,:,j,k,m,n)';
                   value_f{k,j,m,n} = val_f_new(:,:,j,k,m,n)';
                   end
               end
           end
       end
       for j = 1:dim0
           for k = 1:dim1
               for l = 1:net_size+1
               lambda_h{j,k,l} = lh_new(:,:,j,k,l)'; 
               value_h{j,k,l} = val_h_new(:,:,j,k,l)'; 
               end
           end
       end
        
        %Also reshape expected payoffs to a match
        %Rows are productivity, columns are macro shocks
        pi_f = reshape(pi_f,n_x,n_phi)';
        pi_h = reshape(pi_h,n_x,n_phi)';
        
        %reshape match continuation value functions
        %the cell is over demand shocks
        %each element of the shell is row productivity shock and column macro
        %shock.
        c_val_f_cell = cell(n_z,1);
        c_val_h_cell = cell(n_z,1);
        for k = 1:n_z
           c_val_f_cell{k,1} = reshape(c_val_f(:,k),n_x,n_phi)';
           c_val_h_cell{k,1} = reshape(c_val_h(:,k),n_x,n_phi)';
        end
        c_val_f = c_val_f_cell;
        c_val_h = c_val_h_cell;

    else %IN CASE OF ERROR
        lambda_f = 1;lambda_h = 1;pi_h = 1;pi_f = 1;c_val_h = 1;c_val_f = 1; value_f=1; value_h=1;
        'policy function error in solve_v1.m'
    end

end
