function [V,l_opt,my_flag,frac_nostop_f,worst_error,worst_scaled_error] = val_loop_f(Q0,Q0_d,a,pi,mm)
%                                      (Q0_f,Q0_f_d,a_f,pi_f,mm)
%                                 a_f = makepost(n_size,th_g,af,bf,alp)

%calculates value function and optimal search intensity for foreign

  % fix original unscaled lifetime relationship profit expectation
    pi_orig = pi; %for help debugging.
    
  % read in parameters
    bet         = mm.b;           % cost function (inverse) convexity parameter
    rh          = mm.r;           % discount parameter per period
    dim0        = mm.dim0;        % size of general th_g (now set to 1)
    N           = mm.n_size;      % maximum learning matches
    net         = mm.net_size;    % maximum network effects
    tol         = mm.v_tolerance; % tolerance for loop below
    rel_tol     = mm.v_rel_tolerance; % relative tolerance for loop below
    gam         = mm.gam;         % network effect parameter
    cscale      = mm.cs_f;        % cost scaling term
    c           = mm.cost_f;      % cost function
    l_opt_func  = mm.l_opt_func_f;% policy function
 
     
  % initialize punishment for non-convergence
    punishment = 0; 
    times_nostop = 0;
    iter = 0;
    worst_scaled_error = [0;0;0]; %err, max_pi, err/max_pi
    worst_error = [0;0;0]; %err, max_pi, err/max_pi
    
  % Inititialize value function and policy functions
    V = zeros(size(Q0,1),N+1,N+1,dim0,net+1); 
    % [macro & prod. state,trials,successes,general type(not used),max network effect]
    
    l_opt = zeros(size(V)); %optimal search intensities
    % [macro & prod. state,trials,successes,general type(not used),max network effect]
    
  % scaling
    %scl = max(pi_orig);
    scl = 1;
  %  c   = c / scl; % won't work: can't scale a function
    
  % scaled cost function
    %c = @(x,n) ((1+x).^(1+1/bet)-1)/((1+1/bet)*(n^gam)*scl/cscale);
%   aexp = 1+1/bet;
% 	c = @(x,n) ((cscale/scl)/aexp)*((1+x).^aexp-(1 + aexp*x))/(n^gam);

%	l_opt_func = @(a,net,pi,V_succ,V_fail,V_orig) max(max((net^gam*(a*(pi + V_succ) + (1 - a) * V_fail - V_orig)*scl/cscale),0).^bet-1,0);
                %l_opt(:,N+1,j,k,net+1) = max(((net+1)^gam*a(N+1,j,k)*pi*scl/cscale).^bet-1,0);  
                
    
  % scaled expectation of profit from a relationship
    pi = pi_orig/scl;

  % make sure intensity matrix diagonals are positive
    diag_Q = abs(diag(Q0));
    
  % % solver options
  %   options = optimset('Display','off','Jacobian','on','GradObj','on',...
  %           'GradConstr','on','TolX',tol/100/scl,'TolFun',tol/scl,...
  %           'DerivativeCheck','off','maxiter',10000,'Algorithm','trust-region-dogleg');

  % Setup break flag
    my_flag= 0;
    
  %start timer
  tic;

  % backward induction to solve value function
    for k = 1:dim0  % k indexes general type (not used: dim0=1)
        try
            %Solve for value function, last learning step.             
            for j = 1:N+1 % j-1 is the number of successes plus 1 
                % policy function when network effect maxed out and trials > N
                %l_opt(:,N+1,j,k,net+1) = max(((net+1)^gam*a(N+1,j,k)*pi*scl/cscale).^bet-1,0);  
				l_opt(:,N+1,j,k,net+1) = l_opt_func(a(N+1,j,k),net+1,pi,0,0,0); %max(((net+1)^gam*a(j,k)*pi*scl/cscale + 1).^bet-1,0);
                Q0_diag = -Q0_d;
                Q0_diag(1:size(Q0_d,1)+1:end) = rh+diag_Q;
                perm_flow = a(N+1,j,k)*l_opt(:,N+1,j,k,net+1).*pi-c(l_opt(:,N+1,j,k,net+1),net+1); % scaled expected net rev. flow
                % next eq. follows from V = (rh + diag_Q)^-1 * (perm_flow + Q0_d * V)
                V(:,N+1,j,k,net+1) = max(Q0_diag^-1*perm_flow,0);
     
              % learning maxed out
                for m = 1:net+1-j

                    % initial guess
                    init = V(:,N+1,j,k,net+2-m);

%                     % solve for value function
                     %[V(:,N+1,j,k,net+1-m),~,flag_update] = fsolve(@(x) sim_solve_h(x,bet,a(N+1,j,k),pi,net+1-m,size(Q0,1),rh,diag_Q,Q0_d,V(:,N+1,j,k,net+2-m),gam,scl,cscale),init,options);
% 
%                     %check for convergence
                     %if flag_update<1 
                     %    display('WARNING: foreign value function convergence issue!');
                     %    my_flag = 1;
                     %    V = 0;
                     %    L_opt = 0;
                     %    return
                     %end


                    
                    %value function iteration
                    new_v = val_func_iter(init,bet,a(N+1,j,k),pi,net+1-m,Q0,rh,diag_Q,Q0_d,V(:,N+1,j,k,net+2-m),1,gam,scl,cscale,tol,rel_tol,1,[],[],[],c,l_opt_func,pi_orig);
                    V(:,N+1,j,k,net+1-m) = new_v;
                    
                    % get policy function
                    [~,~,l_opt(:,N+1,j,k,net+1-m)] = sim_solve_h(V(:,N+1,j,k,net+1-m),...
                    bet,a(N+1,j,k),pi,net+1-m,size(Q0,1),rh,diag_Q,Q0_d,V(:,N+1,j,k,net+2-m),gam,scl,cscale,c,l_opt_func);
                    
                    if toc > 1000
                        disp("ERROR:Foreign value function iteration taking too long")
                        my_flag = 1;
                        return
                    end 
                
                end
            end
            
          % the rest of the value surface, 
            % N+1-m is current trial number, j is successes.
            for m = 1:N
                for j = 1:N+1-m

                    %use the previous iterations solution as initial guess
                    init = V(:,N+2-m,j,k,j); 

                      % solve for value;
                      %[res,~,flag_update] = fsolve(@(x) sim_solve(x,bet,a(:,:,k),pi,N+1-m,j,size(Q0,1),rh,diag_Q,Q0_d,V(:,N-m+2,j,k,j),V(:,N-m+2,j+1,k,j+1),gam,scl,cscale),init,options);
% 
%                     %check for convergence
%                     if flag_update<1 
%                         display('WARNING: foreign value function convergence issue!');
%                         my_flag = 1;
%                         V = 0;
%                         L_opt = 0;
%                         return
%                     end

                     %value function iteration
                     [new_v,nostop,differ] = val_func_iter(init,bet,a(:,:,k),pi,N+1-m,Q0,rh,diag_Q,Q0_d,V(:,N-m+2,j+1,k,j+1),V(:,N-m+2,j,k,j),gam,scl,cscale,tol,rel_tol,0,N,m,j,c,l_opt_func,pi_orig);
                     times_nostop = times_nostop + nostop;
                     iter = iter + 1;
                     %max(abs(new_v - res)) DEBUG: Difference between val
                     %func iteration and fsolve.  Small (depends on
                     %tolerance)
                     
                    if differ/max(pi_orig) > worst_scaled_error(3)
                        worst_scaled_error = [differ,max(pi_orig),differ/max(pi_orig)];   
                    end
                    if differ > worst_error(1)
                        worst_error = [differ,max(pi_orig),differ/max(pi_orig)];
                    end
                     
                    V(:,N+1-m,j,k,j) = new_v;
                    
                    % policy function
                    [~,~,l_opt(:,N+1-m,j,k,j)] = sim_solve(V(:,N+1-m,j,k,j),bet,a(:,:,k),pi,N+1-m,j,size(Q0,1),rh,diag_Q,Q0_d,V(:,N-m+2,j,k,j),V(:,N-m+2,j+1,k,j+1),gam,scl,cscale,c,l_opt_func);
                                                 
                    if toc > 1000
                        disp("ERROR:Foreign value function iteration taking too long")
                        my_flag = 1;
                        return
                    end 
                end
            end

            % add to punishment for convergence problems
        catch myerr% catch general errors
            display(myerr)
            display('WARNING: foreign value function convergence issue!');
            my_flag = 1;
            V = 0;
            L_opt = 0;
            return

        end
        
    end

    % end timer
    frac_nostop_f = times_nostop / iter;
    %display(frac_nostop_f)
    
end %end function

