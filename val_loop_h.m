function [V,l_opt,my_flag,frac_nostop_h,worst_error,worst_scaled_error] = val_loop_h(Q0,Q0_d,a,pi,mm)
%This function solves for the value function and optimal search intensity 

    % fix original unscaled lifetime relationship profit expectation
    pi_orig = pi;
    
    % Read in parameters
    bet         = mm.b;               %cost parameter
    rh          = mm.r;               %discount factor
    dim0        = mm.dim0;            %size of th_g common
    dim1        = mm.dim1;            %size of th_h home
    net_size    = mm.net_size;        %max number of network effects
    gam         = mm.gam;             %network effect parameter
    tol         = mm.v_tolerance;     %tolerance for value function loop
    rel_tol     = mm.v_rel_tolerance; %relative tolerance for value function loop
    cscale      = mm.cs_h;            %counterfactual scaling parameter    
    c           = mm.cost_h;          %cost function
    l_opt_func  = mm.l_opt_func_h;    %search policy function
 
    
    % Initialize value function and policy function
    V = zeros(size(Q0,1),dim0,dim1,net_size+1); %[macro-prod states, theta_g,theta_h, max network effect]
    l_opt = zeros(size(V)); %optimal search intensities
    
    % Scale the profit expectation 
    %scl = max(pi_orig);
    scl=1;
%   c  = c / scl; won't work: can't scale a function
    pi = pi_orig/scl;
    
    % Scaled cost function
    %c = @(x,n) ((1+x).^(1+1/bet)-1)/((1+1/bet)*n^gam*scl/cscale);
%   aexp = 1+1/bet;
%   c = @(x,n) ((cscale/scl)/aexp)*((1+x).^aexp-(1 + aexp*x))/(n^gam);
    
%   l_opt_func = @(a,net,pi,V_succ,V_fail,V_orig) max(max((net^gam*(a*(pi + V_succ) + (1 - a) * V_fail - V_orig)*scl/cscale),0).^bet-1,0);       

	%l_opt_func = @(a,net_size,pi,V_succ,V_fail) max(((net_size+1)^gam*(a*(pi + V_succ - V_fail))*scl/cscale + 1).^bet-1,0);
                %l_opt(:,N+1,j,k,net+1) = max(((net+1)^gam*a(N+1,j,k)*pi*scl/cscale).^bet-1,0); 
                
    % Initialize punishment for non-convergence
    punishment = 0; 
    times_nostop = 0;
    iter = 0;
    worst_scaled_error = [0;0;0]; %err, max_pi, err/max_pi
    worst_error = [0;0;0]; %err, max_pi, err/max_pi
    
    % Make sure that the Q matrix (intensity matrix for events) is positive 
    diag_Q = abs(diag(Q0));
    
    % Options for solver
     % options = optimset('Display','off','Jacobian','on','GradObj','on',...
     %         'GradConstr','on','TolX',tol/100/scl,'TolFun',tol/scl^2,...
     %         'DerivativeCheck','off','maxiter',10000,'Algorithm','trust-region-dogleg');
    % Setup break flag
    my_flag= 0;
    
    % Start timer
    tic;

    % Backward induction to solve for ;value
    for j = 1:dim0
        tic
        for k = 1:dim1
            
            % first backwards induction 
        	%l_opt(:,j,k,net_size+1) = max(((net_size+1)^gam*a(j,k)*pi*scl/cscale).^bet-1,0);
            l_opt(:,j,k,net_size+1) = l_opt_func(a(j,k),net_size,pi,0,0,0);
            den = rh+diag_Q;
            perm_flow = -c(l_opt(:,j,k,net_size+1),net_size+1)+a(j,k)*l_opt(:,j,k,net_size+1).*pi;
            Q0_diag = -Q0_d;
            Q0_diag(1:size(Q0_d,1)+1:end) = den;
                % next eq. follows from V = (rh + diag_Q)^-1 * (perm_flow + Q0_d * V)
            V(:,j,k,net_size+1) = max(Q0_diag^-1*perm_flow,0);     
            
            % subsequent steps 
            for l = 1:net_size

                % use the previous iterations solution as initial guess 
                %init = ones(size(V(:,j,k,net_size-l+2)));
                init = V(:,j,k,net_size-l+2); 

%                 % solve for value using fsolve
                 %[res,~,flag_update] = fsolve(@(x) sim_solve_h(x,bet,a(j,k),pi,net_size+1-l,size(Q0,1),rh,diag_Q,Q0_d,V(:,j,k,net_size-l+2),gam,scl,cscale,c,l_opt_func),init,options);
                                  %V(:,j,k,net_size-l+1) = res;
                 %                 [~,~,l_opt_fsolve] = sim_solve_h(res,bet,a(j,k),pi,net_size+1-l,size(Q0,1),rh,diag_Q,Q0_d,V(:,j,k,net_size-l+2),gam,scl,cscale,c,l_opt_func);
%                 
                 % check for convergence
%                 if flag_update<1
%                     display('WARNING: home value function convergence issue!');
%                     my_flag = 1;
%                     V = 0;
%                     L_opt = 0;
%                     return
%                 end

                %value function iteration
                [new_v,nostop,differ] = val_func_iter(init,bet,a(j,k),pi,net_size+1-l,Q0,rh,diag_Q,Q0_d,V(:,j,k,net_size-l+2),1,gam,scl,cscale,tol,rel_tol,1,[],[],[],c,l_opt_func,pi_orig);
                %max(abs(res - new_v)) %DEBUG: any difference between val func iteration and fsolve?  Apparently not.
                times_nostop = times_nostop + nostop;
                iter = iter + 1;
                
                if differ/max(pi_orig) > worst_scaled_error(3)
                    worst_scaled_error = [differ,max(pi_orig),differ/max(pi_orig)];               
                end
                if differ > worst_error(1)
                    worst_error = [differ,max(pi_orig),differ/max(pi_orig)];
                end
            
                V(:,j,k,net_size-l+1) = new_v;
                
                
                
                % get the policy function
                [~,~,l_opt(:,j,k,net_size+1-l)] = sim_solve_h(V(:,j,k,net_size+1-l),bet,a(j,k),pi,net_size+1-l,size(Q0,1),rh,diag_Q,Q0_d,V(:,j,k,net_size-l+2),gam,scl,cscale,c,l_opt_func);
                                                
                if toc > 120
                    disp("ERROR:Home Value function iteration taking too long")
                    my_flag = 1;
                    return
                end 
            
            end
        end
    end
    frac_nostop_h = times_nostop / iter;
    %display(frac_nostop_h);
end

