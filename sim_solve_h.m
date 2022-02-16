function [err,jac,l] = sim_solve_h(x,bet,a,pi,n,ex_siz,rh,diag_Q,Q0_d,V_succ,gam,scl,cscale,c,l_opt_func)
%This function creates and solves a set of ex_siz value functions for the
%search and learning Colombia Project.

  % Scaled cost function
    %c = @(x,n) ((1+x).^(1+1/bet)-1)/((1+1/bet)*n^gam*scl/cscale);
    
  % Search policy, Zero if expected profit is less than zero
    l = l_opt_func(a,n,pi,V_succ,x,x);
 %  l_check = max((max((a*(pi + V_succ) + (1 - a) * x - x),0)*n^gam*scl/cscale).^bet-1,0);
    l_check = max(max(((1+log(n))^gam*(a*(pi+V_succ) + (1-a)*x - x)/cscale)+1,0).^bet-1,0);
    
    if norm(l - l_check) > 1e-9
        disp('Warning:Policy check in sim_solve_h failed')
    end
  % Get error in current value function step 
    denom = rh+l+diag_Q;
    last_term = l.*(a*(pi+V_succ)+(1-a)*x);
    err = denom.^-1.*(-c(l,n)+Q0_d*x+last_term)-x;
    
  % Calculate jacobian
%     dl = -a*bet * n^gam * scl / cscale * (a*(pi+V_succ-x)*n^gam*scl/cscale).^(bet-1).*(l>0); %okay, checked by me
%     dc = (1+l).^(1/bet)/(n^gam*scl/cscale) .* dl;
%     dlast = (a*(pi + V_succ)+(1-a)*x).*dl+(1-a)*l;
%     dden = - denom.^-2.*dl;
%     derr = dden.*(-c(l,n)+Q0_d*x+last_term)+denom.^-1.*(-dc+dlast)-1;
%     jac = bsxfun(@times,denom.^-1,Q0_d); 
%     jac(1:ex_siz+1:end) = derr; %note that the diagonals of Q0_d are zero
%     %jac = jac .* repmat(2 * x',size(x,1),1); % when using the squared method of ensuring x>0
%     jac(isnan(jac)==1 | abs(jac) == inf) = 0; %filter out weird values
jac = 0;

end
