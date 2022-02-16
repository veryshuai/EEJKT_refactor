function [err,jac,l] = sim_solve(x,bet,a,pi,n,j,ex_siz,rh,diag_Q,Q0_d,V_fail,V_succ,gam,scl,cscale,c,l_opt_func)
%This function creates and solves a set of ex_siz value functions for the search and learning Colombia Project.

    % scaled cost function
    %c = @(x,n) ((1+x).^(1+1/bet)-1)/((1+1/bet)*n^gam*scl/cscale);
    
    % Search policy, Zero if expected profit is less than zero
    l = l_opt_func(a(n,j),j,pi,V_succ,V_fail,x); 
%    l_check = max(max(((a(n,j)*(pi + V_succ)+(1-a(n,j))*V_fail-x)*j^gam*scl/cscale),0).^bet-1,0);
    l_check = max(max(((1+log(j))^gam*(a(n,j)*(pi+V_succ) + (1-a(n,j))*V_fail - x)/cscale)+1,0).^bet-1,0);
    
    if norm(l - l_check) > 1e-9
        disp('Warning:Policy check in sim_solve failed')
    end
    
    % Get error in current value function step 
    denom = rh+l+diag_Q;
    last_term = l.*(a(n,j)*(pi+V_succ)+(1-a(n,j))*V_fail);
    err = denom.^-1.*(-c(l,j)+Q0_d*x+last_term)-x;
    
    % Get jacobian
%     dl = (-bet*j^gam*((a(n,j)*(pi+V_succ)+(1-a(n,j))*V_fail-x)*j^gam*scl/cscale).^(bet-1)).*(l>0);
%     dc = (1+l).^(1/bet)/(j^gam*scl/cscale) .* dl;
%     dlast = (a(n,j)*(pi + V_succ)+(1-a(n,j))*V_fail).*dl;
%     dden = denom.^-2.*dl;
%     derr = dden.*(-c(l,j)+Q0_d*x+last_term)+denom.^-1.*(-dc+dlast)-1;
%     jac = bsxfun(@times,denom.^-1,Q0_d)';
%     jac(1:ex_siz+1:end) = derr'; 
%     jac = jac';
%     jac(isnan(jac)==1 | abs(jac) == inf) = 0; %filter out weird values
    jac = 0;
end %end function
