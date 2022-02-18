function [err,jac,l] = sim_solve_h(x,bet,a,pi,n,ex_siz,rh,diag_Q,Q0_d,V_succ,gam,scl,cscale,c,l_opt_func)

    l = l_opt_func(a,n,pi,V_succ,x,x);  
    denom = rh+l+diag_Q;
    last_term = l.*(a*(pi+V_succ)+(1-a)*x);
    err = denom.^-1.*(-c(l,n)+Q0_d*x+last_term)-x;
    
jac = 0;

end
