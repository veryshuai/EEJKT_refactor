function [err,jac,l] = iterateValueFunction_f(x,a,pi,tri,succ,rh,diag_Q,Q0_d,V_fail,V_succ,c,l_opt_func)

    l = l_opt_func(a(tri,succ),succ,pi,V_succ,V_fail,x);  
    denom = rh+l+diag_Q;
    last_term = l.*(a(tri,succ)*(pi+V_succ)+(1-a(tri,succ))*V_fail);
    err = denom.^-1.*(-c(l,succ)+Q0_d*x+last_term)-x;
    
    jac = 0;
end 
