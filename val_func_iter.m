function [new_v] = val_func_iter(init,bet,a,pi,net,Q0,rh,diag_Q,Q0_d,V_succ,V_fail,gam,scl,cscale,tol,rel_tol,no_learning,N,m,j,c,l_opt_func,pi_orig)

     differ = 2;
     differ_rel = 2;
     old_v = init;
     max_damp = 0;
     iter = 0;
     max_iter = 5e4;
     max_pi = max(pi_orig);

     while differ > tol/scl && differ_rel > rel_tol/scl && iter < max_iter

         if no_learning == 1
            new_v = old_v + iterateValueFunction_h(old_v,a,pi,net,rh,diag_Q,Q0_d,V_succ,c,l_opt_func);
         else
            new_v = old_v + iterateValueFunction_f(old_v,a,pi,N+1-m,j,rh,diag_Q,Q0_d,V_fail,V_succ,c,l_opt_func);
         end
         
         if isnan(new_v)
             return;
         end

         differ = sum(sum((new_v - old_v).*(new_v - old_v))); 
         differ_rel = differ / max_pi;
         damp = min(max_damp,1 - 1/iter);
         old_v = damp * old_v + (1 - damp) * new_v;
         iter = iter + 1;

         if isnan(old_v)
             return;
         end
     end

end
