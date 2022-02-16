function [new_v,nostop,differ] = val_func_iter(init,bet,a,pi,net,Q0,rh,diag_Q,Q0_d,V_succ,V_fail,gam,scl,cscale,tol,rel_tol,no_learning,N,m,j,c,l_opt_func,pi_orig)
%This function performs value function iteration, which must be done several times in the EEJKT code 

     differ = 2;
     differ_rel = 2;
     old_v = init;
     max_damp = 0;
     iter = 0;
     nostop = 0;
     max_iter = 5e4;
     max_pi = max(pi_orig);

     while differ > tol/scl && differ_rel > rel_tol/scl && iter < max_iter
         
         iter = iter + 1;
         if iter == max_iter
             %display('WARNING: maximum iterations reached in value function iteration.  Results unreliable.  Value function difference stuck at:')
             %differ
             nostop = 1;
         end
         
         if no_learning == 1
            new_v = sim_solve_h(old_v,bet,a,pi,net,size(Q0,1),rh,diag_Q,Q0_d,V_succ,gam,scl,cscale,c,l_opt_func) + old_v;
         else
            new_v = sim_solve(old_v,bet,a,pi,N+1-m,j,size(Q0,1),rh,diag_Q,Q0_d,V_fail,V_succ,gam,scl,cscale,c,l_opt_func) + old_v;
         end
         
         %differ = max(abs(new_v - old_v));
         differ = sum(sum((new_v - old_v).*(new_v - old_v))); %this is the fsolve definition
         differ_rel = differ / max_pi;
         %differ = max(abs(new_v - old_v)./old_v)
         %damp = rand * damp;
         damp = min(max_damp,1 - 1/iter);
         old_v = damp * old_v + (1 - damp) * new_v;
     end

end
