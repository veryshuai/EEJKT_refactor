function [new_v] = val_func_iter_f(init,a,net,V_succ,V_fail,mm,policy)
   
     differ = 2;
     differ_rel = 2;
     old_v = init;
     iter = 0;
     max_pi = max(pi);

     while differ > mm.v_tolerance && differ_rel > mm.v_rel_tolerance && iter < mm.max_iter
         
        iter = iter + 1;
        new_v = old_v + iterationError_f(old_v,a,net,V_fail,V_succ,mm,policy);
        differ = sum(sum((new_v - old_v).*(new_v - old_v))); 
        differ_rel = differ / max_pi;
        old_v = new_v;
         
     end

end
