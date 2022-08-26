function [new_v] = val_func_iter(learn_state,init,V_succ,V_fail,a,pi,net,mm,policy,country)

     differ = 2;
     differ_rel = 2;
     old_v = init;
     iter = 0;
     max_pi = retrieveCountrySpecificParameters(country,mm,policy);

     while differ > mm.v_tolerance && differ_rel > mm.v_rel_tolerance && iter < mm.max_iter
         
        iter = iter + 1;
        if learn_state == "no_learning"
            new_v = old_v + iterationError_h(old_v,a,pi,net,V_succ,mm,policy,country);
        elseif learn_state == "learning"
            new_v = old_v + iterationError_f(old_v,a,pi,net,V_fail,V_succ,mm);
        end
        differ = sum(sum((new_v - old_v).*(new_v - old_v))); 
        differ_rel = differ / max_pi;
        old_v = new_v;
         
     end

end
