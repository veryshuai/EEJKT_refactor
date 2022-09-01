function [new_v] = val_func_iter(learn_state,init,V_succ,V_fail,a,pi,net,mm,policy,country)

     differ = 2;
     old_v = init;
     iter = 0;

     while differ > mm.v_tolerance && iter < mm.max_iter
         
        iter = iter + 1;
        if learn_state == "no_learning"
            new_v = old_v + iterationError_h(old_v,a,pi,net,V_succ,mm,policy,country);
        elseif learn_state == "learning"
            new_v = old_v + iterationError_f(old_v,a,pi,net,V_fail,V_succ,mm);
        end
        differ = norm(new_v - old_v)/norm(old_v); 
        old_v = 0.0 * old_v + 1.0 * new_v;
         
     end

end
