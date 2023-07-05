function out = policy_h(a,net,pi,V_succ,V_fail,V_orig,mm)
      out = max(max(((1+log(net))^mm.gam*(a*(pi+V_succ) + (1-a)*V_fail - V_orig)/mm.cs_h)+1,0).^(1/(mm.kappa1-1))-1,0);    
 end