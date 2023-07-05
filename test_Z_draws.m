test_erg = mm.erg_pz;
Nmatch = 4;
fail_1 = 0;
fail_2 = 0;

n = 1e5;
for i = 1:n
     test_1= new_vec_C(Nmatch,size(mm.Z,1),cumsum(test_erg));
     flag = sum(test_1,2)~=Nmatch;
     if flag == 1
          fail_1 = fail_1 + flag;
          test_2 = new_vec_C(Nmatch,size(mm.Z,1),cumsum(test_erg));
          fail_2 = fail_2 + sum(test_2,2)~=Nmatch;
          [test_1;test_2]
     end

end

fprintf('\r\n %3.0f problems out of %6.0f trials\n',[fail_1,n]);
fprintf('\r %3.0f failures on second try\n', fail_2);