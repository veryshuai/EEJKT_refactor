function [Q0,Q0_d,st] = makebigq(Q_p,Q_f,n_p,n_f,Phi,X_f)
%get the exogenous state transition matrix

%make correct sizes
n_p = 2*n_p+1;
n_f = 2*n_f+1;

st = zeros(2,n_p*n_f);

for k = 1:n_p
    for m = 1:n_f
            st(:,(k-1)*n_f+m) = [Phi(k);X_f(m)];
    end
end

Q0 = zeros(n_p*n_f); % square matrix of zeros

for rp = 1:n_p
    for rf = 1:n_f
        for cp = 1:n_p
            for cf = 1:n_f
                if (rp == cp && rf == cf)
                    Q0((rp-1)*n_f + rf,(cp-1)*n_f + cf) = Q_p(rp,cp)+Q_f(rf,cf);
                elseif rp == cp
                    Q0((rp-1)*n_f + rf,(cp-1)*n_f + cf) = Q_f(rf,cf);
                elseif rf == cf
                    Q0((rp-1)*n_f + rf,(cp-1)*n_f + cf) = Q_p(rp,cp);
                else
                    Q0((rp-1)*n_f + rf,(cp-1)*n_f + cf) = 0;
                end
            end
        end
    end
end 

%Create Q0 but with zeros on the diagonal for later use 
Q0_d = Q0;
j = size(Q0,1);
Q0_d(1:(j+1):end) = 0; 

end
