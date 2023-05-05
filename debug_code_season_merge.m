pickID = [ 0 0 0 0 1 0 0 0 0 0];
pickIDyr = kron(ones(1,12),pickID);
test2a = all_seas.*(ones(size(all_seas,1),1)*pickIDyr);
firmIDa = max(test2a,[],2);
test2s = som_seas.*(ones(size(all_seas,1),1)*pickIDyr);
firmIDs = max(test2s,[],2);