Cmax = max(gridHaz(:,3));
for j = 1:Cmax
Cfinder = find(gridHaz(:,3)==j); % for j cum. successes 
margDist = gridHaz(Cfinder,:);
fprintf('\r cumulative successes: %0.1f\n', j)
figure(j+6)
hist(margDist); 
end
figure(7)
           surf(succ_grid,csucc_grid,tmod_haz)
        title('model-based predicted log match hazard')
        ylabel('log(cum. success)')
        xlabel('log(1+success rate)')
        zlabel('log match hazard')     