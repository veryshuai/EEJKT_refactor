function null = policy_plot(str,lambda,diff,lambda_new)
% This function plots foreign policy functions for the counterfactual simulations

    null = 0;
    
    % The original lambda isn't quite plotable as it is
    if diff == 0
        %successes vs failures
        plotable = zeros(size(lambda,1),size(lambda,2)) * NaN;
        for j = 1:size(lambda,1)
            for k = 1:j
                plotable(j - k + 1,k) = lambda{k,j,2,k}(10,9);
            end
        end
    
        fig = figure(1);
        surf(plotable)
        title(str)
        zlabel('search intensity')
        xlabel('succeses')
        ylabel('failures')
        print(fig,'-dpng',['results/', str, '_policy.png'])
    
        % Alternative 2d figure
        fig = figure(2);
        plot(plotable(1,:)','b');
        hold on;
        plot(plotable(6,:)','r--');
        hold on;
        plot(plotable(11,:)','g-.');
        hold on;
        plot(plotable(16,:)','k*');
        hold off;
        xlabel('succeses')
        ylabel('search intensity')
        set(gca,'XTick',[1:3:25]);
        set(gca,'XTickLabel',{'0','3','6','9','12','15','18','21','24'});
        legend('0 failures','5 failures','10 failures','15 failures','Location','SouthEast')
        print(fig,'-dpng',['results/', str, '_policy_2d.png'])
        
        %productivity vs trials, no successes
        plotable = zeros(size(lambda,1),size(lambda{1,1,1,1},2)) * NaN;
        for j = 1:size(lambda,1)
            for k = 1:size(lambda{1,1,1,1},2)
                suc = 1;
                plotable(j,k) = lambda{suc,j,1,suc}(k,10);
            end
        end
    
        fig = figure(3);
        surf(plotable)
        title([str,' no successes'])
        zlabel('search intensity')
        xlabel('productivity')
        ylabel('trials')
        az = 145;
        el = 50;
        view(az,el);
        print(fig,'-dpng',['results/', str, '_policy_no_suc.png'])
    
        % Alternative 2d figure
        fig = figure(4);
        plot(plotable(:,8),'b');
        hold on;
        plot(plotable(:,10),'r--');
        hold on;
        plot(plotable(:,12),'g-.');
        hold on;
        plot(plotable(:,14),'k*');
        hold off;
        xlabel('trials')
        ylabel('search intensity')
        legend('low productivity','mid low productivity','mid high productivity','high productivity','Location','SouthEast')
        print(fig,'-dpng',['results/', str, '_policy_no_suc_2d.png'])
    
    
    elseif diff == 1
        plotable = zeros(size(lambda,1),size(lambda,2)) * NaN;
        for j = 1:size(lambda,1)
            for k = 1:j
                plotable(j - k + 1,k) = lambda_new{k,j,2,k}(8,8) - lambda{k,j,2,k}(8,8);
            end
        end
        fig = figure(1);
        surf(plotable)
        title([str, ' policy difference'])
        zlabel('search intensity')
        xlabel('succeses')
        ylabel('failures')
        print(fig,'-dpng',['results/', str, '_policy_diff'])
    
        % Alternative 2d figure
        fig = figure(2);
        plot(plotable(1,:)','b');
        hold on;
        plot(plotable(6,:)','r--');
        hold on;
        plot(plotable(11,:)','g-.');
        hold on;
        plot(plotable(16,:)','k*');
        hold off;
        xlabel('succeses')
        ylabel('search intensity')
        set(gca,'XTick',[1:3:25]);
        set(gca,'XTickLabel',{'0','3','6','9','12','15','18','21','24'});
        legend('0 failures','5 failures','10 failures','15 failures','Location','SouthEast')
        print(fig,'-dpng',['results/', str, '_policy_diff_2d.png'])
    
        %productivity vs trials, no successes
        plotable = zeros(size(lambda,1),size(lambda{1,1,1,1},2)) * NaN;
        for j = 1:size(lambda,1)
            for k = 1:size(lambda{1,1,1,1},2)
                suc = 1;
                plotable(j,k) = lambda_new{suc,j,1,suc}(6,k)-lambda{suc,j,1,suc}(6,k);
            end
        end
    
        fig = figure(3);
        surf(plotable)
        title([str,' policy difference, no successes'])
        zlabel('search intensity')
        xlabel('productivity')
        ylabel('trials')
        az = 145;
        el = 50;
        view(az,el);
        print(fig,'-dpng',['results/', str, '_policy_no_suc_diff.png'])
    
        % Alternative 2d figure
        fig = figure(4);
        plot(plotable(:,1),'b');
        hold on;
        plot(plotable(:,5),'r--');
        hold on;
        plot(plotable(:,10),'g-.');
        hold on;
        plot(plotable(:,15),'k*');
        hold off;
        xlabel('trials')
        ylabel('search intensity')
        legend('low productivity','mid low productivity','mid high productivity','high productivity','Location','SouthEast')
        print(fig,'-dpng',['results/', str, '_policy_no_suc_2d.png'])
    end
end
