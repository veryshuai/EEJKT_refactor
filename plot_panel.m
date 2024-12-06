function plot_panel(subplot_idx, data1, data2, title_text)
    subplot(2, 2, subplot_idx); % Create a subplot

    % Plot first dataset
    scatter(data1(:, 2), data1(:, 1), 'filled', 'MarkerFaceColor', 'b'); % Blue filled circles
    hold on;

    % Plot second dataset
    scatter(data2(:, 2), data2(:, 1), 'd', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none'); % Red diamonds

    % Calculate the combined range of x and y data
    all_data = [data1; data2];
    x_min = min(all_data(:, 2));
    x_max = max(all_data(:, 2));
    y_min = min(all_data(:, 1));
    y_max = max(all_data(:, 1));
    overall_min = min(x_min, y_min);
    overall_max = max(x_max, y_max);

    % Add padding as a fraction of the range
    padding = 0.1; % Adjust if needed
    range = overall_max - overall_min;
    axis_min = overall_min - padding * range;
    axis_max = overall_max + padding * range;

    % Set equal limits for both x and y axes
    xlim([axis_min, axis_max]);
    ylim([axis_min, axis_max]);

    % Plot the red dotted line with slope 1
    plot([axis_min, axis_max], [axis_min, axis_max], 'r--', 'LineWidth', 1.5);

    hold off;
    title(title_text);
    xlabel('Simulated Moments');
    ylabel('Data Moments');
end
