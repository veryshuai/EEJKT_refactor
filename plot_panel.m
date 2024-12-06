function plot_panel(subplot_idx, data, title_text)
    subplot(2, 2, subplot_idx); % Create a subplot
    scatter(data(:, 2), data(:, 1), 'filled'); % Scatter plot
    hold on;

    % Calculate the combined range of x and y data
    x_min = min(data(:, 2));
    x_max = max(data(:, 2));
    y_min = min(data(:, 1));
    y_max = max(data(:, 1));
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
