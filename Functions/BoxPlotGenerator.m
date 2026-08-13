function BoxPlotGenerator(NTdata,ASDdata, titleName, labelName)
% BoxPlotGenerator: This function generates boxplots with the NT group is
% blue and the ASD group in red. The raw data is overlayed with jitter.

    % Input: 
        % NT: The input data should be a column vector
        % ASD: The input data should be a column vector
        % titleName: String containing the figure title
        % labelName: String containing the figure ylabel

    boxplot_data = [NTdata; ASDdata];
    group_labels = [ones(length(NTdata), 1); 1.5*ones(length(ASDdata), 1)];
    
    % figure; hold on;
    scatter(ones(length(NTdata),1), NTdata, 30, [0.6 0.8 1], 'filled', 'jitter', 'on', 'jitterAmount', 0.1, 'MarkerFaceAlpha', 0.7);
    scatter(2*ones(length(ASDdata),1), ASDdata, 30, [1 0.6 0.6], 'filled', 'jitter', 'on', 'jitterAmount', 0.1, 'MarkerFaceAlpha', 0.7);
    
    
    h1 = boxplot(boxplot_data, group_labels, 'Symbol', '');
    set(h1(:,1), 'Color', 'b', 'LineWidth', 1.5); 
    set(h1(:,2), 'Color', 'r', 'LineWidth', 1.5);
    boxes = findobj(gca, 'Tag', 'Box');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'NT', 'ASD'}, 'fontsize',12);
    ylabel(labelName);
    title(titleName);

end