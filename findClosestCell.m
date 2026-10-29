function [closestIdx, closestCoord, clickXY] = findClosestCell(CellCalled)
    % Prompts user to click on the current plot and finds the nearest cell.
    % 
    % Assumptions:
    %   CellCalled.CellYX is an N x 2 matrix: [Y, X]
    %
    % Outputs:
    %   closestIdx   - Row index of the nearest cell
    %   closestCoord - Coordinates of the nearest cell: [Y, X]
    %   clickXY      - Coordinates of the user click: [X, Y]

    disp('Click on the figure to select a location...');
    
    % 1. Get click coordinates from the active plot
    [xClick, yClick] = ginput(1);
    clickXY = [xClick, yClick];

    % 2. Extract cell coordinates (assuming Col 1 = Y, Col 2 = X)
    cellY = CellCalled.CellYX(:, 1);
    cellX = CellCalled.CellYX(:, 2);

    % 3. Calculate Euclidean distances
    squaredDistances = (cellX - xClick).^2 + (cellY - yClick).^2;

    % 4. Find the minimum distance index
    [minDistSq, closestIdx] = min(squaredDistances);
    closestCoord = CellCalled.CellYX(closestIdx, :);

    % Print the result to command window
    fprintf('Clicked at: (X=%.2f, Y=%.2f)\n', xClick, yClick);
    fprintf('Closest cell #%d at: (Y=%.2f, X=%.2f) [Distance: %.2f]\n', ...
        closestIdx, closestCoord(1), closestCoord(2), sqrt(minDistSq));

    % Optional: Highlight on the current plot
    hold on;
    plot(xClick, yClick, 'rx', 'MarkerSize', 12, 'LineWidth', 2); % Red X = Click
    plot(closestCoord(2), closestCoord(1), 'go', 'MarkerSize', 10, 'LineWidth', 2); % Green Circle = Nearest cell
end