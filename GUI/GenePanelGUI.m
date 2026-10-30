function GenePanelGUI(geneList)
    % GENEPANELGUI Interactive tool to assign markers and colors to a gene panel.
    % Features: Double-click color cell to open color picker, robust Excel loading, 
    % alphabetical table ordering, multi-column legend preview map.

    if nargin < 1 || isempty(geneList)
        defaultRaw = change_gene_symbols();
        geneList = defaultRaw(:, 1)';
    else
        if isstring(geneList)
            geneList = cellstr(geneList);
        end
    end

    try
        masterData = change_gene_symbols();
    catch
        masterData = {}; 
    end

    numGenes = length(geneList);
    tempData = cell(numGenes, 3);

    for i = 1:numGenes
        gName = char(geneList{i});
        tempData{i, 1} = gName;
        
        matchedIdx = [];
        if ~isempty(masterData)
            matchLog = strcmp(masterData(:, 1), gName);
            if any(matchLog), matchedIdx = find(matchLog, 1); end
        end
        
        if ~isempty(matchedIdx)
            rgb = masterData{matchedIdx, 2};
            tempData{i, 2} = sprintf('%.3f, %.3f, %.3f', rgb(1), rgb(2), rgb(3));
            tempData{i, 3} = char(masterData{matchedIdx, 3});
        else
            tempData{i, 2} = '0.700, 0.700, 0.700'; 
            tempData{i, 3} = 'o';                   
        end
    end

    % Sort rows by color family grouping initially
    colorStrings = tempData(:, 2);
    [~, sortIdx] = sort(colorStrings);
    tableData = tempData(sortIdx, :);

    % --- UI Figure Setup ---
    fig = uifigure('Name', 'Gene Panel Legend Customizer', ...
                   'Position', [50 50 1350, 800]);

    mainGrid = uigridlayout(fig, [1, 2]);
    mainGrid.ColumnWidth = {'2.5x', '3.5x'}; 

    % Left Sub-grid for Table & Control Buttons
    leftGrid = uigridlayout(mainGrid, [3, 3]);
    leftGrid.RowHeight = {'fit', '1x', 'fit'};
    leftGrid.ColumnWidth = {'1x', '1x', '1x'};

    lblTitle = uilabel(leftGrid, 'Text', 'Double-click any Color cell to open Color Picker', ...
        'FontSize', 12, 'FontWeight', 'bold');
    lblTitle.Layout.Row = 1;
    lblTitle.Layout.Column = [1 3];

    % UITable
    uit = uitable(leftGrid);
    uit.Layout.Row = 2;
    uit.Layout.Column = [1 3];
    uit.ColumnName = {'Gene Symbol', 'Color (R, G, B)', 'Marker Type'};
    uit.ColumnEditable = [true, true, true];
    uit.ColumnFormat = {[], [], {'+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h'}};
    uit.Data = tableData;

    % Control Buttons row
    btnReload = uibutton(leftGrid, 'Text', 'Reload & Sort');
    btnReload.Layout.Row = 3;
    btnReload.Layout.Column = 1;
    btnReload.ButtonPushedFcn = @(~,~) reloadDefaults();

    btnLoad = uibutton(leftGrid, 'Text', 'Load Excel');
    btnLoad.Layout.Row = 3;
    btnLoad.Layout.Column = 2;
    btnLoad.ButtonPushedFcn = @(~,~) loadFromExcel();

    btnSave = uibutton(leftGrid, 'Text', 'Save Excel');
    btnSave.Layout.Row = 3;
    btnSave.Layout.Column = 3;
    btnSave.ButtonPushedFcn = @(~,~) saveToExcel();

    % --- Right Side: Legend Plot with Black Title ---
    rightPanel = uipanel(mainGrid, 'Title', 'Full Legend Layout Preview (Grouped by Color)', 'FontSize', 12, 'FontWeight', 'bold');
    rightGrid = uigridlayout(rightPanel, [1, 1]);

    ax = uiaxes(rightGrid);
    ax.Layout.Row = 1;
    ax.Box = 'on';
    ax.Color = [0 0 0];
    ax.XColor = [0.2 0.2 0.2];
    ax.YColor = [0.2 0.2 0.2];
    
    tObj = title(ax, 'Complete Gene Panel Visual Map');
    tObj.Color = [0 0 0]; 
    tObj.FontSize = 13;

    % Hook cell selection / edits for color picking
    uit.CellEditCallback = @(~,~) handleManualEdit();
    uit.CellSelectionCallback = @(~,event) checkColorPicker(event);

    % Initial render
    redrawFullLegend();

    % --- Helper & Callback Functions ---

    function checkColorPicker(event)
        % If user clicks on a cell in Column 2 (Color column)
        if ~isempty(event.Indices) && event.Indices(2) == 2
            r = event.Indices(1);
            currentValStr = char(string(uit.Data{r, 2}));
            currentRgb = sscanf(currentValStr, '%f,%f,%f')';
            if length(currentRgb) ~= 3, currentRgb = [0.7 0.7 0.7]; end
            
            % Trigger MATLAB's color selection palette dialog
            selectedColor = uisetcolor(currentRgb, sprintf('Pick Color for: %s', char(string(uit.Data{r, 1}))));
            
            if numel(selectedColor) == 3
                uit.Data{r, 2} = sprintf('%.3f, %.3f, %.3f', selectedColor(1), selectedColor(2), selectedColor(3));
                handleManualEdit();
            end
        end
    end

    function handleManualEdit()
        currentData = uit.Data;
        cleanData = cell(size(currentData));
        for r = 1:size(currentData, 1)
            cleanData{r, 1} = char(string(currentData{r, 1}));
            cleanData{r, 2} = char(string(currentData{r, 2}));
            cleanData{r, 3} = char(string(currentData{r, 3}));
        end
        
        [~, sIdx] = sort(cleanData(:, 2));
        uit.Data = cleanData(sIdx, :);
        redrawFullLegend();
    end

    function redrawFullLegend()
        cla(ax);
        d = uit.Data;
        totalGenes = size(d, 1);
        if totalGenes == 0; return; end

        safeD = cell(size(d));
        for k = 1:totalGenes
            safeD{k, 1} = char(string(d{k, 1}));
            safeD{k, 2} = char(string(d{k, 2}));
            safeD{k, 3} = char(string(d{k, 3}));
        end

        colorStrings = safeD(:, 2);
        [~, colorSortIdx] = sort(colorStrings);
        plotData = safeD(colorSortIdx, :);

        rowsPerCol = max(10, ceil(totalGenes / 5)); 
        numCols = ceil(totalGenes / rowsPerCol);

        xSpacing = 2.8;  
        ySpacing = 1.1;  

        maxLenX = numCols * xSpacing + 0.5;
        maxLenY = rowsPerCol * ySpacing + 0.5;

        ax.XLim = [0 maxLenX];
        ax.YLim = [0 maxLenY];
        ax.XTick = []; 
        ax.YTick = [];
        ax.Color = [0 0 0];
        hold(ax, 'on');

        fSize = max(10, min(14, 17 - floor(totalGenes / 25)));

        idx = 1;
        for c = 1:numCols
            for r = rowsPerCol:-1:1 
                if idx > totalGenes; break; end
                
                gName = plotData{idx, 1};
                rgbStr = plotData{idx, 2};
                markerSymbol = plotData{idx, 3};
                
                rgbVal = sscanf(rgbStr, '%f,%f,%f')';
                if length(rgbVal) ~= 3, rgbVal = [1 1 1]; end
                rgbVal = min(max(rgbVal, 0), 1);

                xPos = (c - 1) * xSpacing + 0.4;
                yPos = (r - 1) * ySpacing + 0.5;

                plot(ax, xPos, yPos, markerSymbol, ...
                    'Color', rgbVal, ...
                    'MarkerFaceColor', rgbVal, ...
                    'MarkerSize', 8.5, ...
                    'LineWidth', 1.5);

                text(ax, xPos + 0.3, yPos, gName, ...
                    'Color', [1 1 1], ...
                    'FontSize', fSize, ...
                    'FontWeight', 'bold', ...
                    'Interpreter', 'none');

                idx = idx + 1;
            end
        end
        hold(ax, 'off');
        
        tObj = title(ax, sprintf('Complete Gene Panel Visual Map (%d Genes, %d Columns)', totalGenes, numCols));
        tObj.Color = [0 0 0];
        tObj.FontSize = 13;
    end

    function reloadDefaults()
        try
            raw_symbols = change_gene_symbols();
            currentData = uit.Data;
            for r = 1:size(currentData, 1)
                gName = char(string(currentData{r, 1}));
                matchLog = strcmp(raw_symbols(:, 1), gName);
                if any(matchLog)
                    matchedIdx = find(matchLog, 1);
                    rgb = raw_symbols{matchedIdx, 2};
                    currentData{r, 2} = sprintf('%.3f, %.3f, %.3f', rgb(1), rgb(2), rgb(3));
                    currentData{r, 3} = char(raw_symbols{matchedIdx, 3});
                end
            end
            
            cleanData = cell(size(currentData));
            for r = 1:size(currentData, 1)
                cleanData{r, 1} = char(string(currentData{r, 1}));
                cleanData{r, 2} = char(string(currentData{r, 2}));
                cleanData{r, 3} = char(string(currentData{r, 3}));
            end
            
            [~, sIdx] = sort(cleanData(:, 2));
            uit.Data = cleanData(sIdx, :);
            redrawFullLegend();
            uialert(fig, 'Refreshed and regrouped known genes by color!', 'Success', 'Icon', 'success');
        catch ME
            uialert(fig, ['Error loading function: ' ME.message], 'Function Error', 'Icon', 'error');
        end
    end

    function loadFromExcel()
        [file, path] = uigetfile({'*.xlsx';'*.xls'}, 'Select Gene Panel Excel File');
        if isequal(file, 0); return; end
        fullPath = fullfile(path, file);
        try
            T = readtable(fullPath, 'FileType', 'spreadsheet', 'ReadVariableNames', true);
            rawCell = table2cell(T);
            
            cleanData = cell(size(rawCell));
            for r = 1:size(rawCell, 1)
                for c = 1:size(rawCell, 2)
                    cleanData{r, c} = char(string(rawCell{r, c}));
                end
            end
            
            [~, sIdx] = sort(cleanData(:, 2));
            uit.Data = cleanData(sIdx, :);
            redrawFullLegend();
            uialert(fig, 'Successfully loaded and sorted gene table from Excel!', 'Success', 'Icon', 'success');
        catch ME
            uialert(fig, ['Could not read Excel file: ' ME.message], 'File Error', 'Icon', 'error');
        end
    end

    function saveToExcel()
        [file, path] = uiputfile({'*.xlsx', 'Excel Spreadsheet (*.xlsx)'}, 'Save Gene Panel As', 'GenePanelConfig.xlsx');
        if isequal(file, 0); return; end
        fullPath = fullfile(path, file);
        try
            currentData = uit.Data;
            cleanData = cell(size(currentData));
            for r = 1:size(currentData, 1)
                cleanData{r, 1} = char(string(currentData{r, 1}));
                cleanData{r, 2} = char(string(currentData{r, 2}));
                cleanData{r, 3} = char(string(currentData{r, 3}));
            end
            
            T = table(cleanData(:,1), cleanData(:,2), cleanData(:,3), ...
                'VariableNames', {'GeneSymbol', 'Color_RGB', 'MarkerType'});
            writetable(T, fullPath);
            uialert(fig, ['Successfully saved configurations to: ' fullPath], 'Export Success', 'Icon', 'success');
        catch ME
            uialert(fig, ['Error saving file: ' ME.message], 'Save Error', 'Icon', 'error');
        end
    end
end