Spotnb = 121229;
SpotScore = o.ompSpotScore;
SpotI = o.ompSpotIntensity2;
NeighbNonZeros = o.ompNeighbNonZeros;

ThisSpotScore = SpotScore(Spotnb);
ThisSpotI = SpotI(Spotnb);
ThisNeighbNonZeros = NeighbNonZeros(Spotnb);

figure
subplot(1,3,1)
histogram(SpotScore,100)
hold on
YL = ylim;
plot([ ThisSpotScore ThisSpotScore],YL,'-r')
subplot(1,3,2)
histogram(log10(1+SpotI),1000)
hold on
YL = ylim;
plot([ log10(1+ThisSpotI) log10(1+ThisSpotI)],YL,'-r')
subplot(1,3,3)
histogram(NeighbNonZeros,40)
hold on
YL = ylim;
plot([ ThisNeighbNonZeros ThisNeighbNonZeros],YL,'-r')

figure
scatter(SpotScore,SpotI,'.')
%% we could do a plot which would give the average color code of spots for a given gene, across spots inside a range of thresholds
figure
for i = 1:100
    Gene = o.GeneNames{i};
    ThisSpots = strcmp(o.GeneNames(o.ompSpotCodeNo),Gene);
    ThisSpotsScore = SpotScore(ThisSpots);
    ThisSpotsI = SpotI(ThisSpots);
    ThisSpotsNeighbNonZeros = NeighbNonZeros(ThisSpots);

    subplot(7,15,i)
    scatter(ThisSpotsScore,ThisSpotsI,'.')
    title(Gene)
    xlim([-10 20])
    ylim([-0.5 2.5])
end
%%
Gene = 'Chodl';
ThisSpots = strcmp(o.GeneNames(o.ompSpotCodeNo),Gene);

figure
R = -5:1:20;
for i=1:length(R)
    if i==1
        AboveS = ThisSpots & o.ompSpotScore<R(i) ;
    else
        AboveS = ThisSpots & o.ompSpotScore<R(i) & o.ompSpotScore>R(i-1);
    end

    C = squeeze(nanmean(o.ompSpotColors(AboveS,:,:),1));
    subplot(5,6,i)
    imagesc(C)
    caxis([-500 15000])
end


figure
R = -0.1:0.2:2.2;
for i=1:length(R)
    if i==1
        AboveS = ThisSpots & o.ompSpotIntensity2<R(i) ;
    else
        AboveS = ThisSpots & o.ompSpotIntensity2<R(i) & o.ompSpotIntensity2>R(i-1);
    end

    C = squeeze(nanmean(o.ompSpotColors(AboveS,:,:),1));
    subplot(3,5,i)
    imagesc(C)
    caxis([-500 20000])
end

%%
nMax = 50;
X = [o.ompSpotScore,o.ompSpotIntensity2];
figure
scatter(X(ThisSpots,1), X(ThisSpots,2), '.k');
grid on;
title('Click and drag to draw a rectangle');
xlabel('Spot Score')
ylabel('Spot Intensity')

% 3. Let the user draw a rectangle on the plot
% drawrectangle returns a ROIRectangle object
h = drawrectangle();

% Wait for the user to double-click or finish the rectangle (optional, but good practice)
pause(0.1);

% 4. Get the bounding box limits [xmin, ymin, width, height]
pos = h.Position;
xmin = pos(1);
ymin = pos(2);
xmax = pos(1) + pos(3);
ymax = pos(2) + pos(4);

% 5. Find the points inside the rectangle using logical indexing
in_x = (X(:,1) >= xmin) & (X(:,1) <= xmax);
in_y = (X(:,2) >= ymin) & (X(:,2) <= ymax);
inside_idx = in_x & in_y;

% Extract the coordinates of the points inside the rectangle
X_inside = X(inside_idx, :);

% 6. Highlight the selected points on the plot
hold on;
scatter(X_inside(:,1), X_inside(:,2), 50, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
legend('All Points', 'Selected Points', 'Location', 'best');
title(sprintf('Selected %d points inside the rectangle', size(X_inside, 1)));

keepSp = find(ThisSpots & inside_idx);
length(keepSp)
if length(keepSp)>nMax
    keepSp = keepSp(randperm(length(keepSp),nMax));
end
C = squeeze(nanmean(o.ompSpotColors(keepSp,:,:),1));
figure
imagesc(C)
caxis([-500 20000])

%% spot diag plot
QualityScore = get_auto_omp_score(o);
Gene = 'Wfs1';
ThisSpots = find(strcmp(o.GeneNames(o.ompSpotCodeNo),Gene));

R = randperm(length(ThisSpots), length(ThisSpots));

% Optional: Create a structure or array to store the user's choices
userResponses = cell(100, 1);

for i = 1:100
    % Display the plot
    iss_view_omp(o, 234321, [], [], ThisSpots(R(i)));

    % Prompt the user for input
    % \n adds a clean line break in the command window
    resp = input(sprintf('Spot %d/100 (Index: %d). Is this a good spot? (y/n): ', i, ThisSpots(R(i))), 's');

    % Store the response
    userResponses{i} = resp;

    % Optional: Add logic based on the response
    if strcmpi(resp, 'q')
        disp('Exiting loop early.');
        break; % Allows the user to quit early by typing 'q'
    end
end

X = [o.ompSpotScore(ThisSpots),QualityScore(ThisSpots)];
figure
scatter(X(:,1), X(:,2), '.k');
hold on
scatter(X(R(strcmp(userResponses,'n')),1), X(R(strcmp(userResponses,'n')),2), '.r');
hold on
scatter(X(R(strcmp(userResponses,'y')),1), X(R(strcmp(userResponses,'y')),2), '.g');

last = min([100 length(ThisSpots)]);
inspected_idx = ThisSpots(R(1:last));
labels = strcmp(userResponses(1:last),'y');

%%
% --- 0. Ensure base variables are computed and cast to full double ---
nCodes = length(o.CharCodes);
MaxCoef = full(double(max(abs(o.ompCoefs(:, 1:nCodes)), [], 2)));
SpotInd = sub2ind(size(o.ompCoefs), (1:length(o.ompSpotCodeNo))', double(o.ompSpotCodeNo));
SpotCoef = full(double(o.ompCoefs(SpotInd)));

if isprop(o, 'ompNeighbNearPosNeighbMultiplier') && size(o.ompNeighbNonZeros, 2) == 2
    NeighbNonZeros = double(o.ompNeighbNearPosNeighbMultiplier) * double(o.ompNeighbNonZeros(:, 1)) + ...
                     double(o.ompNeighbNonZeros(:, 2));
else
    NeighbNonZeros = double(o.ompNeighbNonZeros(:, 1));
end
NeighbNonZeros = full(NeighbNonZeros);

% --- 1. Construct X_train safely ---
f1_train = full(double(o.ompSpotScore(inspected_idx)));
f2_train = log1p(max(0, full(double(o.ompSpotIntensity2(inspected_idx)))));
f3_train = full(double(NeighbNonZeros(inspected_idx)));
f4_train = full(double(SpotCoef(inspected_idx))) ./ max(full(double(MaxCoef(inspected_idx))), eps);

% Column-align and concatenate (all are regular dense double vectors)
X_train = [f1_train(:), f2_train(:), f3_train(:), f4_train(:)];

% Ensure labels are a column vector of doubles (1 for 'y', 0 for 'n')
y_train = double(labels(:));

% --- 2. Fit the Logistic Regression Model ---
model = fitglm(X_train, y_train, 'Distribution', 'binomial', 'Link', 'logit');

% --- 3. Construct X_all safely for the entire dataset ---
f1_all = full(double(o.ompSpotScore));
f2_all = log1p(max(0, full(double(o.ompSpotIntensity2))));
f3_all = full(double(NeighbNonZeros));
f4_all = full(double(SpotCoef)) ./ max(full(double(MaxCoef)), eps);

X_all = [f1_all(:), f2_all(:), f3_all(:), f4_all(:)];

% --- 4. Predict spot quality probability (0.0 to 1.0) ---
ompProbGood = predict(model, X_all);

% Filter spots
QualOK = ompProbGood > 0.5;

figure
scatter(ompProbGood(ThisSpots),o.ompSpotScore(ThisSpots))
scatter(ompProbGood(ThisSpots),o.ompSpotIntensity2(ThisSpots))
%%
I = imadjust(imread(fullfile(o.OutputDirectory,'anchor_image.tif'))); % background image = Anchor
Roi = round([1, max(o.ompSpotGlobalYX(:,2)), ...
1, max(o.ompSpotGlobalYX(:,1))]);
o.plot(I,Roi,'OMP',QualOK);
daspect([1 1 1])

%%

% 1. Extract and Clean Data
nSpots = length(o.ompSpotCodeNo);
scores = full(double(o.ompSpotScore(:)));
intensities = full(double(o.ompSpotIntensity2(:)));
geneIDs = full(double(o.ompSpotCodeNo(:)));

% Avoid negative or invalid values for intensity
cleanIntensities = max(0, intensities);

% 2. Divide SpotScore & SpotIntensity2 into 4 Quantile Bins (16 grid cells)
nBins = 4;
qScore = quantile(scores, linspace(0, 1, nBins + 1));
qScore(1) = -Inf; qScore(end) = Inf;

qInt = quantile(cleanIntensities(cleanIntensities > 0), linspace(0, 1, nBins + 1));
qInt = unique(qInt); % guard against identical quantiles
qInt(1) = -Inf; qInt(end) = Inf;

binScore = discretize(scores, qScore);
binInt   = discretize(cleanIntensities, qInt);

% 3. Stratified Sampling Across Bins & Genes
% Target ~100 spots total (~6 spots per grid bin)
spotsPerBin = 15;
sampledIndices = [];

for s_bin = 1:nBins
    for i_bin = 1:(length(qInt) - 1)
        % Spots in this specific (Score, Intensity) bin
        inCell = find(binScore == s_bin & binInt == i_bin);
        
        if isempty(inCell)
            continue;
        end
        
        % Ensure gene diversity: sort candidates by unique genes first
        [~, uniqueGeneOrder] = sort(rand(length(inCell), 1));
        inCellShuffled = inCell(uniqueGeneOrder);
        [~, firstOfGene] = unique(geneIDs(inCellShuffled), 'stable');
        
        % Pick spots from distinct genes first, fill remainder if needed
        distinctGeneSpots = inCellShuffled(firstOfGene);
        if length(distinctGeneSpots) >= spotsPerBin
            picked = distinctGeneSpots(1:spotsPerBin);
        else
            remaining = setdiff(inCellShuffled, distinctGeneSpots, 'stable');
            picked = [distinctGeneSpots; remaining(1:min(end, spotsPerBin - length(distinctGeneSpots)))];
        end
        
        sampledIndices = [sampledIndices; picked]; %#ok<AGROW>
    end
end

% Randomize the order of presentation so you don't review in sorted order
sampledIndices = sampledIndices(randperm(length(sampledIndices)));
totalToReview = length(sampledIndices);

fprintf('Selected %d representative spots spanning %d genes and score/intensity bins.\n', ...
    totalToReview, length(unique(geneIDs(sampledIndices))));

%% 4. Interactive Review Loop
% Controls: 
%   'y' -> Good spot (1)
%   'n' -> Bad spot (0)
%   's' -> Skip (ambiguous)
%   'q' -> Quit early and save what you have

reviewResults = table('Size', [totalToReview, 5], ...
    'VariableTypes', {'double', 'string', 'double', 'double', 'double'}, ...
    'VariableNames', {'SpotIndex', 'Gene', 'Score', 'Intensity', 'Label'});

for i = 1:totalToReview
    idx = sampledIndices(i);
    geneName = string(o.GeneNames{geneIDs(idx)});
    spScore  = scores(idx);
    spInt    = intensities(idx);
    
    % Display plot
    iss_view_omp(o, 234321, [], [], idx);
    
    % Prompt user
    prompt = sprintf('[%d/%d] Gene: %s | Score: %.3f | Intensity: %.2f | Good? (y/n, s=skip, q=quit): ', ...
        i, totalToReview, geneName, spScore, spInt);
    resp = lower(strtrim(input(prompt, 's')));
    
    if strcmp(resp, 'q')
        fprintf('Exiting early at spot %d.\n', i);
        reviewResults = reviewResults(1:i-1, :);
        break;
    elseif strcmp(resp, 'y')
        label = 1;
    elseif strcmp(resp, 'n')
        label = 0;
    else
        label = NaN; % skipped
    end
    
    reviewResults(i, :) = {idx, geneName, spScore, spInt, label};
end

% Filter out skipped spots
validLabels = reviewResults(~isnan(reviewResults.Label), :);
fprintf('\nDone! Collected %d labels (%d Good, %d Bad).\n', ...
height(validLabels), sum(validLabels.Label == 1), sum(validLabels.Label == 0));

inspected_idx = validLabels.SpotIndex;
labels = validLabels.Label;
%%
% --- 0. Ensure base variables are computed and cast to full double ---
nCodes = length(o.CharCodes);
MaxCoef = full(double(max(abs(o.ompCoefs(:, 1:nCodes)), [], 2)));
SpotInd = sub2ind(size(o.ompCoefs), (1:length(o.ompSpotCodeNo))', double(o.ompSpotCodeNo));
SpotCoef = full(double(o.ompCoefs(SpotInd)));

if isprop(o, 'ompNeighbNearPosNeighbMultiplier') && size(o.ompNeighbNonZeros, 2) == 2
    NeighbNonZeros = double(o.ompNeighbNearPosNeighbMultiplier) * double(o.ompNeighbNonZeros(:, 1)) + ...
                     double(o.ompNeighbNonZeros(:, 2));
else
    NeighbNonZeros = double(o.ompNeighbNonZeros(:, 1));
end
NeighbNonZeros = full(NeighbNonZeros);

% --- 1. Construct X_train safely ---
f1_train = full(double(o.ompSpotScore(inspected_idx)));
f2_train = log1p(max(0, full(double(o.ompSpotIntensity2(inspected_idx)))));
f3_train = full(double(NeighbNonZeros(inspected_idx)));
f4_train = full(double(SpotCoef(inspected_idx))) ./ max(full(double(MaxCoef(inspected_idx))), eps);

% Column-align and concatenate (all are regular dense double vectors)
X_train = [f1_train(:), f2_train(:), f3_train(:), f4_train(:)];

% Ensure labels are a column vector of doubles (1 for 'y', 0 for 'n')
y_train = double(labels(:));

% --- 2. Fit the Logistic Regression Model ---
model = fitglm(X_train, y_train, 'Distribution', 'binomial', 'Link', 'logit');

% --- 3. Construct X_all safely for the entire dataset ---
f1_all = full(double(o.ompSpotScore));
f2_all = log1p(max(0, full(double(o.ompSpotIntensity2))));
f3_all = full(double(NeighbNonZeros));
f4_all = full(double(SpotCoef)) ./ max(full(double(MaxCoef)), eps);

X_all = [f1_all(:), f2_all(:), f3_all(:), f4_all(:)];

% --- 4. Predict spot quality probability (0.0 to 1.0) ---
ompProbGood = predict(model, X_all);

% Filter spots
QualOK = ompProbGood > 0.5;

figure
scatter(ompProbGood(ThisSpots),o.ompSpotScore(ThisSpots))
figure
scatter(ompProbGood(ThisSpots),o.ompSpotIntensity2(ThisSpots))
%%
I = imadjust(imread(fullfile(o.OutputDirectory,'anchor_image.tif'))); % background image = Anchor
Roi = round([1, max(o.ompSpotGlobalYX(:,2)), ...
1, max(o.ompSpotGlobalYX(:,1))]);
o.plot(I,Roi,'OMP',QualOK);
daspect([1 1 1])

% Strip internal training cache to keep it lightweight (~a few KB)
compactModel = compact(model);

% Also save the feature names / beta coefficients directly
beta_coeffs = model.Coefficients.Estimate;

% Save to a .mat file
save('D:\ISS\omp_logistic_model.mat', 'compactModel', 'beta_coeffs');
fprintf('Model saved successfully to omp_logistic_model.mat\n');

%%
clear
load('D:\ISS\Round_-\outputB\1\oCall_spots_OMP.mat')
[probGood, QualOK] = predict_omp_quality(o, 'D:\ISS\omp_logistic_model.mat');
I = imadjust(imread(fullfile('Z:\shared\Stephane\Thomas\humansection-20260901\output','anchor_image.tif'))); % background image = Anchor
Roi = round([1, max(o.ompSpotGlobalYX(:,2)), ...
1, max(o.ompSpotGlobalYX(:,1))]);
o.plot(I,Roi,'OMP',QualOK);
daspect([1 1 1])

