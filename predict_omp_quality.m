function [probGood, QualOK] = predict_omp_quality(o, modelPathOrObj, threshold, pf)
% PREDICT_OMP_QUALITY Applies a saved logistic model to an ISS OMP dataset.
%
% Usage:
%   [probGood, QualOK] = predict_omp_quality(o, 'omp_logistic_model.mat');
%   [probGood, QualOK] = predict_omp_quality(o, 'omp_logistic_model.mat', 0.5);

    if nargin < 3 || isempty(threshold)
        threshold = 0.5; % Default probability cutoff
    end
    if nargin < 4 || isempty(pf)
        pf = 'omp';
    end

    % --- 1. Load the Model ---
    if ischar(modelPathOrObj) || isstring(modelPathOrObj)
        loaded = load(modelPathOrObj);
        if isfield(loaded, 'compactModel')
            model = loaded.compactModel;
        else
            error('The .mat file must contain "compactModel".');
        end
    else
        model = modelPathOrObj;
    end

    % --- 2. Extract Features on New Dataset ---
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

    % Feature 1: SpotScore
    f1 = full(double(o.([pf, 'SpotScore'])));
    
    % Feature 2: Intensity (Log-transformed)
    rawInt = full(double(o.([pf, 'SpotIntensity2'])));
    f2 = log1p(max(0, rawInt));
    
    % Feature 3: Neighborhood count
    f3 = full(double(NeighbNonZeros));
    
    % Feature 4: Coefficient Uniqueness Ratio
    f4 = full(double(SpotCoef)) ./ max(full(double(MaxCoef)), eps);

    % Assemble matrix (all dense column vectors)
    X_new = [f1(:), f2(:), f3(:), f4(:)];

    % --- 3. Predict Probabilities ---
    probGood = predict(model, X_new);
    QualOK = probGood >= threshold;
end