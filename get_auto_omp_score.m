function QualityScore = get_auto_omp_score(o, pf)
    if nargin < 2; pf = 'omp'; end

    % --- 1. Gene Uniqueness Ratio (U in [0, 1]) ---
    nCodes = length(o.CharCodes);
    MaxCoef = double(max(abs(o.ompCoefs(:, 1:nCodes)), [], 2));
    SpotInd = sub2ind(size(o.ompCoefs), (1:length(o.ompSpotCodeNo))', double(o.ompSpotCodeNo));
    SpotCoef = double(o.ompCoefs(SpotInd));
    
    % Ensure dense double and clip to [0, 1]
    U = full(max(0, min(1, SpotCoef ./ max(MaxCoef, eps))));

    % --- 2. Spectral Match (S) ---
    % SpotScore is cosine similarity (clip negatives to 0)
    S = full(max(0, double(o.([pf, 'SpotScore']))));

    % --- 3. Self-Calibrated Intensity / SNR (I) ---
    rawIntensity = full(double(o.([pf, 'SpotIntensity2'])));
    
    % Clamp negative intensities to 0 so log1p never produces complex numbers
    cleanIntensity = max(0, rawIntensity);
    
    % Calculate noise floor from positive intensities
    posIntensities = cleanIntensity(cleanIntensity > 0);
    if isempty(posIntensities)
        noiseFloor = 1.0;
    else
        noiseFloor = median(posIntensities);
    end
    
    SNR = cleanIntensity ./ max(noiseFloor, eps);
    I_score = log1p(SNR); 

    % --- 4. Spatial Neighborhood Support (N in [0, 1]) ---
    if isprop(o, 'ompNeighbNearPosNeighbMultiplier') && size(o.ompNeighbNonZeros, 2) == 2
        Neighb = double(o.ompNeighbNearPosNeighbMultiplier) * double(o.ompNeighbNonZeros(:, 1)) + ...
                 double(o.ompNeighbNonZeros(:, 2));
    else
        Neighb = double(o.ompNeighbNonZeros(:, 1));
    end
    Neighb = full(Neighb);

    % Safe exponential saturating curve
    N_score = 1 - exp(-Neighb ./ 3);

    % --- Final Clean Dense Double Vector ---
    QualityScore = full(real(S .* I_score .* U .* N_score));
end