function o = reduce_memory_omp(o, nNeighbThresh, RemovePixelInfo)
%% o = o.reduce_memory_omp(nNeighbThresh, RemovePixelInfo);
% This reduces the memory of the o object by only keeping ompSpots for
% which o.ompNeighbNonZeros > nNeighbThresh.
% The quality thresholding requires o.ompNeighbNonZeros>o.ompNeighbThresh2
% for any spot to be accepted so must have nNeighbThresh below this.
% If RemovePixelInfo is true then will also remove pixel based data.
% This is worth doing if pixel based data is saved in a previous version.

if nargin<2 || isempty(nNeighbThresh)
    nNeighbThresh = o.ompNeighbThresh2-2;
end
if nargin<3 || isempty(RemovePixelInfo)
    RemovePixelInfo = false;
end

if nNeighbThresh>o.ompNeighbThresh2
    error('Keep nNeighbThresh below o.ompNeighbThresh2 = %.0f',o.ompNeighbThresh2);
end

o.ompCoefs = sparse(o.ompCoefs);  % Mainly zeros so use sparse

if nNeighbThresh>0
    ompKeep = o.ompNeighbNonZeros > nNeighbThresh;
    o.ompSpotColors = o.ompSpotColors(ompKeep,:,:);
    o.ompSpotCodeNo = o.ompSpotCodeNo(ompKeep);
    o.ompSpotGlobalYX = o.ompSpotGlobalYX(ompKeep,:);
    o.ompLocalTile = o.ompLocalTile(ompKeep);
    o.ompCoefs = o.ompCoefs(ompKeep,:);
    o.ompNeighbNonZeros = o.ompNeighbNonZeros(ompKeep);
    o.ompSpotIntensity2 = o.ompSpotIntensity2(ompKeep);
    if isempty(o.ompSpotScore)
        o.ompSpotScore = o.get_omp_score;
    else
        o.ompSpotScore = o.ompSpotScore(ompKeep);
    end
end

if RemovePixelInfo
    % Pixel Info
    o.pxSpotColors = [];
    o.pxSpotCodeNo = [];
    o.pxSpotGlobalYX = [];
    o.pxLogProbOverBackground = [];
    o.pxSpotScore = [];
    o.pxSpotScoreDev = [];
    o.pxSpotIntensity = [];
    o.pxLocalTile = [];
    o.pxSpotBestGene = [];
    o.pxSpotScore2 = [];
    % Prob Info
    o.pSpotColors = [];
    o.pSpotCodeNo = [];
    o.pSpotGlobalYX = [];
    o.pLogProbOverBackground = [];
    o.pSpotScore = [];
    o.pSpotScoreDev = [];
    o.pSpotIntensity = [];
    o.pSpotIntensity2 = [];
    o.pLocalTile = [];
    o.pSpotIsolated = [];
end
end

