function o = get_initial_bled_codes(o)
% o = o.get_initial_bled_codes
% This finds o.z_scoreBleedMatrix(o.nBP,o.nBP) using the initial
% dpSpotColors from all rounds. It then uses this to produce
% iompBledCodes(nGenes,o.nBP*o.nRounds). 
%
% Only difference between this bleed matrix and o.BleedMatrix is that here
% we force o.BleedMatrixType = 'Single'. 
% 
% Using o.UseChannels and o.UseRounds, you can do spot calling
% without using certain rounds and colour channels.
%
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end
%%

%Only using channels and rounds given by o.UseChannels and o.UseRounds
if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end
    
if isempty(o.UseRounds)
    o.UseRounds = 1:o.nRounds;
end

nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);
%o.cSpotColors = o.cSpotColors(:,o.UseChannels,o.UseRounds);

%FILTER OUT REALLY HIGH VALUES
%Good = all(o.cSpotColors(:,:)<10000,2);         
%o.cSpotColors = o.cSpotColors(Good,:,:);
%o.cSpotIsolated = o.cSpotIsolated(Good);
%o.SpotGlobalYX = o.SpotGlobalYX(Good,:);

% First divide by o.ExtractScale so comparable between experiments

%Normalise each colour channel by a percentile as to correct for weaker
%colour channels
p = zeros(1,o.nBP,o.nRounds);
for b = 1:o.nBP
    bSpotColors = o.dpSpotColors(:,b,:);
    p(:,b,:) = prctile(bSpotColors(:), o.SpotNormPrctile);
end
o.z_scoreSCALE = p;
o.z_scoreSHIFT = zeros(size(o.z_scoreSCALE));

if strcmpi(o.ompBleedMatrixEigMethod,'Mean')
    SpotIsolated = o.dpSpotIsolated;
elseif strcmpi(o.ompBleedMatrixEigMethod,'Median')
    %Use all spots if taking median as robust to outliers.
    SpotIsolated = true(size(o.dpSpotIsolated));
end
DiagMeasure = 0;
nTries = 1;
nIter = 100;
BleedMatrixScoreThresh = o.BleedMatrixScoreThresh;

while DiagMeasure<nChans && nTries<nIter
    SpotColors = bsxfun(@rdivide, o.dpSpotColors, p);
    [BleedMatrix,DiagMeasure] = ...
        get_bleed_matrix(o,SpotColors,SpotIsolated,BleedMatrixScoreThresh,nTries);
    %If bleed matrix not diagonal, try only using spots that are more like
    %the initial diagonal bleed matrix. 
    if DiagMeasure<nChans
        BleedMatrixScoreThresh = BleedMatrixScoreThresh + o.BleedMatrixScoreThreshStep;
        if BleedMatrixScoreThresh>o.BleedMatrixScoreThreshMax
            if sum(SpotIsolated)~=size(SpotIsolated,1)
                SpotIsolated = true(size(o.dpSpotIsolated));
                BleedMatrixScoreThresh = 0;
                warning('Bleed matrix not diagonal - Now using all spots, not just isolated');
            elseif strcmpi(o.ompBleedMatrixEigMethod,'Median')
                SpotIsolated = o.dpSpotIsolated;
                BleedMatrixScoreThresh = 0;
                o.ompBleedMatrixEigMethod = 'Mean';
                warning('Bleed matrix not diagonal - Setting o.ompBleedMatrixEigMethod = Mean');
            else
                break;
            end
        end       
        warning('Bleed matrix not diagonal - Setting BleedMatrixScoreThresh = %.1f',...
            BleedMatrixScoreThresh);
    elseif DiagMeasure>=nChans && nTries>1
        fprintf('Bleed matrix now diagonal\n');
    end
    nTries = nTries+1;
end
if DiagMeasure<nChans
    error('Bleed matrix not diagonal')
end
o.z_scoreBleedMatrix = BleedMatrix;

% now load in the code book and apply bleeds to it
CharCode = o.CharCodes;
nCodes = length(CharCode);

% create numerical code (e.g. 33244 for CCGAA)
NumericalCode = zeros(nCodes, o.nRounds);
for r=1:o.nRounds
    if r<=o.nRounds-o.nRedundantRounds
        for c=1:nCodes
            try
                [~, NumericalCode(c,r)] = ismember(CharCode{c}(r), o.bpLabels);
            catch
                error('Code %s has no channel for round %.0f.\nCheck for missing leading zeros in CodeFile:\n%s.',GeneName{c},r,o.CodeFile);
            end
        end
    else
        % redundant round - compute codes automatically
        % find pseudobases for this code
        for c=1:nCodes
            PseudoCode = repmat('0',1,o.nRounds-o.nRedundantRounds);
            for p = 1:length(o.RedundantPseudobases)
                PseudoCode(1,ismember(CharCode{c}, o.RedundantPseudobases{p}))=('0'+p);
            end
            % now match them to the redundant codes
            for cc=1:nChans
                rrn = r-o.nRounds+o.nRedundantRounds;
                if ~isempty(regexp(PseudoCode, o.RedundantCodes{rrn,cc}, 'once'))
                    NumericalCode(c,r)=cc;
                end
            end
        end
    end
end

BledCodes = nan(nCodes, o.nBP,o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:nRounds
        ChannelNo = NumericalCode(i,o.UseRounds(r));
        if any(o.UseChannels == ChannelNo) == 0 continue; end        
        BledCodes(i,o.UseChannels,o.UseRounds(r)) = BleedMatrix(o.UseChannels,ChannelNo);
    end
end
o.iompBledCodes = BledCodes(:,:);
