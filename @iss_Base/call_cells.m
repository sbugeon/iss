function [o,cellCallData] = call_cells(o, gSet, Boundaries, CellMap, Method)
% o = o.call_cells(gSet)
%
% Cell calling via negative binomial model
%
% input gSet: GeneSet structure containing results of scRNAseq clustering
%
% Creates outputs:
% pCellClass: posterior probability of each cell to be in each class (nCells*nClasses)
% pSpotCell: posterior probability of each spot to be in top 5 neighboring cells (nSpots * nCells, sparse)
% note that the last class is a zero-expressing cell; the last cell is background
%%
if nargin<2 || isempty(Method)
    Method = 'OMP';
    warning('Method not provided, using OMP.');
end
if ~ismember({Method},o.CallMethods)
    error('Method invalid, must be member of o.CallMethods.');
end
%
y0 = min(o.CellCallRegionYX(:,1));
x0 = min(o.CellCallRegionYX(:,2));

pf = o.CallMethodPrefix(Method);
%% load properties of local region of interest
if nargin<3 | isempty(CellMap)
    load(o.CellMapFile); % CellMap and y0, y1, x0, x1 that are its coords in full image
    % get info about cells
    rp = regionprops(CellMap);
    CellYX = fliplr(vertcat(rp.Centroid)) + [y0 x0]; % convert XY to YX
    CellArea0 = vertcat(rp.Area);
else
    CellYX = zeros(length(Boundaries),2);
    CellArea0 = zeros(length(Boundaries),1);
    
    for iB = 1 : length(Boundaries)
        if ~ isempty(Boundaries{iB}) & ~iscell(Boundaries{iB})
            Centroids = mean(Boundaries{iB},1);
            CellYX(iB,:) = fliplr(Centroids) ; % convert XY to YX
            CellArea0(iB,:) = polyarea(Boundaries{iB}(:,1),Boundaries{iB}(:,2));
        end
    end
end
MeanCellRadius = mean(sqrt(CellArea0/pi))*.5; % the dapi part is only half of the typical radius
RelCellRadius = [sqrt(CellArea0/pi)/MeanCellRadius; 1]; % but here we want the whole thing
%% diagnostic parameters in local coordinates XY not YX!!!
o.CellCallShowCenter = [18836 21762];
o.CellCallShowRad = 100;
o.ExampleCellCenter = o.CellCallShowCenter;

%% include correct spots
AllGeneNames = o.GeneNames(o.([pf,'SpotCodeNo']));
IncludeSpot = ~ismember(AllGeneNames, o.ExcludeGenes) ...
    & inpolygon(o.([pf,'SpotGlobalYX'])(:,1), o.([pf,'SpotGlobalYX'])(:,2), o.CellCallRegionYX(:,1), o.CellCallRegionYX(:,2)) ...
    & o.quality_threshold(Method);
% SpotYX is only the spots we are bothered with, in global coordinates
SpotYX = round(o.([pf,'SpotGlobalYX'])(IncludeSpot,:));
SpotGeneName = AllGeneNames(IncludeSpot);
%% get arrays ready

% SpotGene(nS): which gene is each spot
% MeanClassExp(nK,nG): mean expression of each gene in each class
% Neighbors(nS, nN): closest neighboring cells for each spot
% D(nS, nN): distance penalty for each of these
% GeneNames(nG): name of each gene
% ClassNames(nK): name of each class

[GeneNames, ~, SpotGeneNo] = unique(SpotGeneName);
TotGeneSpots = accumarray(SpotGeneNo,1);
ClassNames = vertcat(unique(gSet.Class, 'stable'), {'Zero'});

nG = length(GeneNames);
nK = length(ClassNames); % last is zero-expression
nC = size(CellYX,1)+1; % last is misreads
nS = size(SpotYX,1);
nN = o.nNeighbors+1; % last is misreads (always a neighbor)

ClassPrior = [.5*ones(1,nK-1)/nK .5];

ClassDisplayNames = ClassNames;

MeanClassExp = zeros(nK, nG);
gSub = gSet.GeneSubset(GeneNames);
for k=1:nK-1 % don't include last since it is zero-expression class
    MeanClassExp(k,:) = o.Inefficiency * mean(gSub.ScaleCell(0).CellSubset(ClassNames{k}).GeneExp,2)';
end
lMeanClassExp = log(MeanClassExp + o.SpotReg);

% now find each spot's neighboring cells and distances (nS, nN)
[Neighbors, Dist] = knnsearch(CellYX, SpotYX, 'K', nN);
Neighbors(:,end) = nC; % set last neighbor to misreads

D = -Dist.^2./(2*MeanCellRadius^2) - log(2*pi*MeanCellRadius^2); % don't normalize: bigger cells express more
D(:,end) = log(o.MisreadDensity); % this is log likelihood of misread

% any inside cell radius given a bonus
SpotInCell = IndexArrayNan(CellMap, (SpotYX - [y0 x0])');
if Neighbors(SpotInCell>0,1)~=SpotInCell(SpotInCell>0)
    error('a spot is in a cell not closest neighbor!');
end
D(SpotInCell>0, 1) = D(SpotInCell>0, 1) + o.InsideCellBonus;

LogClassPrior = log(ClassPrior);

% this is area factor relative to that of the average cell
CellAreaFactor = (exp(-RelCellRadius.^2/2)*(1-exp(o.InsideCellBonus)) + exp(o.InsideCellBonus)) ...
    / (exp(-1/2)*(1-exp(o.InsideCellBonus)) + exp(o.InsideCellBonus));
%% initialize variables for main loop

pSpotNeighb = zeros(nS, nN); % prob each spot goes to each neighboring cell: last assigned to noise
pCellClass = zeros(nC, nK); % prob each cell goes to each class: last has zero expression

% start a spot in cell it is in, otherwise misread
pSpotNeighb(Neighbors==SpotInCell)=1;
pSpotNeighb(SpotInCell==0,end)=1;

% gammas start off as priors
eSpotGamma = ones(nC, nK, nG);
elSpotGamma = ones(nC, nK, nG)*psi(1); % start with r=1 prior, no spots

eGeneGamma = ones(nG,1); % start with just 1

% this is to check convergence
pSpotNeighbOld = zeros(nS, nN);

%% now main loop
for i=1:o.CellCallMaxIter
    % CellGeneCount(nC, nG): number of copies of each gene in each cell
    CellGeneCount = zeros(nC,nG);
    for n=1:nN-1
        c = Neighbors(:,n);
        CellGeneCount = CellGeneCount + accumarray([c, SpotGeneNo], pSpotNeighb(:,n), [nC,nG]);
    end
    
    %% call cell gammas
    
    % eSpotGamma(nC, nK, nG); expected gamma parameter
    % elSpotGamma(nC, nK, nG); expected log gamma parameter
    ScaledMean = CellAreaFactor.*reshape(MeanClassExp,[1 nK nG]);
    eSpotGamma = (o.rSpot+reshape(CellGeneCount,[nC 1 nG]))./(o.rSpot + ScaledMean);
    elSpotGamma = psi(o.rSpot+reshape(CellGeneCount,[nC 1 nG])) - log(o.rSpot + ScaledMean); % expectation of log gamma
    %% call cells
    
    % ScaledExp(nC, nK, nG): expected expression under current parameters
    ScaledExp = reshape(MeanClassExp,[1 nK nG]) .* reshape(eGeneGamma,[1 1 nG]) .* CellAreaFactor   + o.SpotReg;
    
    % pNegBin(nC, nK, nG): negbin parameter
    pNegBin = ScaledExp ./ (o.rSpot + ScaledExp);
    
    % wCellClass(nC, nK): summed log likelihoods
    wCellClass = sum(reshape(CellGeneCount,[nC 1 nG]).*log(pNegBin) + o.rSpot*log(1-pNegBin),3) + LogClassPrior;
    
    % pCellClass(nC, nK): probabilities
    pCellClass = LogLtoP(wCellClass')';
    
    %% call spots
    % wSpotCell(nS, nN)
    aSpotCell = zeros(nS, nN);
    for n=1:nN-1 % don't include misread possibility
        c = Neighbors(:,n);
        aSpotCell(:,n) = sum(pCellClass(c,:) .* lMeanClassExp(:,SpotGeneNo)',2) + ...
            sum(pCellClass(c,:) .* bi(elSpotGamma, c, 1:nK, SpotGeneNo), 2);
    end
    wSpotCell = aSpotCell + D ;
    
    pSpotNeighb = LogLtoP(wSpotCell')';
    MeanProbChanged = max(abs(pSpotNeighb(:)-pSpotNeighbOld(:)));
    fprintf('Iteration %d, mean prob change %f\n', i, MeanProbChanged)
    Converged = ( MeanProbChanged<o.CellCallTolerance);
    pSpotNeighbOld = pSpotNeighb;
    
    %% call gene gammas (etas)
    % to count non-background expression of each gene first compute background
    TotPredictedB = accumarray(SpotGeneNo, pSpotNeighb(:,end), [nG 1]);
    % and total spots in zero cells:
    pCellZero = pCellClass(:,nK); % prob a cell is class zero (nC)
    pSpotZero = sum(pSpotNeighb(:,1:nN-1).*pCellZero(Neighbors(:,1:nN-1)),2); % prob a spot comes from cell class zero (nS)
    TotPredictedZ = accumarray(SpotGeneNo, pSpotZero);
    
    % total counts predicted by all cells of each class (nK, nG)
    ClassTotPredicted = shiftdim(sum(eSpotGamma.*pCellClass.*CellAreaFactor,1),1).*(MeanClassExp + o.SpotReg);
    % total of each gene (nG):
    TotPredicted = sum(ClassTotPredicted(1:nK-1,:),1)';
    
    eGeneGamma = (o.rGene + TotGeneSpots - TotPredictedB - TotPredictedZ)./(o.rGene + TotPredicted);
    if 0
        for gg=1:nG
            fprintf('%s:\t%f\n', GeneNames{gg}, eGeneGamma(gg));
        end
    end
    
    
    if Converged; break; end
end

% save(o.CellMapFile, 'ScaledExp', 'eGeneGamma', 'IncludeSpot', '-append');

%% make dense array output

o.pSpotCell = sparse(repmat(1:nS,1,nN)', Neighbors(:), pSpotNeighb(:));
o.CellYX = CellYX;
o.pCellClass = pCellClass;
o.ClassNames = ClassNames;

cellCallData.GeneNames = GeneNames;
cellCallData.ClassNames = ClassNames;
cellCallData.IncludeSpot = IncludeSpot;

cellCallData.NegBin = sum(reshape(CellGeneCount,[nC 1 nG]).*log(pNegBin) + o.rSpot*log(1-pNegBin),3) - LogClassPrior;
cellCallData.pCellClass = pCellClass;
cellCallData.CellGeneCount = CellGeneCount;
cellCallData.Neighbors = Neighbors;
cellCallData.pSpotNeighb = pSpotNeighb;
cellCallData.pSpotCell = o.pSpotCell;
cellCallData.CellYX = o.CellYX;
cellCallData.CellMap = CellMap;
cellCallData.CellGeneCount = CellGeneCount;
cellCallData.eGeneGamma = eGeneGamma;
cellCallData.ScaledExp = ScaledExp;
cellCallData.ClassTotPredicted = ClassTotPredicted;
cellCallData.GeneContrib = sum(reshape(CellGeneCount,[nC 1 nG]).*log(pNegBin) + o.rSpot*log(1-pNegBin),3);
end