function SpotNo = iss_view_spot(o, FigNo, ImSz, SpotLocation, ScoreMethod,...
    IncludeGT, Filter, Norm, SpotNum)
%% SpotNo = iss_view_spot(o, FigNo, ImSz, SpotLocation, ScoreMethod,...
% IncludeGT, Filter, Norm, SpotNum)
%
% Check PCR by plotting location of spot in each round and color channel
%
% o: iss object.
% FigNo: o.plot figure number (default, current figure)
% ImSz: radius of image that is plotted for each round and channel.
%   Default value is 7 pixels.
% SpotLocation: logical,  if true, will use location of spot closest to
%   crosshair, otherwise will use actual position of crosshair. Default is false.
% ScoreMethod: The set of spots to consider e.g. ScoreMethod = 'OMP'
%   would find spot from the set o.ompSpotGlobalYX. Defaults to value used
%   in the plot image.
% IncludeGT: if true, will also plot the ground truth rounds.
% Filter: true to get the SpotColors from the filtered images i.e. the
%   tiles in o.TileDirectory. false to get the raw images i.e. directly
%   from the nd2 file (and then focus-stack).
% Norm: true to normalise spot colors by round/channel. May be different
%   for ScoreMethod = 'OMP' or other. false to see raw spot colors.
% SpotNum: index of spot that you want to look at.
% SpotNo: returns the index of the spot analyzed.

%%
if nargin<3 || isempty(ImSz)
    ImSz = 7;
end
if ImSz>100
    warning('ImSz too large, setting to 7');
    ImSz = 7;
end

if nargin<6 || isempty(IncludeGT)
    IncludeGT = false;
end

if nargin<7 || isempty(Filter)
    Filter = true;
end

if nargin<8 || isempty(Norm)
    Norm = false;
end

if ~Filter
    Norm = false;
end

if nargin<4 || isempty(SpotLocation)
    SpotLocation = false;
end

if nargin<5
    ScoreMethod = [];
end

if nargin<9
    SpotNum = [];
end

[xy, SpotLocation, ScoreMethod, SpotNo, Dist]  = ...
    get_crosshair_location(o, FigNo, SpotLocation, ScoreMethod, SpotNum);
pf = o.CallMethodPrefix(ScoreMethod);
SpotCodeNo = o.([pf,'SpotCodeNo'])(SpotNo);


[SpotColors, PointCorrectedLocalYX] = get_spot_colors_grid(o, pf, xy, ImSz, SpotNo,...
    SpotLocation, IncludeGT, Filter);
if IncludeGT
    Clim = zeros(2,o.nBP,max([o.nRounds,o.gtRounds]));
else
    Clim = zeros(2,o.nBP,o.nRounds);
end
climThresh = 100;
Clim(1,:,:) = min(min(SpotColors),-climThresh);
Clim(2,:,:) = max(max(SpotColors),climThresh);

if strcmpi(ScoreMethod,'OMP')
    NormFactor = o.z_scoreSCALE;
else
    NormFactor = o.BledCodesPercentile;
end

try
    clf(27642)
    figure(27642)
catch
    figure(27642)
end
set(gcf,'Position',[164,108,1621,805]);

if Norm
    Clim(:,:,1:o.nRounds) = Clim(:,:,1:o.nRounds)./NormFactor;
    Clim(1,:,:) = min(min(min(Clim(:,:,1:o.nRounds))));
    Clim(2,:,:) = max(max(max(Clim(:,:,1:o.nRounds))));
    SpotColors(:,:,1:o.nRounds) = double(SpotColors(:,:,1:o.nRounds))./NormFactor;
end
plot_spot_colors_grid(o, SpotColors, PointCorrectedLocalYX,...
    ImSz, Dist, SpotCodeNo, Clim, IncludeGT, Filter);

iss_view_plot_title(o, ScoreMethod, SpotLocation, SpotNo);

end



