function [xy, SpotLocation, ScoreMethod, SpotNo, Dist]  = ...
    get_crosshair_location(o, FigNo, SpotLocation, ScoreMethod, SpotNum)
%% [xy, SpotLocation, ScoreMethod, SpotNo, Dist]  = ...
%    get_crosshair_location(o, FigNo, SpotLocation, ScoreMethod, SpotNum)
% Functions to view specific spots e.g. iss_view_spot, require a spot to be
% found using a cross-hair on the o.plot image. This function returns the
% spot and other information found.
% Input
%   o: iss object.
%   FigNo: o.plot figure number
%   SpotLocation: true to get nearest spot to crosshair, false returns
%       exact position of crosshair.
%   ScoreMethod: The set of spots to consider e.g. ScoreMethod = 'OMP'
%       would find spot from the set o.ompSpotGlobalYX. Defaults to value used
%       in the plot image.
%   SpotNum: if know spot index already, can specify and cross-hair won't
%       be required.
% Output
%   xy: location of cross-hair or spot
%   SpotLocation: same as input but might be set to true if click exactly
%   on spot.
%   SpotNo: index of nearest spot to cross-hair.
%   Dist: distance of nearest spot to cross-hair, set to 0 if SpotLocation
%   = true.

%%
if isempty(ScoreMethod)
    try
        S = evalin('base', 'issPlot2DObject');
        ScoreMethod = S.CallMethod;
    catch
        ScoreMethod = [];
    end
end

if ~isempty(SpotNum)
    if length(SpotNum)==2
        SpotLocation = false;
        xy = [SpotNum(2),SpotNum(1)];
        S.SpotYX = o.([o.CallMethodPrefix(ScoreMethod),'SpotGlobalYX']);
        [Dist,SpotNo] = min(sum(abs(S.SpotYX-[xy(2),xy(1)]),2));
        if round(Dist)==0
            SpotLocation=true;
        end
    else
        SpotLocation = true;
        SpotNo = SpotNum;
        xy = o.([o.CallMethodPrefix(ScoreMethod),'SpotGlobalYX'])(SpotNo,[2,1]);
        Dist = 0;
    end
else
    figure(FigNo);
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    try
        S = evalin('base', 'issPlot2DObject');
        if ~strcmpi(S.CallMethod,ScoreMethod)
            S.QualOK = 1;
        end
    catch
        S.QualOK = 1;
        S.Roi = [1,inf,1,inf];
    end
    S.SpotYX = o.([o.CallMethodPrefix(ScoreMethod),'SpotGlobalYX']);
    if size(S.SpotYX,1)~=size(S.QualOK,1)
        S.QualOK = 1;
    end
    %Only consider spots that can be seen in current plot
    InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
    PlotSpots = find(InRoi & S.QualOK);        
    [Dist,SpotIdx] = min(sum(abs(S.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
    SpotNo = PlotSpots(SpotIdx);
    if SpotLocation || round(Dist)==0
        SpotLocation = true;
        xy = S.SpotYX(SpotNo,[2,1]);
        Dist = 0;
    end   
end

if ~ismember({ScoreMethod},o.CallMethods)
    error('Method invalid, must be member of o.CallMethods.');
end

end

