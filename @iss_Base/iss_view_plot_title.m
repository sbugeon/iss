function iss_view_plot_title(o, ScoreMethod, SpotLocation, SpotNo)
%% iss_view_plot_title(o, ScoreMethod, SpotLocation, SpotNo)
% This just adds a title to some of the iss_view_spot type plots 
% (where you use cross-hair to select location in o.plot) indicating
% the spot number and some of the variables associated with it.
% o: iss object.
% ScoreMethod: The set of spots to consider e.g. ScoreMethod = 'OMP'
%   would find spot from the set o.ompSpotGlobalYX.
% SpotLocation: whether pixel clicked on contains a spot.
% SpotNo: index of nearest spot to cross-hair.
pf = o.CallMethodPrefix(ScoreMethod);

if SpotLocation
    if strcmpi(ScoreMethod,'DotProduct')
        if o.dpSpotScore(SpotNo)>o.CombiQualThresh
            c1 = [0,0.7,0]; else; c1 = [0,0,0];end
        if o.dpSpotScoreDev(SpotNo)<o.CombiDevThresh
            c2 = [1,0,0]; else; c2 = [0,0,0];end
        if o.dpSpotIntensity(SpotNo)<o.CombiIntensityThresh
            c3 = [1,0,0]; else; c3 = [0,0,0];end
        figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
        figtitle.String = sprintf('Spot %.0f is %s: %s{%f %f %f}Score = %.1f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
            SpotNo,o.GeneNames{o.dpSpotCodeNo(SpotNo)},'\color[rgb]',c1,o.dpSpotScore(SpotNo),'\color[rgb]',c2, o.dpSpotScoreDev(SpotNo),...
            '\color[rgb]',c3,o.dpSpotIntensity(SpotNo));
    elseif strcmpi(ScoreMethod,'Prob') || strcmpi(ScoreMethod,'Pixel')
        %Color different parameters depending if over threshold
        if o.([pf,'SpotScore'])(SpotNo)>o.pScoreThresh
            c1 = [0,0.7,0]; else; c1 = [0,0,0];end
        if o.([pf,'LogProbOverBackground'])(SpotNo)<o.pLogProbThresh
            c2 = [1,0,0]; else; c2 = [0,0,0];end
        if o.([pf,'SpotScore'])(SpotNo)+o.([pf,'SpotScoreDev'])(SpotNo)<o.pDevThresh
            c3 = [1,0,0]; else; c3 = [0,0,0];end
        if o.([pf,'SpotIntensity'])(SpotNo)<o.pIntensityThresh
            c4 = [1,0,0]; else; c4 = [0,0,0];end
        figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
        figtitle.String = sprintf('Spot %.0f is %s: %s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProb = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
            SpotNo,o.GeneNames{o.([pf,'SpotCodeNo'])(SpotNo)},'\color[rgb]',c1,o.([pf,'SpotScore'])(SpotNo),'\color[rgb]',c2, o.([pf,'LogProbOverBackground'])(SpotNo),...
            '\color[rgb]',c3,o.([pf,'SpotScoreDev'])(SpotNo),'\color[rgb]',c4,o.([pf,'SpotIntensity'])(SpotNo));
    elseif strcmpi(ScoreMethod,'OMP')
        SpotCodeNo = o.([pf,'SpotCodeNo'])(SpotNo);
        SpotCoefs = full(o.([pf,'Coefs'])(SpotNo,:));
        SpotScore = o.([pf,'SpotScore'])(SpotNo);
        SpotIntensity = o.([pf,'SpotIntensity2'])(SpotNo);
        if isprop(o,'ompNeighbNearPosNeighbMultiplier') && size(o.ompNeighbNonZeros,2)==2
            SpotNeighbNonZero = o.ompNeighbNearPosNeighbMultiplier*...
                o.([pf,'NeighbNonZeros'])(SpotNo,1) + o.([pf,'NeighbNonZeros'])(SpotNo,2);
        else
            SpotNeighbNonZero = o.([pf,'NeighbNonZeros'])(SpotNo);
        end
        figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
        figtitle.String = sprintf('Spot %d: code %d, %s. Coef = %.2f, NeighbNonZero = %d, Score = %.2f, Intensity = %.3f',...
            SpotNo, SpotCodeNo, o.GeneNames{SpotCodeNo}, SpotCoefs(SpotCodeNo), SpotNeighbNonZero,SpotScore,SpotIntensity);
    end
end
end

