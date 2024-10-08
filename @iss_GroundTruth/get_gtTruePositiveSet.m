function [gtTruePositiveSet,thresh] = get_gtTruePositiveSet(o,thresh,...
    OnlyHighOneChannel,MaxSpots,FractHighChannelThresh)
%% o = o.get_gtTruePositiveSet;
% Returns a logical array that is true for all spots that should be Gad
% i.e. that have high intensity in ground truth round, ground truth
% channel and not in any other channels of ground truth round
% If OnlyHighOneChannel == true then:
% Has to be closer to Line2 (OtherGeneChannelColor=0) than 
% Line1 (OtherGeneChannelColor=GeneChannelOfInterestColor).
% This has to be the case for all but a fraction as indicated by 
% FractHighChannelThresh to be included in the set.
% I.e. want GeneChannelOfInterestColor high and all OtherGeneChannelColor
% low in that imaging round only.
% Else: just require GeneChannelOfInterestColor>thresh.
% thresh(r,b): raw color intensity threshold for round r, channel b.
% If given MaxSpots will alter thresh(r,b) so at most MaxSpots will be
% selected.
if nargin<4 || isempty(MaxSpots)
    MaxSpots = inf;
end

if nargin<5 || isempty(FractHighChannelThresh)
    FractHighChannelThresh = 0;
end

Line1 = [0,0;1,1];
Line2 = [0,0;1,0];
gtTruePositiveSet = cell(o.nRounds+o.nExtraRounds,o.nBP);
for r=o.gtRounds
    for b=1:o.nBP
        if o.gtGeneNo(r,b)==0; continue; end
        if OnlyHighOneChannel==true
            z_scoreSCALE = sqrt(sum(o.HistValues.^2'.*o.gtHistCounts)/sum(o.gtHistCounts(:,b,r)));
            z_scoreColor = o.gt_gtColor{r,b}./z_scoreSCALE;
            CurrentRoundChannel = sub2ind(size(z_scoreColor,2:3),b,r);
            NonGTChannels = [];
            NonGTRounds = [];
            for r2=o.gtRounds
                rNonGTChannels = setdiff(1:o.nBP,o.gtAnchorChannel(r2));
                if r2==o.AnchorRound
                    rNonGTChannels = setdiff(rNonGTChannels,o.DapiChannel);
                end
                if r2==r
                    rNonGTChannels = setdiff(rNonGTChannels,b);
                end
                NonGTChannels = [NonGTChannels,rNonGTChannels];  
                NonGTRounds = [NonGTRounds,ones(size(rNonGTChannels))*r2];
            end
            NonGTRoundChannel = sub2ind(size(z_scoreColor,2:3),NonGTChannels,NonGTRounds);
            nRoundChannel = length(NonGTRoundChannel);
            ToUseAll = false(length(o.gt_gtColor{r,b}),nRoundChannel);
            for r2b2=1:nRoundChannel
                pts = z_scoreColor(:,[CurrentRoundChannel,NonGTRoundChannel(r2b2)]);
                %pts = o.gt_gtColor{r,b}(:,[b,non_gt_channels(b2)],r);
                dist1 = point_to_line(pts,Line1(1,:),Line1(2,:));
                %dist2 = point_to_line(pts,Line2(1,:),Line2(2,:));
                %dist2 = abs(o.gt_gtColor{r,b}(:,non_gt_channels(b2),r));
                dist2 = abs(z_scoreColor(:,NonGTRoundChannel(r2b2)));
                ToUseAll(:,r2b2) = dist2<dist1;
            end
            %ToUse = all(ToUseAll,2);
            ToUse = sum(~ToUseAll,2)<=nRoundChannel*FractHighChannelThresh;
        elseif OnlyHighOneChannel==false
            ToUse = true(length(o.gt_gtColor{r,b}),1);
        end
        gtTruePositiveSet{r,b} = ToUse & o.gt_gtColor{r,b}(:,b,r) > thresh(r,b);
        if sum(gtTruePositiveSet{r,b})>MaxSpots
            SortSpots = sort(o.gt_gtColor{r,b}(ToUse,b,r),'descend');
            thresh(r,b) = SortSpots(MaxSpots);
            gtTruePositiveSet{r,b} = ToUse & o.gt_gtColor{r,b}(:,b,r) > thresh(r,b);
        end
    end
end

end


function d = point_to_line(pt, v1, v2)
% pt should be nx2 
% v1 and v2 are vertices on the line (each 1x2)
% d is a nx1 vector with the orthogonal distances
v1 = [v1,0];
v2 = [v2,0];
pt = [pt,zeros(size(pt,1),1)];
v1 = repmat(v1,size(pt,1),1);
v2 = repmat(v2,size(pt,1),1);
a = v1 - v2;
b = pt - v2;
d = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
end