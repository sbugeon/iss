function Imgs = plot_spot_colors_grid(o, SpotColors, PointCorrectedLocalYX, ImSz, Dist,...
    SpotCodeNo, Clim, IncludeGT, Filter)
%% plot_spot_colors_grid(o, SpotColors, PointCorrectedLocalYX, ImSz, Dist,...
%    SpotCodeNo, Clim, IncludeGT, Filter)
% Once SpotColors found for a subset of pixels by get_spot_colors_grid,
% this produces the plot of each round and channel.
% o: iss object
% SpotColors(s,b,r): gives intensity for spot s in channel b, round r.
% PointCorrectedLocalYX(s,:,r,b): gives local coordinate of spot s but
%   transformed to round r, channel b according to o.D(:,:,t,r,b).
% ImSz: radius of image that is plotted for each round and channel.
% Dist: distance of nearest spot to cross-hair
% SpotCodeNo: gene assigned to nearest spot
% Clim(1,b,r) gives the minimum intensity to be shown by the image of
%   channel b, round r. Clim(2,b,r) gives the maximum. If not given, then
%   will use min/max of SpotColors(:,b,r).
% IncludeGT: true if also want to plot the SpotColors in the ground
%   truth rounds as indicated by o.gtRounds.
% Filter: true if SpotColors from the filtered images i.e. the
%   tiles in o.TileDirectory. false if from the raw images i.e. directly
%   from the nd2 file (and then focus-stacked).

if nargin<7
    Clim = [];
end
if nargin<8 || isempty(IncludeGT)
    IncludeGT = false;
end
if nargin<9 || isempty(Filter)
    Filter = true;
end

numCharCode = str2double(regexp(cell2mat(o.CharCodes(SpotCodeNo)),'\d','match'))+1;
if IncludeGT
    o.UseRounds = [o.UseRounds,o.gtRounds];
    numCharCode = [numCharCode,990*o.gtRounds];
end
nRounds = max(o.UseRounds);

%If Dist to nearest spot is less than MaxDist, then the crosshair will
%be red on the squares corresponding to the gene the nearest spot
%matches to.
MaxDist = 10;
Ylegends = {o.bpLabels{:}};
Xlegends = string(1:nRounds);
ColorCh = cellfun(@str2num,o.bpLabels);

Imgs = cell(nRounds,o.nBP);
for r=1:nRounds
    for b=1:o.nBP
        c = ColorCh(b)+1;
        h = subplot(o.nBP, nRounds, (b-1)*nRounds + r);
        if r == 1 && b == 1
            Pos1 = get(h,'position');
        end
        if r == 1 && b == o.nBP
            Pos2 = get(h,'position');
        end
        if r == nRounds && b == o.nBP
            Pos3 = get(h,'position');
        end
        if ~ismember(b,o.UseChannels) || ~ismember(r,o.UseRounds)
            xticks([]);
            yticks([]);
            continue;
        end
        
        rbYX = round(PointCorrectedLocalYX(1,:,r,b));
        y0 = rbYX(1);
        x0 = rbYX(2);
        if y0>o.TileSz || y0<1 || x0>o.TileSz || x0<1
            continue;
        end
        y1 = max(1,y0 - ImSz);
        y2 = min(o.TileSz,y0 + ImSz);
        x1 = max(1,x0 - ImSz);
        x2 = min(o.TileSz,x0 + ImSz);
        BaseIm = reshape(SpotColors(:,b,r),[ImSz*2+1,ImSz*2+1]);
        Imgs{b,r} = BaseIm;
        imagesc([x1 x2], [y1 y2], BaseIm); hold on
        if Filter
            if ~isempty(Clim) && sum(size(Clim) == [2,o.nBP,nRounds])==3
                if IncludeGT && ismember(r,o.gtRounds)
                else
                    caxis([Clim(1,b,r),Clim(2,b,r)]);
                end
            end
            colormap(gca,bluewhitered);
        end
        axis([x0-ImSz, x0+ImSz, y0-ImSz, y0+ImSz]);
        colorbar;
        if numCharCode(r)==c && Dist<MaxDist 
            ax = gca;
            ax.XColor = 'r';
            ax.YColor = 'r';
            plot(xlim, [y0 y0], 'g'); plot([x0 x0], ylim, 'g');
        elseif IncludeGT && o.gtGeneNo(r,b)>0
            ax = gca;
            if Dist<MaxDist && SpotCodeNo==o.gtGeneNo(r,b)
                ax.XColor = 'r';
                ax.YColor = 'r';
                plot(xlim, [y0 y0], 'g'); plot([x0 x0], ylim, 'g');
            else
                ax.XColor = 'y';
                ax.YColor = 'y';
                plot(xlim, [y0 y0], 'y'); plot([x0 x0], ylim, 'y');
            end
        else
            plot(xlim, [y0 y0], 'k'); plot([x0 x0], ylim, 'k');
        end
        if r==1; ylabel(Ylegends{b},'Color',[0.15 0.15 0.15]); end
        if b==o.nBP; xlabel(Xlegends(r),'Color',[0.15 0.15 0.15]); end
        set(gca, 'YDir', 'normal');
        hold off
    end
end
PosDev = 0.02;
SuperAxisPos = [Pos2(1:2)-PosDev,Pos3(1)+Pos3(2)-Pos2(1)+PosDev*2,Pos1(2)+Pos1(4)-Pos3(2)+PosDev*2];
hSuper=axes('position',SuperAxisPos,'visible','off');
hSuper.XLabel.Visible='on';
hSuper.YLabel.Visible='on';
axes(hSuper);
ylabel('Channel');
xlabel('Round');
end



