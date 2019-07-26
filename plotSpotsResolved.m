function [] = plotSpotsResolved(o,GlobalYXZ,Good,Tiles,name)
% Used for debugging in detectSpots. Plots image with spots on top
% Plot different plots according to slider location.
S.fh = figure('units','pixels','position',[500 200 800 600],'name',name,'numbertitle','off');  hold on; set(gca, 'color', 'k');
SquareX1 = [0, 0, o.TileSz];
SquareY1 = [o.TileSz, 0, 0];
SquareX2 = [o.TileSz, o.TileSz, 0];
SquareY2 = [0, o.TileSz, o.TileSz];

SquareColors = hsv2rgb([(o.UseRounds)'/size(o.UseRounds,2), [.5, .6] .*ones(size(o.UseRounds,2),1)]);
SquareColors(o.ReferenceRound,:)=1.0;
for r=[o.UseRounds,o.ReferenceRound]            %Add anchor round in as treated as extra round. Wont be in UseRounds
        for t=Tiles
            MyOrigin = o.TileOrigin(t,:,r);
            plot(SquareX1 + MyOrigin(2), SquareY1 + MyOrigin(1),...
                '--', 'Color', SquareColors(r,:));
            plot(SquareX2 + MyOrigin(2), SquareY2 + MyOrigin(1),...
                ':', 'Color', SquareColors(r,:));

            text(MyOrigin(2), MyOrigin(1),...
                sprintf('T%d r%d', t, r), 'color', SquareColors(r,:)); 
        end
end   


S.Data = GlobalYXZ;
S.Good = Good;

title(['Z Plane ' num2str(1)]);
XLim = max(S.Data(:,2))-mod(max(S.Data(:,2)),500)+500;
YLim = max(S.Data(:,1))-mod(max(S.Data(:,1)),500)+500;
xlim([0 XLim]);
ylim([0 YLim]);

ToUseGood = round(S.Data(:,3)) == 1 & S.Good;
S.xGood = S.Data(ToUseGood,2);  % For plotting.  
S.yGood = S.Data(ToUseGood,1);
if isempty(S.xGood)
    S.xGood = XLim+1;
    S.yGood = YLim+1;
end
S.LNGood = plot(S.xGood,S.yGood, 'b.', 'markersize', 1);    %Format dots at this point

ToUseBad = round(S.Data(:,3)) == 1 & ~S.Good;
S.xBad = S.Data(ToUseBad,2);  % For plotting.  
S.yBad = S.Data(ToUseBad,1);
if isempty(S.xBad)      %If 1st graph has no data, slider doesn't work
    S.xBad = XLim+1;
    S.yBad = YLim+1;
end
S.LNBad = plot(S.xBad,S.yBad, 'r.', 'markersize', 1);


S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[60 8 693 18],...
                 'min',1,'max',o.nZ,'val',1,...
                 'sliderstep',[1/(o.nZ-1) 1/(o.nZ-1)],...
                 'callback',{@sl_call,S});  
set( findall( S.fh, '-property', 'Units' ), 'Units', 'Normalized' )




function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
ZPlane = round(get(h,'value'));

ToUseGood = round(S.Data(:,3)) == ZPlane & S.Good;
set(S.LNGood,'xdata',S.Data(ToUseGood,2));
set(S.LNGood,'ydata',S.Data(ToUseGood,1));

%S.LNGood = plot(S.xGood,S.yGood, 'b.', 'markersize', 1);    %Format dots at this point

ToUseBad = round(S.Data(:,3)) == ZPlane & ~S.Good;
set(S.LNBad,'xdata',S.Data(ToUseBad,2));
set(S.LNBad,'ydata',S.Data(ToUseBad,1));
%S.LNBad = plot(S.xBad,S.yBad, 'r.', 'markersize', 1);

title(['Z Plane ' num2str(ZPlane)]);
drawnow;

