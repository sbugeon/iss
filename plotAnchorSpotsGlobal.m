function [] = plotAnchorSpotsGlobal(GlobalYXZ,name)
% Given the global coordinates GlobalYXZ, this plots the positions of all
% the spots found
% Plot different plots according to slider location.
%MAY NEED TO PUT AN ISEMPTY THING LIKE IN plotSpotsResolved.m
S.fh = figure('units','pixels','position',[500 200 800 600],'name',name,'numbertitle','off');  %Left, Bottom, Width, Height
S.Data = GlobalYXZ;
ToUse = round(S.Data(:,3)) == 1;
S.x = S.Data(ToUse,2);  % For plotting.  
S.y = S.Data(ToUse,1);
S.LN = scatter(S.x,S.y);    %Format dots at this point
title(['Z Plane ' num2str(1)]);
XLim = max(S.Data(:,2))-mod(max(S.Data(:,2)),500)+500;
YLim = max(S.Data(:,1))-mod(max(S.Data(:,1)),500)+500;
xlim([0 XLim]);
ylim([0 YLim]);
nZ = round(max(GlobalYXZ(:,3)));
S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[60 8 693 18],...
                 'min',1,'max',nZ,'val',1,...
                 'sliderstep',[1/(nZ-1) 1/(nZ-1)],...
                 'callback',{@sl_call,S});  
set( findall( S.fh, '-property', 'Units' ), 'Units', 'Normalized' )


function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
ZPlane = round(get(h,'value'));
ToUse = round(S.Data(:,3)) == ZPlane;
set(S.LN,'xdata',S.Data(ToUse,2));
set(S.LN,'ydata',S.Data(ToUse,1));
title(['Z Plane ' num2str(ZPlane)]);