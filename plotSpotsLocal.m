function [] = plotSpotsLocal(LocalYXZ,Image,nZ,TileSz,name)
% Used for debugging in detectSpots. Plots image with spots on top
% Plot different plots according to slider location.
S.fh = figure('units','pixels','position',[500 200 800 600],'name',name,'numbertitle','off');  %Left, Bottom, Width, Height
S.Image = Image;
S.Background = imagesc(S.Image(:,:,1)); hold on; colormap hot;
S.Data = LocalYXZ;
ToUse = round(S.Data(:,3)) == 1;
S.x = S.Data(ToUse,2);  % For plotting.  
S.y = S.Data(ToUse,1);
S.LN = scatter(S.x,S.y,'wx');    %Format dots at this point
title(['Z Plane ' num2str(1)]);
xlim([0 TileSz]);
ylim([0 TileSz]);
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
S.Background = imagesc(S.Image(:,:,ZPlane)); hold on; colormap hot;
ToUse = round(S.Data(:,3)) == ZPlane;
S.x = S.Data(ToUse,2);  
S.y = S.Data(ToUse,1);
S.LN = scatter(S.x,S.y,'wx');
title(['Z Plane ' num2str(ZPlane)]);

