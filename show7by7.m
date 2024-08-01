function show7by7(o,Size,spath,I)
%% this scrip will produce a 7 color image for each round after alignment, 
% and the final spot calling results around a given point of interest selected with the cross-hair

% I background image
% Size : pixels to show around point of interest
% spath: where to save the data
%%
o.ompScoreThresh = 5; % more stringent threshold, need to be true for only one of the three
o.ompScoreThresh2 = 5; % less stringent threshold, need to be true for all three
o.ompIntensityThresh = 0.01; % more stringent threshold, need to be true for only one of the three
o.ompIntensityThresh2 = 0.01; % less stringent threshold, need to be true for all three
o.ompNeighbThresh = 12; % more stringent threshold, need to be true for only one of the three
o.ompNeighbThresh2 = 12; % less stringent threshold, need to be true for all three

o.MarkerSize = 5;
o.PlotLineWidth = 1.2;
Roi = round([1, max(o.dpSpotGlobalYX(:,2)), ...
    1, max(o.dpSpotGlobalYX(:,1))]);
o.plot(I,Roi,'OMP');
daspect([1 1 1])
camroll(90)
[~,Imgs,YX] = iss_view_spot_omp3(o, 234321,Size);
YX = round(YX);

figure(234321)
xlim([YX(1)-Size+0.5 YX(1)+Size+0.5])
ylim([YX(2)-Size+0.5 YX(2)+Size+0.5])
set(gca,'Clipping','on')
axis off
set(gcf,'Position',[0.147265625,0.359722222222222,0.550390625,0.503472222222222])


F = getframe(gca);
Image = frame2im(F);
imwrite(Image, fullfile(spath,'ISS_result.png'))
S = size(Image);

% AF405: rouge
% AF488: vert
% Dy485xl: jaune
% Cy3: bleu
% TexasRED: cyan
% Cy5: violet
% Atto425: white

Ca = [0.3 0.3 0.4 0.4 0.5 0.2 0.3];
% BigN = 10^100;
% Ca = [BigN 0.3 0.4 BigN BigN BigN 0.3];
sub = [95 70 70 70 70 70 40];
clear options;
options.message   = false;
options.append = true;
Imgs2 = Imgs;

for i=1:size(Imgs,2)
    for j=1:size(Imgs,1)
        Imgs2{j,i}=Imgs2{j,i} - prctile(Imgs2{j,i}(:),sub(j));
        Imgs2{j,i} = Imgs2{j,i}/(Ca(j)*max(Imgs2{j,i}(:)));
    end
    RGB_im = zeros(size(Imgs2{j,i},1),size(Imgs2{j,i},2),3);
    RGB_im(:,:,1) = Imgs2{1,i} +Imgs2{3,i} + Imgs2{7,i}  +  Imgs2{6,i};% red
    RGB_im(:,:,2) = Imgs2{2,i} + Imgs2{3,i} + Imgs2{7,i}  + Imgs2{5,i};% green
    RGB_im(:,:,3) =  Imgs2{4,i}  + Imgs2{5,i}  +Imgs2{6,i} + Imgs2{7,i} ;% blue
   
    figure
     imagesc(fliplr(rot90(RGB_im,1)))
    daspect([1 1 1])
    axis off
    set(gcf,'Position',[616,600,944,738])
    F = getframe(gca);
    Image = frame2im(F);S2 =size(Image);
    Image = imresize(Image,S(1)/S2(1));
    imwrite(Image, fullfile(spath,['Round',num2str(i),'.png']))

end