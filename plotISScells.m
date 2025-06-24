function plotISScells(DapiImg,o,CellCalled,BoundaryID,SaveF,j)

[ClassCollapse,ColorTable] = ClassProperties('me',0,0);
figure(234321);
o.PlotLineWidth = 1.5;
o.MarkerSize = 7;
o.ompScoreThresh = 5;
o.ompScoreThresh2 = 2;
o.ompIntensityThresh = 0.01;
o.ompIntensityThresh2 = 0.005;
o.ompNeighbThresh = 12;
o.ompNeighbThresh2 = 10;


Roi = round([1, max(o.ompSpotGlobalYX(:,2)), ...
    1, max(o.ompSpotGlobalYX(:,1))]);

o.plot(DapiImg,Roi,'OMP');
daspect([1 1 1])
camroll(90)
o.iss_change_plot('OMP',[],CellCalled.GeneNames )
[h1,h2] = showISScell(BoundaryID,ColorTable,CellCalled,o,50,0);
saveas(h1,fullfile(SaveF,['Cell',num2str(j),'spots_allgenes.png']))
saveas(h2,fullfile(SaveF,['Cell',num2str(j),'pieplot.png']))