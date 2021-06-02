function track_residual(o,TrackInfo)
%% o.track_residual;
% Gives plots of residual spot color, residual in each round and gene
% coefficients at each step of the omp iteration. 
% o: iss object
% TrackInfo contains:
% thresh - omp stops when reduction in residual (norm(r)) falls below
% thresh.
% n_nonzero_coefs - a maximum of n_nonzero_coefs genes can be assigned to each pixel;
% r - r(i,:) is the residual after iteration i.
% i=1 is the initial residual i.e. the SpotColor, 
% i=2 is the residual after all the background vectors have been removed. 
% i=2+j is the residual after the jth gene has been removed. 
% r_norm(i) - the norm of the residual after iteration i. 
% coefs(i,g) - the coefficient for gene g in the iteration i.
% A_omega(:,1:j) contains the gene codes for the genes removed after
% iteration i=2+j.
% omega(1:j) contains the gene numbers for the genes removed after
% iteration i=2+j. 
%% Note that the final iteration is iteration i such that r_norm(i)-r_norm(i+1)<thresh
% or if r_norm(i)-r_norm(i+1)>thresh ∀i, then it is i=n_nonzero_coefs+2. 
ImShape = [o.nBP,o.nRounds];
nGenes = size(TrackInfo.omega,1);
nIters = 2+nGenes;
if TrackInfo.r_norm(nIters-1)-TrackInfo.r_norm(nIters)<TrackInfo.thresh
    FinalIter = nIters-1;
else
    FinalIter = nIters;
end
IterTitle = cell(nIters,1);
IterTitle{1} = 'Initial';
IterTitle{2} = 'Post Background';
for i=3:nIters
    IterTitle{i} = ['Post ', o.GeneNames{TrackInfo.omega(i-2)}];
end
ResClims = [min(TrackInfo.r(:)),max(TrackInfo.r(:))];
BledCodeClims = [min([TrackInfo.A_omega(:);0]),max(TrackInfo.A_omega(:))];
TextColor = 'k';

%% Plots residual spot color at each step of iteration.
if ishandle(35821)
    fig = figure(35822);
else 
    fig = figure(35821);
end
fig.Position = [227,83,1037,669];
j=1;
for i=1:nIters
    subplot(2, nIters, i);
    imagesc(reshape(TrackInfo.r(i,:),ImShape));
    caxis(ResClims);
    colormap(gca,bluewhitered);
    if j<nGenes+1
        gUnbled = reshape(o.UnbledCodes(TrackInfo.omega(j),:),ImShape);
        gSquares = zeros(o.nRounds,4);
        for r=1:o.nRounds
            try
                gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
            end
        end
        for r=1:o.nRounds
            rectangle('Position',gSquares(r,:),'EdgeColor','g',...
                'LineWidth',2,'LineStyle',':')
        end
        hold off
    end
    title(sprintf('%s\nRes = %.2f',IterTitle{i},TrackInfo.r_norm(i)),'Color',TextColor);
    if i==nIters
        colorbar('Position', [0.93  0.55  0.01  0.4]);
    end
    
    if i==nIters-1 && FinalIter==nIters-1
        CoefIter = FinalIter+1;
        TextColor = 'r';
    else
        CoefIter = FinalIter;
        TextColor = 'k';
    end
    
    if i==1
        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'YTickLabel', o.bpLabels);
        set(gca, 'xtick', 1:o.nRounds);
        ylabel('Channel');
        xlabel('Round');
    else
        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'YTickLabel', []);
        set(gca, 'xtick', 1:o.nRounds);
        set(gca, 'XTickLabel', []);
    end
    
    if i>=2 && j<nGenes+1
        subplot(2, nIters, i+nIters);
        imagesc(reshape(TrackInfo.A_omega(:,j),ImShape));
        caxis(BledCodeClims);
        hold on
        for r=1:o.nRounds
            rectangle('Position',gSquares(r,:),'EdgeColor','r',...
                'LineWidth',2,'LineStyle',':')
        end
        hold off
        GeneNo = TrackInfo.omega(j);
        title(sprintf('%.0f: %s\nCoef = %.3f',GeneNo,...
            o.GeneNames{GeneNo},TrackInfo.coefs(CoefIter,GeneNo)),'Color',TextColor);
        if j==1
            set(gca, 'ytick', 1:o.nBP);
            set(gca, 'YTickLabel', o.bpLabels);
            set(gca, 'xtick', 1:o.nRounds);
            ylabel('Channel');
            xlabel('Round');
        else
            set(gca, 'ytick', 1:o.nBP);
            set(gca, 'YTickLabel', []);
            set(gca, 'xtick', 1:o.nRounds);
            set(gca, 'XTickLabel', []);
        end
        j=j+1;
        if j==nGenes+1
            colorbar('Position', [0.93  0.1  0.01  0.4]);
        end
    end
end
sgtitle(sprintf('Residual at each iteration of OMP, Threshold Reduction = %.3f',TrackInfo.thresh));
    
%% Track squared residual in each round at each step of iteration
if ishandle(35823)
    fig2 = figure(35824);
else 
    fig2 = figure(35823);
end
fig2.Position = [227,83,1037,669];
InitialGlyph = 'sq';
InitialColor = 'w';
BackgroundGlyph = 'diamond';
BackgroundColor = [0.6,0.6,0.6];
subplot(2,1,1);
hold on
for i=1:nIters
    Res = reshape(TrackInfo.r(i,:),ImShape);
    rSqNorm = zeros(7,1);
    for r=1:o.nRounds
        rSqNorm(r) = vecnorm(Res(:,r),2,1)'^2;
    end    
    if i>2
        h(i-2) = plot(1:o.nRounds,rSqNorm,'.');
    elseif i==1
        scatter(1:o.nRounds,rSqNorm,'Marker',InitialGlyph,...
            'MarkerFaceColor',InitialColor,'MarkerEdgeColor',InitialColor);
    elseif i==2
        scatter(1:o.nRounds,rSqNorm,'Marker',BackgroundGlyph,...
            'MarkerFaceColor',BackgroundColor,'MarkerEdgeColor',BackgroundColor);
    end
end
if nGenes>0
    legend(h,o.GeneNames(TrackInfo.omega'));
    legend off
    change_gene_symbols(0);
end
allAxesInFigure = findall(fig2,'type','axes');
LegAx = allAxesInFigure(1);
delete(LegAx);
lgnd = legend(IterTitle);
set(lgnd,'color','none','TextColor','w','Location','northeastoutside');
xlabel('Round');
ylabel('Squared Residual');
title('Evolution of Residual in each round with Iteration','Color','w');
xlim([0.5,o.nRounds+0.5]);
xticks(1:o.nRounds);
set(gca,'XColor','w');
set(gca,'YColor','w');
set(gca, 'color', 'k');
hold off

%% Track coefficients of each gene at each step of iteration
subplot(2,1,2);
nCodes = length(o.CharCodes);
GeneNumbers = [TrackInfo.omega',nCodes+1:size(TrackInfo.coefs,2)];
X_values = 1:length(GeneNumbers)+1;
hold on
for i=1:nIters
    UseGene = TrackInfo.coefs(i,GeneNumbers(1:nGenes))~=0;
    UseBackground = TrackInfo.coefs(i,GeneNumbers(nGenes+1:end))~=0;
    UseX = [UseGene,false,UseBackground];
    UseY = [UseGene,UseBackground];
    if i>2
        h2(i-2) = plot(X_values(UseX),TrackInfo.coefs(i,GeneNumbers(UseY)),'.');
    elseif i==1
        scatter(X_values(UseX),TrackInfo.coefs(i,GeneNumbers(UseY)),...
            'Marker',InitialGlyph,'MarkerFaceColor',InitialColor,...
            'MarkerEdgeColor',InitialColor);
    elseif i==2
        scatter(X_values(UseX),TrackInfo.coefs(i,GeneNumbers(UseY)),...
            'Marker',BackgroundGlyph,'MarkerFaceColor',BackgroundColor,...
            'MarkerEdgeColor',BackgroundColor);
    end
end
if nGenes>0
    legend(h2,o.GeneNames(TrackInfo.omega'));
    legend off
    change_gene_symbols(0);
end
allAxesInFigure = findall(fig2,'type','axes');
LegAx = allAxesInFigure(1);
delete(LegAx);
Labels = o.GeneNames(GeneNumbers(1:nGenes));
UseTick = false(size(UseX));
UseTick(1:nGenes) = true;
%Single label for background
UseTick(ceil(nGenes+1+(length(GeneNumbers)-nGenes)/2)) = true;
xlim([min(X_values)-0.5,max(X_values)+0.5]);
xticks(X_values(UseTick));
xticklabels([Labels;'Background']);
xtickangle(90);
dummy_plot = plot(nan,nan,'Marker','none','Color','None');
lgnd2 = legend(dummy_plot,IterTitle(2));
set(lgnd2,'color','none','TextColor','k','Location','northeastoutside','box','off');
xlabel('Gene');
ylabel('Coefficient');
title('Evolution of Gene Coefficients with Iterations','Color','w');
set(gca,'XColor','w');
set(gca,'YColor','w');
set(gca, 'color', 'k');
hold off
end

