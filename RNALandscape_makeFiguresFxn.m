function figureFxnTime = RNALandscape_makeFiguresFxn(allFreeEnergies,minFE,...
    storeGraphs,sortedProbs,weightMatrixList,numBondsMatrixList,allGraphs,indexSortedProbs,...
    possiblePermutations,sequence,indexSortedGraphProbs,sortedGraphProbs,STable,numBP,whichStructMinFE)

figureTime = tic;
%% histogram of free energies
figure;
maxFEDisplayed = 30; %arbitrary choice
hist(allFreeEnergies(allFreeEnergies<maxFEDisplayed),min(500,length(allFreeEnergies(allFreeEnergies<maxFEDisplayed))/4));
xlim([floor(minFE),max(maxFEDisplayed,minFE*2)])
set ( gca, 'xdir', 'reverse','FontSize',18 )
xlabel('Free Energy (kcal/mol)')
ylabel('Structures')

%% histogram of structure probabilities where structures are visualized by topology
if storeGraphs
    figure;
    axes('position',[0.1300 0.1100 0.7050 0.7450]) %default is [0.1300 0.1100 0.7750 0.8150]
    bar(sortedProbs(1:min(length(sortedProbs),11)))
    set(gca,'YScale','log','FontSize',18)
    ylabel('probability')
    for i = 1:min(length(sortedProbs),11)
        axes('position',[0.104+i*0.059 0.88 0.05 0.07])
        plot(graph(weightMatrixList{allGraphs(indexSortedProbs(i))})) %plot graph corresponding to structure i
        set(gca, 'XTickLabelMode', 'Manual')
        set(gca, 'XTick', [])
        set(gca, 'YTickLabelMode', 'Manual')
        set(gca, 'YTick', [])
    end
end


%% histogram of structure probabilities
hf = figure('Color','w');
axes('position',[0.1300 0.2100 0.7050 0.7450]) %default is [0.1300 0.1100 0.7750 0.8150]
bar(sortedProbs(1:min(length(sortedProbs),11)))
set(gca,'YScale','log','FontSize',18)
ylabel('probability')
set(gca, 'XTickLabelMode', 'Manual')
set(gca, 'XTick', [])

for m = 1:min(length(sortedProbs),11)
    ha = axes('position',[0.097+m*0.059 0.1 0.06 0.08],'Tag','Bioinfo:rnaplot:circle');
    whichStructPlot = indexSortedProbs(m);
    
    strucUpperDiagMatrix = zeros(numBP);
    structure = cell(1,length(possiblePermutations{whichStructPlot}));
    
    for i = 1:length(structure)
        structure{i} = STable{possiblePermutations{whichStructPlot}(i),2};
        numBonds = length(structure{i})/2;
        for j = 1:numBonds
            k = structure{i}(j);
            l = structure{i}(j+numBonds);
            strucUpperDiagMatrix(k,l) = 1;
        end
    end
    
    rnaplotEdited(strucUpperDiagMatrix,hf,ha,'Sequence',sequence,'Format','Dotdiagram');
    
    
    set(gca, 'XTickLabelMode', 'Manual')
    set(gca, 'XTick', [])
    set(gca, 'YTickLabelMode', 'Manual')
    set(gca, 'YTick', [])
end


%% look at MFE structure in detail
%which structure do you want to look at in detail?

whichStruct = whichStructMinFE; %plot the RNA structure of the min. free energy structure.

strucUpperDiagMatrix = zeros(numBP); %connectivity matrix of structure (strucUpperDiagMatrix_i,j = 1 if ntds i and j are bonded)
%This matrix is used by rnaplot to show the structure
structure = cell(1,length(possiblePermutations{whichStruct})); %list of bonded bps in structure

for i = 1:length(structure)
    structure{i} = STable{possiblePermutations{whichStruct}(i),2};
    numBonds = length(structure{i})/2;
    for j = 1:numBonds
        %set up strucMatrix
        k = structure{i}(j);
        l = structure{i}(j+numBonds);
        strucUpperDiagMatrix(k,l) = 1;
    end
end

rnaplot(strucUpperDiagMatrix,'Sequence',sequence,'Format','Diagram');

%% histogram of topology probabilities
if storeGraphs
    figure;
    axes('position',[0.1300 0.1100 0.7050 0.7450]) %default is [0.1300 0.1100 0.7750 0.8150]
    bar(sortedGraphProbs(1:min(length(sortedGraphProbs),11)))
    set(gca,'YScale','log','FontSize',18)
    set(gca,'YScale','log')
    ylabel('total probability')
    for i = 1:min(length(sortedGraphProbs),11)
        axes('position',[0.104+i*0.059 0.88 0.05 0.07])
        plot(graph(weightMatrixList{indexSortedGraphProbs(i)})) %plot graph corresponding to structure i
        set(gca, 'XTickLabelMode', 'Manual')
        set(gca, 'XTick', [])
        set(gca, 'YTickLabelMode', 'Manual')
        set(gca, 'YTick', [])
    end
end
    
%% a nicer histogram of topology probabilities

hf = figure('Color','w');
set(gcf, 'Position', [400, 400, 1000, 300])
axes('position',[0.1800 0.2600 0.7050 0.6450]) %default is [0.1300 0.1100 0.7750 0.8150]
bar(sortedGraphProbs(1:min(length(sortedGraphProbs),6)))

set(gca,'YScale','log','FontSize',34)
ylabel({'Topology','probability'})
set(gca, 'XTickLabelMode', 'Manual')
set(gca, 'XTick', [])

for m = 1:min(length(sortedProbs),6)
    ha = axes('position',[0.165+m*0.095 0.03 0.07 0.23]); %for 6
    %ha = axes('position',[0.097+m*0.059 0.1 0.06 0.08]); %for 11
    
    topolMat = numBondsMatrixList{indexSortedGraphProbs(m)};
    for i = 1:length(topolMat)
        for j = 1:length(topolMat)
            if topolMat(i,j) == topolMat(j,i) && topolMat(i,j) == 1 && i~=j
                topolMat(j,i) = 0;
            end
        end
    end
    listDoubleBonds1 = [];
    listDoubleBonds2 = [];
    for i = 1:length(topolMat)
        for j = i:length(topolMat)
            if weightMatrixList{indexSortedGraphProbs(m)}(i,j) == 2 && ...
                    numBondsMatrixList{indexSortedGraphProbs(m)}(i,j)< 2
                listDoubleBonds1 = [listDoubleBonds1,i];
                listDoubleBonds2 = [listDoubleBonds2,j];
            end
        end
    end
    g = digraph(topolMat);
    h = plot(g,'ArrowSize',0,'linewidth',2,'EdgeColor','r','NodeColor','g');
    highlight(h,listDoubleBonds1,listDoubleBonds2,'EdgeColor','b','LineWidth',2)
    h.NodeLabel = {};
    axis off
    set(gca, 'XTickLabelMode', 'Manual')
    set(gca, 'XTick', [])
    set(gca, 'YTickLabelMode', 'Manual')
    set(gca, 'YTick', [])
end

%%
figureFxnTime = toc(figureTime);
disp(strcat('Time for making figures= ',num2str(figureFxnTime)));