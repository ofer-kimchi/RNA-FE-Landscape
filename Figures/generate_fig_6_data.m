
T=300; vs = 0.02;sequenceFrequencies = 1; c1 = 1; c2 = 1; duplexEntropyPenaltyInKB = 17; 
numParForLoops = 5;pairwiseEnergies = 1; storeGraphs = 1; makingFigures = 0;
printProgressUpdate = 1;allowParallelStrands = 0;allowPseudoknots = 1;minNumCompatible = 0;

minBPInRegionVec = [3,4,4,4,3,4]; 
substemsVec = {'all','all','all','all','all','all'}; 
sequencesVec = {'CGCCCGAACUAAGCGCCCGGAAAAAGGCUUAGUUGACGAGGAUGGAGGUUAUCGAAUUUCGCGGAUCCUCCCG',... %top middle
    'UGUGUCUUGGAUCGCGCGGGUCAAAUGUAUAUGGUUCAUAUACAUCCGCAGGCACGUAAUAAAGCGA',...%top left
    'CCGGCUGAGUGUGCAGAUCACAGCCGUAAGAAUUUCUUCAAACCAAGGGGGUGACUCCUUGAACAAAGAGAAAUCACAUGAUCUU',... %bottom left
    'GGAGACCCAUGGAGAGUCUCCAGAGGUUCCUAAGGCCUCGAUCGCGCCUGCCGGGAGGCAGAAUGUCCCGGUUCUCC',... %bottom middle
    'UUGUCAUAUCUGGAUCCAACAGUUAAACCAUGUGAUGGUGUAUACUGUGGUAUGGCGUAAAACAUCGGAG',... %top right
    'GGGAAACAACAGGAGGGGGCCACGUGUGGUGCCGUCCGCGCCCCCUAUGUUGUAACAGAAGCACCACC'}; %bottom right
realDotBracketVec =...
   {'.(((((((((((((............)))))))).........[[[[[....(((...))))))))]]]]]..',...
    '(((((((....[[[[((((.....((((((((.....))))))))))))))))))).......]]]]',...
    '.((((((.......[[[[[[)))))).....(((((((.....(((((((....)))))))......)))))))....]]]]]].',...
    '((((((....[[[[[)))))).(((((.......)))))..........(((((((.........)))))))]]]]]',...
    '.(((((((((....[[[.(((((...((((.....))))....)))))))))))))).........]]].',...
    '.............(((((((.....[[[[[[[.......)))))))..............]]]]]]].'};

for seqIndex = 1:6
    minBPInRegion = minBPInRegionVec(seqIndex); substems = substemsVec{seqIndex};
    sequence = sequencesVec{seqIndex};
    realDotBracket = realDotBracketVec{seqIndex};

    [sortedProbs, indexSortedProbs,sortedFE, indexSortedFE,graphProbs,weightMatrixList,numBondsMatrixList,startAndPermuTime,checkFxnTime,totalTime,STable,possiblePermutations,numRegions,expFEeff] = ...
        RNALandscape_main({sequence},sequenceFrequencies,c1,c2,T,vs,duplexEntropyPenaltyInKB,minBPInRegion,numParForLoops,pairwiseEnergies,storeGraphs,makingFigures,...
        printProgressUpdate,allowParallelStrands, allowPseudoknots, minNumCompatible, substems);
    
    [sortedGraphProbs,indexSortedGraphProbs] =sort(graphProbs,'descend');
    sortedWMList = weightMatrixList(indexSortedGraphProbs);
    sortedNBMList = numBondsMatrixList(indexSortedGraphProbs);

    %find MFE dotbracket notation
    structure = cell(1,length(possiblePermutations{indexSortedProbs(1)})); 
    for i = 1:length(structure)
        structure{i} = STable{possiblePermutations{indexSortedProbs(1)}(i),2};
        numBonds = length(structure{i})/2;
        for j = 1:numBonds 
            k = structure{i}(j);
            l = structure{i}(j+numBonds);
        end
    end

    dotBracketMFE = structBPs2dotBracket2(length(sequence),structure)
    realStrucBPs = dotBracket2structBPs(realDotBracket);
    realStrucString = structBPs2structureString(length(sequence),realStrucBPs);
    [realWM,realNBM] = structureString2graph(realStrucString);
    expWM = realWM; expNBM = realNBM;
    for j = 1:length(sortedGraphProbs)
        if myPseudoIsIsomorphic(expWM,sortedWMList{j},expNBM,sortedNBMList{j})
            probRealTopology = sortedGraphProbs(j);
            whichRealTopology = j %which index of topology is the one found experimentally?
            break
        end
    end


    hf = figure('Color','w');
    set(gcf, 'Position', [400, 400, 1000, 300])
    axes('position',[0.1800 0.2600 0.7050 0.6450]) %default is [0.1300 0.1100 0.7750 0.8150]
    bar(sortedGraphProbs(1:min(length(sortedGraphProbs),6)))

    set(gca,'YScale','log','FontSize',34)
    ylabel('Probability')
    set(gca, 'XTickLabelMode', 'Manual')
    set(gca, 'XTick', [])

    for m = 1:min(length(sortedProbs),6) %plot the topologies (not trivial in MatLab to plot them in a nice way)
        ha = axes('position',[0.165+m*0.095 0.03 0.07 0.23]); %showing top 6 topologies per sequence
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
        h = plot(g,'ShowArrows','off','linewidth',2,'EdgeColor',rgb('salmon'),'NodeColor',rgb('apple green'));
        highlight(h,listDoubleBonds1,listDoubleBonds2,'EdgeColor',rgb('bright blue'),'LineWidth',2)
        h.NodeLabel = {};
        axis off
        set(gca, 'XTickLabelMode', 'Manual')
        set(gca, 'XTick', [])
        set(gca, 'YTickLabelMode', 'Manual')
        set(gca, 'YTick', [])
    end
end
