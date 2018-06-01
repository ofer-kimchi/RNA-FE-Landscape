function [weightMatrix,numBondsMatrix] = structureString2graph(stringStructure) 
%%
%takes as input a string '10 9 8 0 0 0 3 2 1' which means ntd 1 is bonded
%to 10, 2 to 9, 3 to 8, 4-7 are unbonded.
%Returns a weightMatrix and numBondsMatrix representing the topology of
%this structure

strucBPs = structreString2structBPs(stringStructure);
structure = structBPs2structure(strucBPs);

numBP = length(stringStructure(stringStructure == ' '));
if stringStructure(end) ~= ' '
    numBP = numBP + 1;
end

strucMatrix = zeros(numBP);
for i = 1:size(strucBPs,1)
    strucMatrix(strucBPs(i,1),strucBPs(i,2)) = 1;
    strucMatrix(strucBPs(i,2),strucBPs(i,1)) = 1;
end

allPerms = cell(1,4);
for i = 1:4
    allPerms{i} = perms(1:i);
end

%%
[~,weightMatrix,numBondsMatrix] = ...
    calculateEntropy_10_25_17_noGraph(structure,strucMatrix,numBP,zeros(1,numBP),3,1,1,allPerms,1);
