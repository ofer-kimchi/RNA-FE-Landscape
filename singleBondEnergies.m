function[bondEnergyLocal] = singleBondEnergies(structure,sequenceInNumbers)

bondEnergyLocal=0; %the bond energy computed within this function

bondEnergyMatrix = ...
    [[-2.06209  -0.182809  -2.06209   -2.06209];...
    [-0.182809  0          -4.86521   -2.06209];...
    [-2.06209   -4.86521   0          -2.06209];
    [-2.06209   -2.06209   -2.06209   0]];
%bondEnergyMatrix(i,j) gives energy for bond between nucelotides i and j
%where A=1, C=2, G=3, U=4. This comes from Cragnolini's paper
%"Coarse-grained HiRE-RNA model for ab initio RNA folding beyond simple 
%molecules, including noncanonical and multiple base pairings"

for i = 1:length(structure)
    numBonds = length(structure{i})/2;
    for j = 1:numBonds  
        
        k = structure{i}(j);
        l = structure{i}(j+numBonds);
        firstNtd = sequenceInNumbers(k); %what is the 5' ntd?
        firstBP = sequenceInNumbers(l); %what is it bonded to?
        
        if firstNtd ~=0 && firstBP ~= 0
            bondEnergyLocal = bondEnergyLocal + bondEnergyMatrix(firstNtd,firstBP);
        end
    end
end
