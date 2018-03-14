function[bondEnergy,bondEntropy] = deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrix,bondEntropyMatrix,firstNtd,firstBP,secondNtd,secondBP,ntdType)

bondEnergy = 0;
bondEntropy = 0;
if (firstNtd == 0 || firstBP == 0 || secondNtd == 0 || secondBP == 0) 
    %if one of the sequence elements is unknown, don't do anything
elseif ntdType == 0
    %For the bondEnergyMatrix,
    %First index tells you if the first bp of the set is AU (1) CG (2) GC (3)
    %UA (4) GU (5) or UG (6)
    %Second index tells you if the 3' ntd of the second bp is A (1) C (2) G(3)
    %or U(4)
    %Third index tells you if the 5' ntd of the second bp is A (1) C (2) G(3)
    %or U(4)
    if firstNtd == 1 && firstBP == 4
        basePair = 1;
    elseif firstNtd == 2 && firstBP == 3
        basePair = 2;
    elseif firstNtd == 3 && firstBP == 2
        basePair = 3;
    elseif firstNtd == 4 && firstBP == 1
        basePair = 4;
    elseif firstNtd == 3 && firstBP == 4
        basePair = 5;
    elseif firstNtd == 4 && firstBP == 3
        basePair = 6;
    else
        disp('The deltaG function has a problem RNA/RNA')
    end
    bondEnergy = bondEnergy + bondEnergyMatrix(basePair,secondNtd,secondBP);
    bondEntropy = bondEntropy + bondEntropyMatrix(basePair,secondNtd,secondBP);
elseif ntdType == 1
    %For the bondEnergyMatrix,
    %First index tells you if the first bp of the set is AU (1) CG (2) GC (3)
    %UA (4) GU (5) or UG (6)
    %Second index tells you if the 3' ntd of the second bp is A (1) C (2) G(3)
    %or U(4)
    %Third index tells you if the 5' ntd of the second bp is A (1) C (2) G(3)
    %or U(4)
    if firstNtd == 1 && firstBP == 8
        basePair = 1;
    elseif firstNtd == 2 && firstBP == 7
        basePair = 2;
    elseif firstNtd == 3 && firstBP == 6
        basePair = 3;
    elseif firstNtd == 4 && firstBP == 5
        basePair = 4;
    elseif secondNtd == 1 && secondBP == 8
        basePair = 5;
    elseif secondNtd == 2 && secondBP == 7
        basePair = 6;
    elseif secondNtd == 3 && secondBP == 6
        basePair = 7;
    elseif secondNtd == 4 && secondBP == 5
        basePair = 8;
    else
        disp('The deltaG function has a problem RNA/DNA')
    end
    bondEnergy = bondEnergy + bondEnergyMatrix(basePair,secondNtd,secondBP-4);
    bondEntropy = bondEntropy + bondEntropyMatrix(basePair,secondNtd,secondBP-4);
elseif ntdType == 2
    %For the bondEnergyMatrix. We're dealing with DNA/DNA pair
    %First index tells you if the first bp of the set is AT (1) CG (2) GC
    %(3) or TA (4)
    %Second index tells you if the 3' ntd of the second bp is A (1) C (2) G(3)
    %or T(4)
    %Third index tells you if the 5' ntd of the second bp is A (1) C (2) G(3)
    %or T(4)
    if firstNtd == 5 && firstBP == 8
        basePair = 1;
    elseif firstNtd == 6 && firstBP == 7
        basePair = 2;
    elseif firstNtd == 7 && firstBP == 6
        basePair = 3;
    elseif firstNtd == 8 && firstBP == 5
        basePair = 4;
    else
        disp('The deltaG function has a problem DNA/DNA')
    end
    bondEnergy = bondEnergy + bondEnergyMatrix(basePair,secondNtd-4,secondBP-4);
    bondEntropy = bondEntropy + bondEntropyMatrix(basePair,secondNtd-4,secondBP-4);
end