function structure = structBPs2structure(strucBPs)
%%
structure = {};
numStems = 0;
startOfStem = false;
stem1 = [];
stem2 = [];
isParallel = 0;
for i = 1:size(strucBPs,1)
    firstBP = strucBPs(i,1);
    secondBP = strucBPs(i,2);
    
    if i ~= 1
        if firstBP - stem1(length(stem1)) == 1
            if secondBP - stem2(length(stem2)) == 1 && (isParallel == 1 || isParallel == 0)
                startOfStem = false;
                if isParallel == 0 %meaning the stem had been too short to tell if it was parallel or anti-parallel
                    isParallel = 1;
                end
            elseif secondBP - stem2(length(stem2)) == -1 && (isParallel == -1 || isParallel == 0)
                startOfStem = false;
                if isParallel == 0 %meaning the stem had been too short to tell if it was parallel or anti-parallel
                    isParallel = -1;
                end
            else
                startOfStem = true;
            end
        else
            startOfStem = true;
        end
    end
    
    if ~startOfStem 
        stem1 = [stem1, firstBP]; %#ok<*AGROW>
        stem2 = [stem2, secondBP];
    else
        numStems = numStems + 1;
        structure{numStems} = [stem1, stem2];
        stem1 = firstBP;
        stem2 = secondBP;
        isParallel = 0;
    end

    if i == size(strucBPs,1)
        numStems = numStems + 1;
        structure{numStems} = [stem1, stem2];
    end
end