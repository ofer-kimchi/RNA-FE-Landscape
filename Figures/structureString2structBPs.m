function strucBPs = structureString2structBPs(stringStructure) 

bondedBPString = [];
endOfBondedBP = false;
currentBP = 0;
strucBPs = [];
for i = 1:length(stringStructure)
    if stringStructure(i) ~= ' '
        bondedBPString = [bondedBPString, stringStructure(i)];
    else
        endOfBondedBP = true;
    end
    
    if i == length(stringStructure)
        endOfBondedBP = true;
    end

    if endOfBondedBP
        bondedBP = str2double(bondedBPString);
        currentBP = currentBP + 1;
        if ~ismember(currentBP,strucBPs) && bondedBP ~=0
            strucBPs = [strucBPs; currentBP bondedBP]; %#ok<*AGROW>
        end
        bondedBPString = [];
        endOfBondedBP = false;
    end
end