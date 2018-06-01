function structureString = structBPs2structureString(numBPs,strucBPs)

structureString = '';
if isempty(strucBPs)
    structureString = repmat('0 ',1,numBPs);
    structureString = structureString(1:end-1); %to get rid of last space
    return
end

for i = 1:numBPs
    iInFirstCol = find(strucBPs(:,1)==i);
    if ~isempty(iInFirstCol)
        iBP = strucBPs(iInFirstCol,2);
    else 
        iInSecondCol = find(strucBPs(:,2)==i);
        if ~isempty(iInSecondCol)
            iBP = strucBPs(iInSecondCol,1);
        else
            iBP = 0;
        end
    end
    structureString = [structureString, num2str(iBP), ' ']; %#ok<AGROW>
end

structureString = structureString(1:end-1); %to get rid of last space