function [sequence, structure,correctedNC,correctedHairpin] = readCTFile(fileName,correctStructureNonComplementary,correctStructureHairpin) %Convert .ct file to a structure

%fileName must be in .txt format
correctedNC=false;
correctedHairpin=false;

probReadingTable = true;
counter = 0;
while probReadingTable
    if counter==0
        strucTable = readtable(fileName);
    else
        strucTable = readtable(fileName,'HeaderLines',counter+3);
    end
    %From https://rna.urmc.rochester.edu/Text/File_Formats.html:
    % A CT (Connectivity Table) file contains secondary structure information for a sequence. 
    %These files are saved with a CT extension. When entering a structure to calculate the free energy, the following format must be followed.
    % 
    % Start of first line: number of bases in the sequence
    % End of first line: title of the structure
    % Each of the following lines provides information about a given base in the sequence. 
    % Each base has its own line, with these elements in order:
        % Base number: index n
        % Base (A, C, G, T, U, X)
        % Index n-1
        % Index n+1
        % Number of the base to which n is paired. No pairing is indicated by 0 (zero).
        % Natural numbering. RNAstructure ignores the actual value given in natural numbering, so it is easiest to repeat n here.
    % 

    numNtds = strucTable{end,1}';
    if ~isequal(strucTable{:,1}',1:numNtds)
        %disp('problem reading table')
        probReadingTable = true;
        counter = counter + 1;
    else
        probReadingTable = false;
    end
end

if ~isequal(strucTable.Properties.VariableNames,{'Var1','Var2','Var3','Var4','Var5','Var6'})
    disp(['there''s a problem in readCTFile for filename = ',fileName])
    sequence = 'X'; structure = {};
    return
end

if any(any(ismissing(strucTable)))
    disp(['there''s a problem in readCTFile 2 for filename = ',fileName])
    sequence = 'X'; structure = {};
    return
end
    
sequence = '';
for i = 1:numNtds
    if strcmp(strucTable.Var2{i},'a') || strcmp(strucTable.Var2{i},'A')
        sequence = strcat(sequence,'A');
    elseif strcmp(strucTable.Var2{i},'c') || strcmp(strucTable.Var2{i},'C')
        sequence = strcat(sequence,'C');
    elseif strcmp(strucTable.Var2{i},'g') || strcmp(strucTable.Var2{i},'G')
        sequence = strcat(sequence,'G');
    elseif strcmp(strucTable.Var2{i},'u') || strcmp(strucTable.Var2{i},'U')
        sequence = strcat(sequence,'U');
    else
        sequence = strcat(sequence,'O');
    end
end

structure = {};
i=1;
while i <=numNtds
    if strucTable.Var5(i) ~=0 
        %check that we haven't already added this helix (or stem) to structure
        previouslyAdded = false;
        for j = 1:length(structure)
            if any(structure{j} == i) 
                previouslyAdded = true;
            end
        end
        if strucTable.Var5(i) - i == 1 || strucTable.Var5(i) - i == -1
            previouslyAdded = true;
        end
        if correctStructureNonComplementary %if we want to put in corrections to the structure. 
            %We don't allow non-complementary ntds to bond
            if isComplementary(sequence(i),sequence(strucTable.Var5(i)))
                allowedBond = true;
            else
                allowedBond = false;
                correctedNC = true;
            end
        else
            allowedBond = true;
        end
        if ~previouslyAdded && allowedBond
            stemFirst = i;
            stemSecond = strucTable.Var5(i);
            i = i+1;
            if strucTable.Var5(i) == strucTable.Var5(i-1)-1
                isParallel = false;
                continueStem = true;
            elseif strucTable.Var5(i) == strucTable.Var5(i-1)+1
                isParallel = true;
                continueStem = true;
            else 
                continueStem = false;
            end
            if continueStem && correctStructureNonComplementary %if we want to put in corrections to the structure
                %We don't allow non-complementary ntds to bond
                if ~isComplementary(sequence(i),sequence(strucTable.Var5(i)))
                    continueStem = false;
                    correctedNC = true;
                end
            end
            while continueStem
                stemFirst = [stemFirst,i];
                stemSecond = [stemSecond,strucTable.Var5(i)];
                i = i+1;
                if ~isParallel && strucTable.Var5(i) == strucTable.Var5(i-1)-1
                    continueStem = true;
                elseif isParallel && strucTable.Var5(i) == strucTable.Var5(i-1)+1
                    continueStem = true;
                else 
                    continueStem = false;
                end
                if i > strucTable.Var5(i) 
                    continueStem = false;
                end
                for j = 1:length(structure)
                    if any(structure{j} == i) || any(structure{j} == strucTable.Var5(i) )
                        continueStem = false;
                    end
                end
                if continueStem && correctStructureNonComplementary %if we want to put in corrections to the structure
                    %We don't allow non-complementary ntds to bond
                    if ~isComplementary(sequence(i),sequence(strucTable.Var5(i)))
                        continueStem = false;
                        correctedNC = true;
                    end
                end
            end
            stem = [stemFirst,stemSecond];
            structure{length(structure)+1} = stem;
        else
            i = i + 1;
        end
    else 
        i = i + 1;
    end
end



if correctStructureHairpin %if we want to place in corrections by hand to the structure
    %we don't allow hairpin loops to have fewer than three ntds
    minBPInHairpin = 3;
    i = 1;
    while i <= length(structure)
        stem = structure{i};
        if length(stem) > 2
            if stem(end)-stem(end-1) == 1
                isParallel = true;
            elseif stem(end)-stem(end-1) == -1
                isParallel = false;
            else 
                disp('Is this parallel or not??');
            end
            
            checkHairpinLength = true;
            if ~isParallel %if it forms a hairpin
                while checkHairpinLength && ~isempty(stem)
                    stem1 = stem(1:length(stem)/2);
                    stem2 = stem(length(stem)/2+1:end);
                    if stem2(end) - stem1(end) <= minBPInHairpin
                        stem2 = stem2(1:end-1);
                        stem1 = stem1(1:end-1);
                        stem = [stem1,stem2];
                        correctedHairpin = true;
                    else
                        checkHairpinLength = false;
                    end
                end
            else %if it forms an open-net-1
                while checkHairpinLength && ~isempty(stem)
                    stem1 = stem(1:length(stem)/2);
                    stem2 = stem(length(stem)/2+1:end);
                    if stem2(1) - stem1(end) <= minBPInHairpin + length(stem1) 
                        stem2 = stem2(2:end);
                        stem1 = stem1(1:end-1);
                        stem = [stem1,stem2];
                        correctedHairpin = true;
                    else
                        checkHairpinLength = false;
                    end
                end
            end
        end
        if length(stem) == 2
            if stem(2)-stem(1) <= minBPInHairpin
                stem = [];
                correctedHairpin = true;
            end
        end
        if ~isempty(stem)
            structure{i} = stem; %#ok<AGROW>
            i=i+1;
        else
            structure(i) = [];%#ok<AGROW>
        end
    end
end

if correctStructureNonComplementary %check that we really corrected all NC base pairs
    for i = 1:length(structure)
        n = length(structure{i})/2;
        for j = 1:n
            if ~isComplementary(sequence(structure{i}(j)),sequence(structure{i}(j+n)))
                disp(['PROBLEM WITH correctStructureNonComplementary FOR SEQUENCE = ',num2str(sequence)])
            end
        end
    end
end
end


function [allowedBond] = isComplementary(ntd1,ntd2)
allowedBond = false; %can't have a bond forming unless the two ntds are complementary.
    %if ntds are not A,C,G,U, don't allow a bond to form.
    if strcmp(ntd1,'A')
        if strcmp(ntd2,'U')
            allowedBond = true;
        else
            allowedBond = false;
        end
    elseif strcmp(ntd1,'C')
        if strcmp(ntd2,'G')
            allowedBond = true;
        else
            allowedBond = false;
        end 
    elseif strcmp(ntd1,'G')
        if strcmp(ntd2,'C') || strcmp(ntd2,'U')
            allowedBond = true;
        else
            allowedBond = false;
        end 
    elseif strcmp(ntd1,'U')
        if strcmp(ntd2,'G') || strcmp(ntd2,'A')
            allowedBond = true;
        else
            allowedBond = false;
        end 
    end
end