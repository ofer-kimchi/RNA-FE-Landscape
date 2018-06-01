function [sequences, listStructBPs,structuresToPrint] = readBPSEQFile(fileName,minDist)

%get data from .bpseq (pseudobase++ datatype) and convert it to our
%sequences and structure

%fileName = '/Users/Ofer/Desktop/pseudoknotSequences.txt'
strucTable = readtable(fileName,'Delimiter',';','ReadVariableNames',false);

sequences = {};
listStructBPs = {};
structuresToPrint = {};

sequence = [];
structBPs = [];
printedStructString = '';

for line = 1:length(strucTable.Var1)
    a = strucTable.Var1{line};
    for i = 1:length(a)
        if a(i) == ' '
            break
        end
        positionInSequence = str2double(a(1:i));
    end
    
    if positionInSequence == 1 && line ~= 1 %then we've reached the end of one sequence and are about to start the next.
        sequences{length(sequences)+1} = sequence; %#ok<*AGROW>
        listStructBPs{length(listStructBPs)+1} = structBPs;
        
%         printedStructCell = cell(1,length(sequence)); %way to print structure. For now, a cell array where the j'th element is equal to k if ntd j is paired with k, and 0 otherwise
%         numBPs = size(structBPs,1);
%         for j = 1:length(printedStructCell)
%             paired = find(structBPs == j);
%             if isempty(paired)
%                 printedStructCell{j} = '0 ';
%             else
%                 if paired <= numBPs
%                     printedStructCell{j} = [num2str(structBPs(paired+numBPs)),' ']; %don't use strcat to preserve spaces
%                 else
%                     printedStructCell{j} = [num2str(structBPs(paired-numBPs)),' ']; 
%                 end
%             end
%         end
%         printedStructString = '';
%         for j = 1:length(printedStructCell)
%             printedStructString = [printedStructString,printedStructCell{j}]; 
%         end
        
        structuresToPrint{length(structuresToPrint)+1} = printedStructString;
        
        sequence = [];
        structBPs = [];
        printedStructString = '';
    end
        
    sequence = [sequence,a(i+1)];
    
    bondedNtd = str2double(a(i+3:end));
    
    printedStructString = [printedStructString, num2str(bondedNtd),' '];
    
    if bondedNtd ~= 0 && ~ any(any(structBPs==bondedNtd)) %don't put in a bond twice, only once.0
        structBPs = [structBPs; positionInSequence, bondedNtd];
    end
    
end

sequences{length(sequences)+1} = sequence; 
listStructBPs{length(listStructBPs)+1} = structBPs;
structuresToPrint{length(structuresToPrint)+1} = printedStructString;

if minDist >0
    seqD = seqpdist(sequences,'Alphabet','NT'); %measures pairwise Jukes-Cantor distances between sequences
    seqDMat = zeros(length(sequences)); %matrix of pairwise distances between sequences
    counter = 0;
    for j = 1:length(sequences)
        for i = j+1:length(sequences)
            counter = counter + 1;
            seqDMat(i,j) = seqD(counter);
            seqDMat(j,i) = seqD(counter); %to make dMat symmetric
        end
    end
    % 
    % seqLengths = zeros(1,length(sequences));
    % for i = 1:length(sequences)
    %     seqLengths(i) = length(sequences{i});
    % end
    % 
    % for i = 1:length(sequences)
    %     for j = 1:length(sequences)
    %         dMat(i,j) = dMat(i,j) * (seqLengths(i) +seqLengths(j))/2; %turns it from a measure of distance per ntd to a measure of total dist
    %     end
    % end

    %now we want to generate a new set of sequences (and structures) that
    %doesn't include sequences that are very near each other in sequence space
    %and structure space


    indexCulledSequences = []; %sequences to remove
    for i = 1:length(sequences)
        seqDMat(i,i) = 1e6; %so we don't worry about diagonal distances being < minDist
    end

    for i = 1:length(sequences)
        if ~any(indexCulledSequences == i) && any(seqDMat(i,:)<minDist) 
            nearbySequences = setdiff(find(seqDMat(i,:)<minDist),indexCulledSequences); %setdiff(A,B) gives elements of A not in B
            %disp([i,nearbySequences])
            indexCulledSequences = [indexCulledSequences, nearbySequences];
        end
    end

    sequences(indexCulledSequences) = [];
    listStructBPs(indexCulledSequences) = [];
    structuresToPrint(indexCulledSequences) = [];
end

t = table(sequences'); t.Properties.VariableNames = {'sequence'};
t2 = table(structuresToPrint'); t2.Properties.VariableNames = {'structure'};
t = [t, t2];

writetable(t,'/Users/Ofer/Dropbox/Brenner/RNA pseudoknots/compareToTristan/pseudoknotSequences160.txt','Delimiter',',');
