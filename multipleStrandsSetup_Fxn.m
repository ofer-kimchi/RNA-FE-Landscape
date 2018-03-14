function [sequence,sequenceInNumbers,numSequences,numBP,sequence1,numBP1,linkerPos,sequence2,numBP2,sequence3,numBP3] = ...
    multipleStrandsSetup_Fxn(sequences)

numSequences = length(sequences);

sequence1 = sequences{1}; numBP1 = length(sequence1);
sequence2 = ''; sequence3 = '';
if numSequences == 1
    linkerPos = {};
    sequence = sequence1;
    numBP2 = 0; numBP3 = 0; %needed if we're using parfor function. I don't know why.
    isDNA = zeros(1,length(sequence)); %which ntds are DNA (rather than RNA). 1 for ntds which are DNA, 0 for RNA.
    if any(sequence1 == 'T')
        isDNA(1:end) = 1;
    end
elseif numSequences == 2
    sequence2 = sequences{2}; numBP2 = length(sequence2);
    linker = repmat('O',1,length(sequence1)+length(sequence2)); %want linker to not be allowed to bind, and long enough not to constrain anything
    sequence = [sequence1,linker,sequence2];
    linkerPos = {numBP1+1:numBP1+length(linker)}; %ntds of sequence which are actually the linker
    numBP3 = 0;
    isDNA = zeros(1,length(sequence)); %which ntds are DNA (rather than RNA). 1 for ntds which are DNA, 0 for RNA.
    if any(sequence1 == 'T')
        isDNA(1:numBP1) = 1;
    end
    if any(sequence2 == 'T')
        isDNA(linkerPos{1}(end)+1:end) = 1;
    end
    
elseif numSequences == 3
    sequence2 = sequences{2}; numBP2 = length(sequence2);
    sequence3 = sequences{3}; numBP3 = length(sequence3);
    linker = repmat('O',1,length(sequence1)+length(sequence2)+length(sequence3));
    sequence = [sequence1,linker,sequence2,linker,sequence3];
    linkerPos = {numBP1+1:numBP1+length(linker),length([sequence1,linker,sequence2])+1:length([sequence1,linker,sequence2])+length(linker)};
    isDNA = zeros(1,length(sequence)); %which ntds are DNA (rather than RNA). 1 for ntds which are DNA, 0 for RNA.
    if any(sequence1 == 'T')
        isDNA(1:numBP1) = 1;
    end
    if any(sequence2 == 'T')
        isDNA(linkerPos{1}(end)+1:linkerPos{2}(1)-1) = 1;
    end
    if any(sequence3 == 'T')
        isDNA(linkerPos{2}(end)+1:end) = 1;
    end
end

numBP = length(sequence);

%write the sequence as a vector where each instance of A is replaced with 1, 
%C with 2, G with 3, U with 4, and other with 0.
sequenceInNumbers = zeros(1,numBP); 
for i =1:numBP
    if sequence(i) == 'A'
        sequenceInNumbers(i) = 1;
    elseif sequence(i) == 'C'
        sequenceInNumbers(i) = 2;
    elseif sequence(i) == 'G'
        sequenceInNumbers(i) = 3;
    elseif sequence(i) == 'U' || sequence(i) == 'T' 
        sequenceInNumbers(i) = 4;
    end
end
sequenceInNumbers = sequenceInNumbers + 4.*isDNA; %so DNA ntds become 5, 6, 7, 8 for A, C, G, T.
