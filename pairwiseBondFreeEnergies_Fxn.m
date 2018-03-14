function[bondEnergyLocal,bondEntropyLocal] = ...
    pairwiseBondFreeEnergies_Fxn(structure,sequenceInNumbers,linkers,T)
%%
numBP = length(sequenceInNumbers);
bondEnergyLocal=0; %the bond energy computed within this function
bondEntropyLocal=0; %the bond entropy computed within this function

if isempty(structure) 
    return
end

realNtds = 1:length(sequenceInNumbers);
for i = 1:length(linkers)
    realNtds = setdiff(realNtds,linkers{i});
end

firstNtds = 1;
for i = 1:length(linkers)
    firstNtds = [firstNtds,linkers{i}(end)+1];  %#ok<AGROW>
end

lastNtds = numBP;
for i = 1:length(linkers)
    lastNtds = [lastNtds,mod(linkers{i}(1)-1,numBP+1)]; %#ok<AGROW> %returns linkers{i}(1)-1 unless it's 0, in which case it returns numBP];
end


%%
bondFreeEnergyMatrixRNARNA = zeros(6,4,4); %matrix of stacking energies (deltaG). Units are kcal/mol.
bondEnergyMatrixRNARNA = zeros(6,4,4); %matrix of enthalpies (deltaH). Units are kcal/mol.
bondEntropyMatrixRNARNA = zeros(6,4,4); %matrix of entropies (deltaS). Units are initially eu, but then converted to kcal/(mol*K).
%Sources: Table 4 of Xia et al Biochemistry '98
    %Table 4 of Mathews et al. JMB '99
    %Table 3 of Xia, Mathews, Turner "Thermodynamics of RNA secondary structure formation" in book by Soll, Nishmura, Moore
%First index tells you if the first bp of the set is AU (1) CG (2) GC (3)
%UA (4) GU (5) or UG (6)
%Second index tells you if the 3' ntd of the second bp is A (1) C (2) G(3)
%or U(4) (which row of table 1 in Serra & Turner).
%Third index tells you if the 5' ntd of the second bp is A (1) C (2) G(3)
%or U(4) (which column of table 1 in Serra & Turner).

bondEnergyMatrixRNARNA(1,:,:) = [-3.9, 2, -3.5, -6.82; -2.3, 6, -11.4, -0.3; -3.1, -10.48, -3.5, -3.21; -9.38, 4.6, -8.81, -1.7];
bondEnergyMatrixRNARNA(2,:,:) = [-9.1, -5.6, -5.6, -10.44; -5.7, -3.4, -13.39, -2.7; -8.2, -10.64, -9.2, -5.61; -10.48, -5.3, -12.11, -8.6];
bondEnergyMatrixRNARNA(3,:,:) = [-5.2, -4, -5.6, -12.44; -7.2, 0.5, -14.88, -4.2; -7.1, -13.39, -6.2, -8.33; -11.4, -0.3, -12.59, -5];
bondEnergyMatrixRNARNA(4,:,:) = [-4, -6.3, -8.9, -7.69; -4.3, -5.1, -12.44, -1.8; -3.8, -10.44, -8.9, -6.99; -6.82, -1.4, -12.83, 1.4];
bondEnergyMatrixRNARNA(5,:,:) = [3.4, 2, -3.5, -12.83; -2.3, 6, -12.59, -0.3; 0.6, -12.11, -3.5, -13.47; -8.81, 4.6, -14.59, -1.7];
bondEnergyMatrixRNARNA(6,:,:) = [-4.8, -6.3, -8.9, -6.99; -4.3, -5.1, -8.33, -1.8; -3.1, -5.61, 1.5, -9.26; -3.21, -1.4, -13.47, 1.4];

bondEntropyMatrixRNARNA(1,:,:) = [-10.2, 9.6, -8.7, -19; -5.3, 21.6, -29.5, 1.5; -7.3, -27.1, -8.7, -8.6; -26.7, 17.4, -24, -2.7];
bondEntropyMatrixRNARNA(2,:,:) = [-24.5, -13.5, -13.4, -26.9; -15.2, -7.6, -32.7, -6.3; -21.8, -26.7, -24.6, -13.5; -27.1, -12.6, -32.2, -23.9];
bondEntropyMatrixRNARNA(3,:,:) = [-13.2, -8.2, -13.9, -32.5; -19.6, 3.9, -36.9, -12.2; -17.8, -32.7, -15.1, -21.9; -29.5, -2.1, -32.5, -14];
bondEntropyMatrixRNARNA(4,:,:) = [-9.7, -17.1, -25.2, -20.5; -11.6, -14.6, -32.5, -4.2; -8.5, -26.9, -25, -19.3; -19, -2.5, -37.3, 6];
bondEntropyMatrixRNARNA(5,:,:) = [10, 9.6, -8.7, -37.3; -5.3, 21.6, -32.5, 1.5; 0, -32.2, -8.7, -44.9; -24, 17.4, -51.2, -2.7];
bondEntropyMatrixRNARNA(6,:,:) = [-12.1, -17.7, -25.2, -19.3; -11.6, -14.6, -21.9, -4.2; -11.2, -13.5, 2.1, -30.8; -8.6, -2.5, -44.9, 6];

bondEntropyMatrixRNARNA = bondEntropyMatrixRNARNA./1000; %to convert from eu to kcal/(mol*K)

bondFreeEnergyMatrixRNARNA(1,:,:) = [-0.8 -1.0 -0.8 -0.93; -0.6 -0.7 -2.24 -0.7; -0.8 -2.08 -0.8 -0.55; -1.10 -0.8 -1.36 -0.8];
bondFreeEnergyMatrixRNARNA(2,:,:) = [-1.5 -1.5 -1.4 -2.11; -1.0 -1.1 -3.26 -0.8; -1.4 -2.36 -1.6 -1.41; -2.08 -1.4 -2.11 -1.2];
bondFreeEnergyMatrixRNARNA(3,:,:) = [-1.1 -1.5 -1.3 -2.35; -1.1 -0.7 -3.42 -0.5; -1.6 -3.26 -1.4 -1.53; -2.24 -1.0 -2.51 -0.7];
bondFreeEnergyMatrixRNARNA(4,:,:) = [-1.0 -0.8 -1.1 -1.33; -0.7 -0.6 -2.35 -0.5; -1.1 -2.11 -1.2 -1.00; -0.93 -0.6 -1.27 -0.5];
bondFreeEnergyMatrixRNARNA(5,:,:) = [0.3 -1.0 -0.8 -1.27; -0.6 -0.7 -2.51 -0.7; 0.6 -2.11 -0.8 -0.500; -1.36 -0.8 1.29 -0.8]; 
bondFreeEnergyMatrixRNARNA(6,:,:) = [-1.0 -0.8 -1.1 -1.00; -0.7 -0.6 -1.53 -0.5; 0.5 -1.41 0.8 0.30; -0.55 -0.6 -0.500 -0.5];
%the -0.500 was actually measured at +0.47 but the authors claim -0.5 is a better estimate.


bondFreeEnergyMatrixDNADNA = zeros(4,4,4); %matrix of stacking energies (deltaG). Units are kcal/mol.
bondEnergyMatrixDNADNA = zeros(4,4,4); %matrix of enthalpies (deltaH). Units are kcal/mol.
%First index tells you if the first bp of the set is AT (1) CG (2) GC (3)
%or TA (4)
%Second index tells you if the 3' ntd of the second bp is A (1) C (2) G(3)
%or T(4)
%Third index tells you if the 5' ntd of the second bp is A (1) C (2) G(3)
%or T(4) 
%Data is from From: Thermodynamics and NMR of Internal G·T Mismatches in DNA and other papers by Allawi and 
%various SantaLucia publications (cited as 28-32 in Mfold web server for nucleic acid folding and hybridization prediction)
bondFreeEnergyMatrixDNADNA(1,:,:) = [0.61 0.88 0.14 -1.0; 0.77 1.33 -1.44 0.64; 0.02 -1.28 -0.13 0.71; -0.88 0.73 0.07 0.69];
bondFreeEnergyMatrixDNADNA(2,:,:) = [0.43 0.75 0.03 -1.45; 0.79 0.7 -1.84 0.62; 0.11 -2.17 -0.11 -0.47; -1.28 0.4 -0.32 -0.21];
bondFreeEnergyMatrixDNADNA(3,:,:) = [0.17 0.81 -0.25 -1.3; 0.47 0.79 -2.24 0.62; -0.52 -1.84 -1.11 0.08; -1.44 0.98 -0.59 0.45];
bondFreeEnergyMatrixDNADNA(4,:,:) = [0.69 0.92 0.42 -0.58; 1.33 1.05 -1.3 0.97; 0.74 -1.45 0.44 0.43; -1.0 0.75 0.34 0.68];

bondEnergyMatrixDNADNA(1,:,:) = [1.2 2.3 -0.6 -7.9; 5.3 0.0 -8.4 0.7; -0.7 -7.8 -3.1 1.0; -7.2 -1.2 -2.5 -2.7];
bondEnergyMatrixDNADNA(2,:,:) = [-0.9 1.9 -0.7 -8.5; 0.6 -1.5 -8 -0.8; -4 -10.6 -4.9 -4.1; -7.8 -1.5 -2.8 -5];
bondEnergyMatrixDNADNA(3,:,:) = [-2.9 5.2 -0.6 -8.2; -0.7 3.6 -9.8 2.3; 0.5 -8 -6 3.3; -8.4 5.2 -4.4 -2.2];
bondEnergyMatrixDNADNA(4,:,:) = [4.7 3.4 0.7 -7.2; 7.6 6.1 -8.2 1.2; 3 -8.5 1.6 -0.1; -7.9 1.0 -1.3 0.2];

bondEntropyMatrixDNADNA = -(bondFreeEnergyMatrixDNADNA - bondEnergyMatrixDNADNA)./(273.15+37);


bondFreeEnergyMatrixRNADNA = zeros(8,4,4); %matrix of stacking energies (deltaG). Units are kcal/mol.
bondEnergyMatrixRNADNA = zeros(8,4,4); %matrix of enthalpies (deltaH). Units are kcal/mol.
%First index tells you if the set is (putting RNA first) 5'AX3'/3'TY5' (1)
%5CX3/3GY5 (2) 5GX3/3CY5 (3) 5UX3/3AY5 (4) 5XA3/3YT5 (5) 5XC3/3YG5 (6)
%5XG3/3YC5 (7) or 5XU3/3YA5.
%Second index tells you if X is A (1) C (2) G (3) or U (4)
%Third index tells you if Y is A (1) C (2) G (3) or T (4)
%Data is from From: Thermodynamic Parameters To Predict Stability of RNA/DNA Hybrid Duplexes & 
%Thermodynamic contributions of single internal rA·dA, rC·dC, rG·dG and rU·dT mismatches in RNA/DNA duplexes
%100 is put in place of elements for which there aren't published parameters
bondFreeEnergyMatrixRNADNA(1,:,:) = [1.07 100 100 -1.0; 100 1.64 -2.1 100; 100 -1.8 0.31 100; -0.9 100 100 0.63];
bondFreeEnergyMatrixRNADNA(2,:,:) = [0.90 100 100 -0.9; 100 1.04 -2.1 100; 100 -1.7 0.14 100; -0.9 100 100 0.49];
bondFreeEnergyMatrixRNADNA(3,:,:) = [0.51 100 100 -1.3; 100 0.96 -2.7 100; 100 -2.9 -0.58 100; -1.1 100 100 0.18];
bondFreeEnergyMatrixRNADNA(4,:,:) = [1.13 100 100 -0.6; 100 1.15 -1.5 100; 100 -1.6 0.44 100; -0.2 100 100 1.07];
bondFreeEnergyMatrixRNADNA(5,:,:) = [1.36 100 100 -1.0; 100 1.70 -0.9 100; 100 -1.3 0.50 100; -0.6 100 100 1.21];
bondFreeEnergyMatrixRNADNA(6,:,:) = [0.19 100 100 -2.1; 100 0.73 -2.1 100; 100 -2.7 -0.83 100; -1.5 100 100 -0.02];
bondFreeEnergyMatrixRNADNA(7,:,:) = [0.21 100 100 -1.8; 100 0.46 -1.7 100; 100 -2.9 -0.33 100; -1.6 100 100 0.14];
bondFreeEnergyMatrixRNADNA(8,:,:) = [1.85 100 100 -0.9; 100 1.88 -0.9 100; 100 -1.1 0.97 100; -0.2 100 100 1.03];

bondEnergyMatrixRNADNA(1,:,:) = [-4.3 100 100 -7.8; 100 -8.8 -5.9 100; 100 -9.1 -3.3 100; -8.3 100 100 0.6];
bondEnergyMatrixRNADNA(2,:,:) = [5.5 100 100 -9.0; 100 10.5 -9.3 100; 100 -16.3 -8.9 100; -7.0 100 100 -0.4];
bondEnergyMatrixRNADNA(3,:,:) = [-1.9 100 100 -5.5; 100 -0.1 -8.0 100; 100 -12.8 -8.0 100; -7.8 100 100 -11.6];
bondEnergyMatrixRNADNA(4,:,:) = [-1.7 100 100 -7.8; 100 -3.3 -8.6 100; 100 -10.4 -5.8 100; -11.5 100 100 -2.2];
bondEnergyMatrixRNADNA(5,:,:) = [3.0 100 100 -7.8; 100 -0.3 -9.0 100; 100 -5.5 1.1 100; -7.8 100 100 -3.3];
bondEnergyMatrixRNADNA(6,:,:) = [-6.0 100 100 -5.9; 100 9.3 -9.3 100; 100 -8.0 -7.0 100; -8.6 100 100 0.1];
bondEnergyMatrixRNADNA(7,:,:) = [-10.5 100 100 -9.1; 100 -11.5 -16.3 100; 100 -12.8 -16.5 100; -10.4 100 100 -13.4];
bondEnergyMatrixRNADNA(8,:,:) = [5.6 100 100 -8.3; 100 0.8 -7.0 100; 100 -7.8 -3.7 100; -11.5 100 100 3.0];

bondEntropyMatrixRNADNA = -(bondFreeEnergyMatrixRNADNA - bondEnergyMatrixRNADNA)./(273.15+37);

for i = 1:size(bondEntropyMatrixRNADNA,1)
    for j = 1:size(bondEntropyMatrixRNADNA,2)
        for k = 1:size(bondEntropyMatrixRNADNA,3)
            if bondEntropyMatrixRNADNA(i,j,k) == 0
                if i <=4
                    bondFreeEnergyMatrixRNADNA(i,j,k) = (bondFreeEnergyMatrixDNADNA(i,j,k) + bondFreeEnergyMatrixRNARNA(i,j,k)) / 2;
                    bondEnergyMatrixRNADNA(i,j,k) = (bondEnergyMatrixDNADNA(i,j,k) + bondEnergyMatrixRNARNA(i,j,k)) / 2;
                    bondEntropyMatrixRNADNA(i,j,k) = (bondEntropyMatrixDNADNA(i,j,k) + bondEntropyMatrixRNARNA(i,j,k)) / 2;
                else
                    bondEnergyMatrixRNADNA(i,j,k) = (bondEnergyMatrixDNADNA(9-i,k,j) + bondEnergyMatrixRNARNA(9-i,k,j)) / 2;
                    bondFreeEnergyMatrixRNADNA(i,j,k) = (bondFreeEnergyMatrixDNADNA(9-i,k,j) + bondFreeEnergyMatrixRNARNA(9-i,k,j)) / 2;
                    bondEntropyMatrixRNADNA(i,j,k) = (bondEntropyMatrixDNADNA(9-i,k,j) + bondEntropyMatrixRNARNA(9-i,k,j)) / 2;
                end
            end
        end
    end
end



%Dangling ends also have an associated free energy. 
%first index tells you if the paired ntd is an A(1) C(2) G(3) or U(4);
%by paired ntd, we mean the paired ntd either 1 before or 1 after the
%dangling ntd
%second index tells you if the dangling ntd is A(1) C(2) G(3) or U(4). Data
%taken from Xia, Mathews, Turner "Thermodynamics of RNA secondary structure 
%formation" in book by Soll, Nishmura, Moore.
%It is the same as the data from Serra & Turner Table 2 (second half of the table).
%GU base pairs are treated here as if they were AU (so G is replaced with A
%if the base pair is AU).

%For a 5' dangling end:
dangling5FreeEnergyRNARNA = [-0.3 -0.3 -0.4 -0.2; -0.5 -0.3 -0.2 -0.1; -0.2 -0.3 -0.0 -0.0; -0.3 -0.1 -0.2 -0.2];
dangling5EnergyRNARNA = [1.6, 2.2, 0.7, 3.1; -2.4, 3.3, 0.8, -1.4; -1.6, 0.7, -4.6, -0.4; -0.5, 6.9, 0.6, 0.6];
dangling5EntropyRNARNA = [6.1, 7.9, 3.4, 10.6; -6.0, 11.8, 3.4, -4.3; -4.5, 3.1, -14.8, -1.2; -0.7, 22.8, 2.7, 2.7];
dangling5EntropyRNARNA = dangling5EntropyRNARNA./1000; %to convert from eu to kcal/(mol*K)
%For a 3' dangling end:
dangling3FreeEnergyRNARNA = [-0.8 -0.5 -0.8 -0.6; -1.7 -0.8 -1.7 -1.2; -1.1 -0.4 -1.3 -0.6; -0.7 -0.1 -0.7 -0.1];
dangling3EnergyRNARNA = [-4.9 -0.9 -5.5 -2.3; -9.0 -4.1 -8.6 -7.5; -7.4 -2.8 -6.4 -3.6; -5.7 -0.7 -5.8 -2.2];
dangling3EntropyRNARNA = [-13.2, -1.2, -15.0, -5.4; -23.4, -10.7, -22.2, -20.4; -20.0, -7.9, -16.6, -9.7; -16.4, -1.8, -16.4, -6.8];
dangling3EntropyRNARNA = dangling3EntropyRNARNA./1000; %to convert from eu to kcal/(mol*K)


%DNA dangling end parameters taken from "Thermodynamic parameters for DNA
%sequences with dangling ends". All energies in kcal/mol, so entropies are
%in kcal/(mol*K).
%For a 5' dangling end:
dangling5FreeEnergyDNADNA = [-0.51 -0.96 -0.58 -0.5; -0.42 -0.52 -0.34 -0.02; -0.62 -0.72 -0.56 0.48; -0.71 -0.58 -0.61 -0.10]';
dangling5EnergyDNADNA = [0.2 -0.62 -3.7 -2.9; 0.6 -4.4 -4.0 -4.1; -1.1 -5.1 -3.9 -4.2; -6.9 -4.0 -4.9 -0.2]'; %transpose because I copied this from excel doc which had matrix transposed
dangling5EntropyDNADNA = -(dangling5FreeEnergyDNADNA - dangling5EnergyDNADNA)./(273.15+37);
%For a 3' dangling end:
dangling3FreeEnergyDNADNA = [-0.12 -0.82 -0.92 -0.48; 0.28 -0.31 -0.23 -0.19; -0.01 -0.01 -0.44 -0.50; 0.13 -0.52 -0.35 -0.29]';
dangling3EnergyDNADNA = [-0.5 -5.9 -2.1 -0.7; 4.7 -2.6 -0.2 4.4; -4.1 -3.2 -3.9 -1.6; -3.8 -5.2 -4.4 2.9]';
dangling3EntropyDNADNA = -(dangling3FreeEnergyDNADNA - dangling3EnergyDNADNA)./(273.15+37);


dangling5FreeEnergy = [dangling5FreeEnergyRNARNA;dangling5FreeEnergyDNADNA];
dangling5Energy = [dangling5EnergyRNARNA;dangling5EnergyDNADNA];
dangling5Entropy = [dangling5EntropyRNARNA;dangling5EntropyDNADNA];

dangling3FreeEnergy = [dangling3FreeEnergyRNARNA;dangling3FreeEnergyDNADNA];
dangling3Energy = [dangling3EnergyRNARNA;dangling3EnergyDNADNA];
dangling3Entropy = [dangling3EntropyRNARNA;dangling3EntropyDNADNA];

%first index tells you if the paired ntd is an RNA A(1) RNA C(2) RNA G(3) RNA U(4) DNA A(5) DNA C(6) DNAG(7) or DNA U(8);
%second index tells you if the dangling ntd is A(1) C(2) G(3) or U/T(4).
%This doesn't deal well with RNA/DNA hybrids (it uses the parameters of whichever type of ntd is dangling) because I couldn't find experimental estimates of the
%relevant parameters.

%from NNDB In the case of mismatch-mediated coaxial stacking, there are two adjacent stacks. 
%The stack of the mismatch on the adjacent helix, where there is no break in the backbone, 
%is approximated using the terminal mismatch parameters. The second stack is the stack of the 
%mismatch on the second helix, where the backbone is not continuous. This stack is approximated 
%using a sequence-independent term of ?2.1 kcal/mol for folding free energy change and 
%?8.46 ± 2.75 kcal/mol for enthalpy change. If the "mismatch" mediating the coaxial stack 
%could form a Watson-Crick or GU pair, a bonuses of ?0.4 or ?0.2 kcal/mol, respectively, 
%are applied to both free energy and enthalpy changes.
MMCSEnergy = -8.46; %kcal/mol
MMCSEntropy = (-2.1+8.46)/(-310.15); 
extraMMCSEnergyWC = -0.4; extraMMCSEnergyGU = -0.2; 


%listOfBondedNtds = zeros(1,numBP); %list of ntds in the structure. Order doesn't matter except for considering bulge loops
%order it so that listOfBondedNtds(i) is bonded to listOfBondedNtds(i+length(listOfBondedNtds)/2)
%bondedCounter = 1;
listOfFirstNtds = zeros(1,floor(numBP/2));
listOfSecondNtds = zeros(1,floor(numBP/2));
halfBondedCounter = 1;
for i = 1:length(structure)
    %listOfBondedNtds(bondedCounter:bondedCounter + length(structure{i})-1) = structure{i};
    listOfFirstNtds(halfBondedCounter:halfBondedCounter + length(structure{i})/2-1) = structure{i}(1:length(structure{i})/2);
    listOfSecondNtds(halfBondedCounter:halfBondedCounter + length(structure{i})/2-1) = structure{i}(length(structure{i})/2+1:end);
    %bondedCounter = bondedCounter + length(structure{i});
    halfBondedCounter = halfBondedCounter + length(structure{i})/2;
end
listOfFirstNtds(halfBondedCounter:end) = [];
listOfSecondNtds(halfBondedCounter:end) = [];
listOfBondedNtds = [listOfFirstNtds,listOfSecondNtds];
%listOfBondedNtds(bondedCounter:end) = [];


numStems = length(structure);
numNodes = numStems*2; %even if a stem has only one bp, you can have a terminal mismatch on either side of it so we count it as two nodes here.
TMs = cell(numNodes,1); %for each node, what possibilities for terminal mismatches it has
FCSs = cell(numNodes,1); %flush coaxial stacks
MMCSs = cell(numNodes,1); %mismatch-mediated coaxial stacks
DEs = cell(numNodes,1); %dangling ends
nodeNtds = zeros(numNodes,2);
TMsFE = zeros(numNodes,1);
FCSsFE = zeros(numNodes,1);
MMCSsFE = zeros(numNodes,1);
DEsFE = zeros(numNodes,1);

%%
for i = 1:length(structure)
    numBonds = length(structure{i})/2;
    isParallel = false;
    nodeNtds(i,:) = [structure{i}(numBonds+1),structure{i}(numBonds*2)];
    nodeNtds(i+numStems,:) = [structure{i}(1),structure{i}(numBonds)];
    
    %first, find if the stem is parallel or antiparallel or not so as to know how
    %nearest neighbor interactions will work (if it's k+1 bonded to l+1 or
    %l-1)
    
    if numBonds == 1 %figuring out if it should be treated as parallel or antiparallel is tricky
        k = structure{i}(1);
        l = structure{i}(2);
        %find whether next bonded ntd in the sequence after k is bonded to
        %a ntd before or after l.
        distFromK = listOfBondedNtds-k;
        distFromK(distFromK<=0) = inf;
        nextNtd = min(distFromK)+k; %there has to be another ntd after k because l>k
        for j = 1:length(structure)
            if any(structure{j} == nextNtd) %ismember(nextNtd,structure{j})
                nextNtdIndex = find(structure{j}==nextNtd);
                if nextNtdIndex > length(structure{j})/2
                    nextNtdBonded = structure{j}(nextNtdIndex - length(structure{j})/2);
                else
                    nextNtdBonded = structure{j}(nextNtdIndex + length(structure{j})/2);
                end
            end
        end
        
        if nextNtdBonded <l
            isParallel = false;
        else
            isParallel = true;
        end
    elseif structure{i}(2+numBonds)-structure{i}(1+numBonds)==1
        isParallel = true;
    elseif structure{i}(2+numBonds)-structure{i}(1+numBonds)==-1
        isParallel = false;
    end
    
    
    %Calculate free energies of bps -- we don't include terminal mismatches at this step.
    for j = 1:numBonds-1 
        k = structure{i}(j);
        l = structure{i}(j+numBonds);
        firstNtd = sequenceInNumbers(k); %what is the 5' ntd?
        firstBP = sequenceInNumbers(l); %what is it bonded to?
        secondNtd = sequenceInNumbers(k+1); %what is the 3' ntd? %this must exist by the rules we've set up which ntds can form bps.
        secondBP = sequenceInNumbers(l+(2*isParallel-1)); %what is it bonded to? %this must exist by the rules we've set up which ntds can form bps.
            %this is l-1 for antiparallel and l+1 for parallel.
            
        
        if firstNtd <=4 && firstBP <= 4 && secondNtd <=4 && secondBP <=4
            ntdType = 0; %RNARNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNARNA,bondEntropyMatrixRNARNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
        elseif firstNtd <=4 && firstBP > 4 && secondNtd <=4 && secondBP > 4
            ntdType = 1; %RNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
        elseif firstNtd > 4 && firstBP <=4 && ~isParallel && secondNtd > 4 && secondBP <=4
            ntdType = 1; %RNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,secondBP,secondNtd,firstBP,firstNtd,ntdType);
        elseif firstNtd > 4 && firstBP <=4 && isParallel && secondNtd > 4 && secondBP <=4
            ntdType = 1; %RNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstBP,firstNtd,secondBP,secondNtd,ntdType);
        elseif firstNtd > 4 && firstBP > 4 && secondNtd > 4 && secondBP > 4
            ntdType = 2; %DNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixDNADNA,bondEntropyMatrixDNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
        end
        bondEnergyLocal = bondEnergyLocal + bondEnergy;
        bondEntropyLocal = bondEntropyLocal + bondEntropy;
    end
    

    %we treat bulge and internal loops the same as multiloops -- therefore,
    %we allow them to have coaxial stacking (if applicable) as well as
    %terminal mismatches and dangling ends (again, as applicable).
    j = numBonds;
    
    %if both terminal ntds are unbound, we can have a terminal mismatch
    if ~any(listOfBondedNtds==structure{i}(j)+1) && ~any(listOfBondedNtds == structure{i}(j+numBonds)+(2*isParallel-1))...
            && any(realNtds == structure{i}(j+numBonds)+(2*isParallel-1)) %ismember(structure{i}(j+numBonds)+(2*isParallel-1),realNtds) 
        k = structure{i}(j);
        l = structure{i}(j+numBonds);
        firstNtd = sequenceInNumbers(k); %what is the 5' ntd?
        firstBP = sequenceInNumbers(l); %what is it bonded to?
        secondNtd = sequenceInNumbers(k+1); %what is the 3' ntd? %this must exist by the rules we've set up which ntds can form bps.
        secondBP = sequenceInNumbers(l+(2*isParallel-1)); %what is it bonded to? %this must exist by the rules we've set up which ntds can form bps.
            %this is l-1 for antiparallel and l+1 for parallel.
        
        %if secondNtd and secondBP could bind but aren't bound in this
        %structure, Lu, Turner, Mathews (NAR, 2006) say that we should
        %treat them as an AC pair (where A replaces the purine and C the
        %pyrimidine)
        realSecondNtd = secondNtd; realSecondBP = secondBP;
        if (secondNtd == 1 && secondBP == 4) || (secondNtd == 4 && secondBP == 1) ||...
                (secondNtd == 2 && secondBP == 3) || (secondNtd == 3 && secondBP == 2) ||...
                (secondNtd == 3 && secondBP == 4) || (secondNtd == 4 && secondBP == 3)
            if secondNtd == 1 || secondNtd == 3 %if secondNtd is a purine
                secondNtd = 1; %the purine is replaced with A 
                secondBP = 2; %and the pyrimidine with C
            else
                secondNtd = 2;
                secondBP = 1;
            end
        end
        
        if firstNtd <=4 && firstBP <= 4 && secondNtd <=4 && secondBP <=4
            ntdType = 0; %RNARNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNARNA,bondEntropyMatrixRNARNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
        elseif firstNtd <=4 && firstBP > 4 && secondNtd <=4 && secondBP > 4
            ntdType = 1; %RNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
        elseif firstNtd > 4 && firstBP <=4 && ~isParallel && secondNtd > 4 && secondBP <=4
            ntdType = 1; %RNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,secondBP,secondNtd,firstBP,firstNtd,ntdType);
        elseif firstNtd > 4 && firstBP <=4 && isParallel && secondNtd > 4 && secondBP <=4
            ntdType = 1; %RNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstBP,firstNtd,secondBP,secondNtd,ntdType);
        elseif firstNtd > 4 && firstBP > 4 && secondNtd > 4 && secondBP > 4
            ntdType = 2; %DNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixDNADNA,bondEntropyMatrixDNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
        end
        
        TMs{i}{1+length(TMs{i})} = {[k,l,k+1,l+(2*isParallel-1)],bondEnergy,bondEntropy};
        TMsFE(i) = min(TMsFE(i),TMs{i}{length(TMs{i})}{2}-T*TMs{i}{length(TMs{i})}{3});
        
        secondNtd = realSecondNtd; secondBP = realSecondBP;
        %we also might have a mismatch-mediated coaxial stack (MMCS)
        if any(listOfBondedNtds==structure{i}(j)+2) || any(listOfBondedNtds == structure{i}(j+numBonds)+2*(2*isParallel-1)) %if either of the next ntds over are bound, 
             if any(listOfBondedNtds==structure{i}(j)+2) && any(listOfBondedNtds == structure{i}(j+numBonds)+2*(2*isParallel-1)) ...
                     && abs(find(listOfBondedNtds==structure{i}(j)+2)-find(listOfBondedNtds == structure{i}(j+numBonds)+2*(2*isParallel-1))) == length(listOfBondedNtds)/2 
                        %but if they're bound to each other, that's a 1x1 internal loop
                %don't do anything if we have a 1x1 internal loop, since we
                %treat that using two terminal mismatches.
             else
                 
                 MMCScouldBindWC = false; MMCScouldBindGU = false;
                 if (secondNtd == 1 && secondBP == 4) || (secondNtd == 4 && secondBP == 1) ||...
                        (secondNtd == 2 && secondBP == 3) || (secondNtd == 3 && secondBP == 2)
                    MMCScouldBindWC = true;
                 elseif (secondNtd == 3 && secondBP == 4) || (secondNtd == 4 && secondBP == 3)
                     MMCScouldBindGU = true;
                 end
                 
                 if any(listOfBondedNtds==structure{i}(j)+2)
                     kPlusTwoBPIndex = max(0,find(listOfBondedNtds==structure{i}(j)+2)-length(listOfBondedNtds)/2); %could be either + or - length(listOfBondedNtds)/2
                     if kPlusTwoBPIndex == 0
                         kPlusTwoBP = listOfBondedNtds(find(listOfBondedNtds==structure{i}(j)+2)+length(listOfBondedNtds)/2);
                     else
                         kPlusTwoBP = listOfBondedNtds(kPlusTwoBPIndex);
                     end
                     MMCSs{i}{1+length(MMCSs{i})} = {[k,l,k+1,l+(2*isParallel-1),k+2,kPlusTwoBP],...
                         bondEnergy+MMCSEnergy + extraMMCSEnergyWC*MMCScouldBindWC + extraMMCSEnergyGU*MMCScouldBindGU,bondEntropy+MMCSEntropy};
                     MMCSsFE(i) = min(MMCSsFE(i),MMCSs{i}{length(MMCSs{i})}{2}-T*MMCSs{i}{length(MMCSs{i})}{3});
                 end
                 if any(listOfBondedNtds == structure{i}(j+numBonds)+2*(2*isParallel-1))
                     lPlusTwoBPIndex = max(0,find(listOfBondedNtds==structure{i}(j+numBonds)+2*(2*isParallel-1))-length(listOfBondedNtds)/2); %could be either + or - length(listOfBondedNtds)/2
                     if lPlusTwoBPIndex == 0
                         lPlusTwoBP = listOfBondedNtds(find(listOfBondedNtds==structure{i}(j+numBonds)+2*(2*isParallel-1))+length(listOfBondedNtds)/2);
                     else
                         lPlusTwoBP = listOfBondedNtds(lPlusTwoBPIndex);
                     end
                     MMCSs{i}{1+length(MMCSs{i})} = {[k,l,k+1,l+(2*isParallel-1),l+2*(2*isParallel-1),lPlusTwoBP],...
                         bondEnergy+MMCSEnergy  + extraMMCSEnergyWC*MMCScouldBindWC + extraMMCSEnergyGU*MMCScouldBindGU,bondEntropy+MMCSEntropy};
                     MMCSsFE(i) = min(MMCSsFE(i),MMCSs{i}{length(MMCSs{i})}{2}-T*MMCSs{i}{length(MMCSs{i})}{3});
                 end
             end
        end
    end  
        
    %if only one terminal ntd is unbound, we might have a dangling end, or
    %we might have a flush coaxial stack (we don't consider the ends of the
    %sequence here, only later).
    %We also consider the possibility of dangling ends even if both
    %terminal ntds are unbound since we might have a multiloop. 
    if ~any(listOfBondedNtds == structure{i}(j+numBonds)+(2*isParallel-1))...% && any(listOfBondedNtds==structure{i}(j)+1) ...
            && any(realNtds == structure{i}(j+numBonds)+(2*isParallel-1)) %ismember(structure{i}(j+numBonds)+(2*isParallel-1),realNtds)
        %this describes a 5' dangling end if it's not parallel since the 3' end of the first
            %strand is bound and the 5' end of the second strand is unbound   

        firstBoundNtd = sequenceInNumbers(structure{i}(j+numBonds));
        danglingNtd = sequenceInNumbers(structure{i}(j+numBonds)+(2*isParallel-1));
        if (firstBoundNtd == 4 && sequenceInNumbers(structure{i}(j))==3) ||...
                (firstBoundNtd == 3 && sequenceInNumbers(structure{i}(j))==4) %if the pair is a GU base pair,
            %then treat it as AU for the purposes of calculating dangling
            %ends (from Lu, Turner, Mathews (NAR, 2006))
            if firstBoundNtd == 3
                firstBoundNtd = 1; %i.e. change the G to an A
            end
        end
        if danglingNtd > 4
            danglingNtd = danglingNtd - 4;
        end
        if (firstBoundNtd == 0 || danglingNtd == 0)
            %if one of the sequence elements is unknown, don't do anything
        else
            if ~isParallel
                DEs{i}{1+length(DEs{i})} = {[structure{i}(j+numBonds),structure{i}(j+numBonds)+(2*isParallel-1),structure{i}(j)],...
                    dangling5Energy(firstBoundNtd,danglingNtd),dangling5Entropy(firstBoundNtd,danglingNtd)};
                DEsFE(i) = min(DEsFE(i),DEs{i}{length(DEs{i})}{2}-T*DEs{i}{length(DEs{i})}{3});
            else
                DEs{i}{1+length(DEs{i})} = {[structure{i}(j+numBonds),structure{i}(j+numBonds)+(2*isParallel-1),structure{i}(j)],...
                    dangling3Energy(firstBoundNtd,danglingNtd),dangling3Entropy(firstBoundNtd,danglingNtd)};
                DEsFE(i) = min(DEsFE(i),DEs{i}{length(DEs{i})}{2}-T*DEs{i}{length(DEs{i})}{3});
            end
        end
    end
    %now consider the flush coaxial stack
    if any(listOfBondedNtds==structure{i}(j)+1)
        
        isBulgeLoopOfLengthGreater1 = false;
        isBulgeLoopOfLength1 = false;
        if ~isParallel %can't have bulge loop if it's parallel
            bondedToNextIndex = find(listOfBondedNtds==structure{i}(j)+1);
            if bondedToNextIndex - length(listOfBondedNtds)/2 <= 0
                bondedToNextIndex = bondedToNextIndex + length(listOfBondedNtds)/2;
            else
                bondedToNextIndex = bondedToNextIndex - length(listOfBondedNtds)/2;
            end
            bondedToNext = listOfBondedNtds(bondedToNextIndex);
            interveningNtds = structure{i}(j+numBonds)+(2*isParallel-1):bondedToNext+1;
            if isempty(interveningNtds)
                interveningNtds = structure{i}(j+numBonds)+(2*isParallel-1):-1:bondedToNext+1;
            end
            if ~isempty(interveningNtds) && all(~any(interveningNtds == listOfBondedNtds')) %if all the intervening ntds are unbound
                isBulgeLoopOfLengthGreater1 = true;
                if length(interveningNtds) == 1
                    isBulgeLoopOfLength1 = true;
                    isBulgeLoopOfLengthGreater1 = false;
                end
            end
        end
        
        if ~isBulgeLoopOfLengthGreater1
            k = structure{i}(j);
            l = structure{i}(j+numBonds);
            firstNtd = sequenceInNumbers(k); %what is the 5' ntd?
            firstBP = sequenceInNumbers(l); %what is it bonded to?
            secondNtd = sequenceInNumbers(k+1); %what is the 3' ntd? %this must exist by the rules we've set up which ntds can form bps.
            ThreePrimeIsBondedToIndex = find(listOfBondedNtds==k+1);
            if ThreePrimeIsBondedToIndex - length(listOfBondedNtds)/2 <= 0
                ThreePrimeIsBondedTo = listOfBondedNtds(ThreePrimeIsBondedToIndex + length(listOfBondedNtds)/2);
            else
                ThreePrimeIsBondedTo = listOfBondedNtds(ThreePrimeIsBondedToIndex - length(listOfBondedNtds)/2);
            end
            secondBP = sequenceInNumbers(ThreePrimeIsBondedTo); %what is it bonded to? %this must exist by the rules we've set up which ntds can form bps.

            if firstNtd <=4 && firstBP <= 4 && secondNtd <=4 && secondBP <=4
                ntdType = 0; %RNARNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNARNA,bondEntropyMatrixRNARNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            elseif firstNtd <=4 && firstBP > 4 && secondNtd <=4 && secondBP > 4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            elseif firstNtd > 4 && firstBP <=4 && ~isParallel && secondNtd > 4 && secondBP <=4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,secondBP,secondNtd,firstBP,firstNtd,ntdType);
            elseif firstNtd > 4 && firstBP <=4 && isParallel && secondNtd > 4 && secondBP <=4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstBP,firstNtd,secondBP,secondNtd,ntdType);
            elseif firstNtd > 4 && firstBP > 4 && secondNtd > 4 && secondBP > 4
                ntdType = 2; %DNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixDNADNA,bondEntropyMatrixDNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            end

            FCSs{i}{1+length(FCSs{i})} = {[k,l,k+1,ThreePrimeIsBondedTo],bondEnergy,bondEntropy};
            FCSsFE(i) = min(FCSsFE(i),FCSs{i}{length(FCSs{i})}{2}-T*FCSs{i}{length(FCSs{i})}{3});
        end
    end

    
    if ~any(listOfBondedNtds==structure{i}(j)+1)... %&& any(listOfBondedNtds == structure{i}(j+numBonds)+(2*isParallel-1))...
            && any(realNtds == structure{i}(j)+1) %ismember(structure{i}(j+numBonds)+(2*isParallel-1),realNtds)
        %this describes a 3' dangling end if it's not parallel since the 5'
        %end of the second strand is bound and the 3' end of the first strand is unbound
        firstBoundNtd = sequenceInNumbers(structure{i}(j));
        danglingNtd = sequenceInNumbers(structure{i}(j)+1);
        if (firstBoundNtd == 4 && sequenceInNumbers(structure{i}(j+numBonds))==3) ||...
            (firstBoundNtd == 3 && sequenceInNumbers(structure{i}(j+numBonds))==4) %if the pair is a GU base pair,
            %then treat it as AU for the purposes of calculating dangling
            %ends (from Lu, Turner, Mathews (NAR, 2006))
            if firstBoundNtd == 3
                firstBoundNtd = 1; %i.e. change the G to an A
            end
        end
        if danglingNtd > 4
            danglingNtd = danglingNtd - 4;
        end
        if (firstBoundNtd == 0 || danglingNtd == 0)
            %if one of the sequence elements is unknown, don't do anything
        else
            if ~isParallel
                DEs{i}{1+length(DEs{i})} = {[structure{i}(j+numBonds),structure{i}(j)+1,structure{i}(j)],dangling3Energy(firstBoundNtd,danglingNtd),dangling3Entropy(firstBoundNtd,danglingNtd)};
                DEsFE(i) = min(DEsFE(i),DEs{i}{length(DEs{i})}{2}-T*DEs{i}{length(DEs{i})}{3});
            else
                DEs{i}{1+length(DEs{i})} = {[structure{i}(j+numBonds),structure{i}(j)+1,structure{i}(j)],dangling5Energy(firstBoundNtd,danglingNtd),dangling5Entropy(firstBoundNtd,danglingNtd)};
                DEsFE(i) = min(DEsFE(i),DEs{i}{length(DEs{i})}{2}-T*DEs{i}{length(DEs{i})}{3});
            end
        end
    end
    %now consider the flush coaxial stack
    if any(listOfBondedNtds == structure{i}(j+numBonds)+(2*isParallel-1))
        isBulgeLoopOfLengthGreater1 = false;
        isBulgeLoopOfLength1 = false;
        if ~isParallel %can't have bulge loop if it's parallel
            bondedToNextIndex = find(listOfBondedNtds==structure{i}(j+numBonds)-1);
            if bondedToNextIndex - length(listOfBondedNtds)/2 <= 0
                bondedToNextIndex = bondedToNextIndex + length(listOfBondedNtds)/2;
            else
                bondedToNextIndex = bondedToNextIndex - length(listOfBondedNtds)/2;
            end
            bondedToNext = listOfBondedNtds(bondedToNextIndex);
            interveningNtds = structure{i}(j)+1:bondedToNext-1;
            if isempty(interveningNtds)
                interveningNtds = structure{i}(j)+1:-1:bondedToNext-1;
            end
            if ~isempty(interveningNtds) && all(~any(interveningNtds == listOfBondedNtds')) %if all the intervening ntds are unbound
                isBulgeLoopOfLengthGreater1 = true;
                if length(interveningNtds) == 1
                    isBulgeLoopOfLength1 = true;
                    isBulgeLoopOfLengthGreater1 = false;
                end
            end
        end
        
        if ~isBulgeLoopOfLengthGreater1
            k = structure{i}(j);
            l = structure{i}(j+numBonds);
            firstNtd = sequenceInNumbers(k); %what is the 5' ntd?
            firstBP = sequenceInNumbers(l); %what is it bonded to?
            threePrimeBondedIndex = find(listOfBondedNtds==l+(2*isParallel-1));
            if threePrimeBondedIndex - length(listOfBondedNtds)/2 <= 0
                threePrime = listOfBondedNtds(threePrimeBondedIndex + length(listOfBondedNtds)/2);
            else
                threePrime = listOfBondedNtds(threePrimeBondedIndex - length(listOfBondedNtds)/2);
            end
            secondNtd = sequenceInNumbers(threePrime); %what is the 3' ntd? %this must exist by the rules we've set up which ntds can form bps.
            secondBP = sequenceInNumbers(l+(2*isParallel-1));


            if firstNtd <=4 && firstBP <= 4 && secondNtd <=4 && secondBP <=4 %don't consider a bond between RNA-RNA on one stem and RNA-DNA on another, (or other combinations)
                ntdType = 0; %RNARNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNARNA,bondEntropyMatrixRNARNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            elseif firstNtd <=4 && firstBP > 4 && secondNtd <=4 && secondBP > 4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            elseif firstNtd > 4 && firstBP <=4 && ~isParallel && secondNtd > 4 && secondBP <=4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,secondBP,secondNtd,firstBP,firstNtd,ntdType);
            elseif firstNtd > 4 && firstBP <=4 && isParallel && secondNtd > 4 && secondBP <=4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstBP,firstNtd,secondBP,secondNtd,ntdType);
            elseif firstNtd > 4 && firstBP > 4 && secondNtd > 4 && secondBP > 4
                ntdType = 2; %DNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixDNADNA,bondEntropyMatrixDNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            end

            FCSs{i}{1+length(FCSs{i})} = {[k,l,l+(2*isParallel-1),threePrime],bondEnergy,bondEntropy};
            FCSsFE(i) = min(FCSsFE(i),FCSs{i}{length(FCSs{i})}{2}-T*FCSs{i}{length(FCSs{i})}{3});
        end
    end
    
    
    
    
    
    %now consider the other end of the stem for a potential terminal mismatch
    if any(realNtds == structure{i}(1+numBonds)-(2*isParallel-1)) &&... %ismember(structure{i}(1+numBonds)-(2*isParallel-1),realNtds) &&...%structure{i}(1+numBonds)-(2*isParallel-1) <= numBP && 
            ~any(firstNtds == structure{i}(1))...%~ismember(structure{i}(1),firstNtds)... %structure{i}(1)~=1  ...
            ...%terminal mismatch 2 doesn't exist if we're dealing with the either end of the sequence.
            && ~any(listOfBondedNtds == structure{i}(1)-1) && ~any(listOfBondedNtds == structure{i}(1+numBonds)-(2*isParallel-1))
        k = structure{i}(1);
        l = structure{i}(1+numBonds);
        
        %going from this direction, we flip 5' to 3' to consider the
        %terminal mismatch
        secondBP = sequenceInNumbers(k-1); %what is the 5' ntd?
        secondNtd = sequenceInNumbers(l-(2*isParallel-1)); %what is it bonded to?
        firstBP = sequenceInNumbers(k); %what is the 3' ntd?
        firstNtd = sequenceInNumbers(l); %what is it bonded to?
        
        %if secondNtd and secondBP could bind but aren't bound in this
        %structure, Lu, Turner, Mathews (NAR, 2006) say that we should
        %treat them as an AC pair (where A replaces the purine and C the
        %pyrimidine)
        realSecondNtd = secondNtd; realSecondBP = secondBP;
        if (secondNtd == 1 && secondBP == 4) || (secondNtd == 4 && secondBP == 1) ||...
                (secondNtd == 2 && secondBP == 3) || (secondNtd == 3 && secondBP == 2) ||...
                (secondNtd == 3 && secondBP == 4) || (secondNtd == 4 && secondBP == 3)
            if secondNtd == 1 || secondNtd == 3 %if secondNtd is a purine
                secondNtd = 1; %the purine is replaced with A 
                secondBP = 2; %and the pyrimidine with C
            else
                secondNtd = 2;
                secondBP = 1;
            end
        end
        
        if firstNtd <=4 && firstBP <= 4 && secondNtd <=4 && secondBP <=4
            ntdType = 0; %RNARNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNARNA,bondEntropyMatrixRNARNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
        elseif firstNtd <=4 && firstBP > 4 && secondNtd <=4 && secondBP > 4
            ntdType = 1; %RNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
        elseif firstNtd > 4 && firstBP <=4 && ~isParallel && secondNtd > 4 && secondBP <=4
            ntdType = 1; %RNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,secondBP,secondNtd,firstBP,firstNtd,ntdType);
        elseif firstNtd > 4 && firstBP <=4 && isParallel && secondNtd > 4 && secondBP <=4
            ntdType = 1; %RNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstBP,firstNtd,secondBP,secondNtd,ntdType);
        elseif firstNtd > 4 && firstBP > 4 && secondNtd > 4 && secondBP > 4
            ntdType = 2; %DNADNA
            [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixDNADNA,bondEntropyMatrixDNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
        end
        
        TMs{i+numStems}{1+length(TMs{i+numStems})} = {[k,l,k-1,l-(2*isParallel-1)],bondEnergy,bondEntropy};
        TMsFE(i+numStems) = min(TMsFE(i+numStems),TMs{i+numStems}{length(TMs{i+numStems})}{2}-T*TMs{i+numStems}{length(TMs{i+numStems})}{3});

        secondNtd = realSecondNtd; secondBP = realSecondBP;
        %we also might have a mismatch-mediated coaxial stack (MMCS)
        if any(listOfBondedNtds == structure{i}(1)-2) || any(listOfBondedNtds == structure{i}(1+numBonds)-2*(2*isParallel-1)) %if either of the next ntds over are bound, 
             if any(listOfBondedNtds == structure{i}(1)-2) && any(listOfBondedNtds == structure{i}(1+numBonds)-2*(2*isParallel-1)) ...
                     && abs(find(listOfBondedNtds==structure{i}(1)-2)-find(listOfBondedNtds == structure{i}(1+numBonds)-2*(2*isParallel-1))) == length(listOfBondedNtds)/2 
                        %but if they're bound to each other, that's a 1x1 internal loop
                %don't do anything if we have a 1x1 internal loop
             else
                 MMCScouldBindWC = false; MMCScouldBindGU = false;
                 if (secondNtd == 1 && secondBP == 4) || (secondNtd == 4 && secondBP == 1) ||...
                        (secondNtd == 2 && secondBP == 3) || (secondNtd == 3 && secondBP == 2)
                    MMCScouldBindWC = true;
                 elseif (secondNtd == 3 && secondBP == 4) || (secondNtd == 4 && secondBP == 3)
                     MMCScouldBindGU = true;
                 end
                 if any(listOfBondedNtds==structure{i}(1)-2)
                     kMinusTwoBPIndex = max(0,find(listOfBondedNtds==structure{i}(1)-2)-length(listOfBondedNtds)/2); %could be either + or - length(listOfBondedNtds)/2
                     if kMinusTwoBPIndex == 0
                         kMinusTwoBP = listOfBondedNtds(find(listOfBondedNtds==structure{i}(1)-2)+length(listOfBondedNtds)/2);
                     else
                         kMinusTwoBP = listOfBondedNtds(kMinusTwoBPIndex);
                     end
                     MMCSs{i+numStems}{1+length(MMCSs{i+numStems})}  = {[k,l,k-1,l-(2*isParallel-1),k-2,kMinusTwoBP],...
                         bondEnergy+MMCSEnergy + extraMMCSEnergyWC*MMCScouldBindWC + extraMMCSEnergyGU*MMCScouldBindGU,bondEntropy+MMCSEntropy};
                     MMCSsFE(i+numStems) = min(MMCSsFE(i+numStems),MMCSs{i+numStems}{length(MMCSs{i+numStems})}{2}-T*MMCSs{i+numStems}{length(MMCSs{i+numStems})}{3});
                 end
                 if any(listOfBondedNtds == structure{i}(1+numBonds)-2*(2*isParallel-1))
                     lMinusTwoBPIndex = max(0,find(listOfBondedNtds==structure{i}(1+numBonds)-2*(2*isParallel-1))-length(listOfBondedNtds)/2); %could be either + or - length(listOfBondedNtds)/2
                     if lMinusTwoBPIndex == 0
                         lMinusTwoBP = listOfBondedNtds(find(listOfBondedNtds==structure{i}(1+numBonds)-2*(2*isParallel-1))+length(listOfBondedNtds)/2);
                     else
                         lMinusTwoBP = listOfBondedNtds(lMinusTwoBPIndex);
                     end
                     MMCSs{i+numStems}{1+length(MMCSs{i+numStems})} = {[k,l,k-1,l-(2*isParallel-1),l-2*(2*isParallel-1),lMinusTwoBP],...
                         bondEnergy+MMCSEnergy + extraMMCSEnergyWC*MMCScouldBindWC + extraMMCSEnergyGU*MMCScouldBindGU,bondEntropy+MMCSEntropy};
                     MMCSsFE(i+numStems) = min(MMCSsFE(i+numStems),MMCSs{i+numStems}{length(MMCSs{i+numStems})}{2}-T*MMCSs{i+numStems}{length(MMCSs{i+numStems})}{3});
                 end
             end
        end
    end
    
    
    %if only one terminal ntd is unbound, we might have a dangling end, or
    %we might have a flush coaxial stack (we don't consider the ends of the
    %sequence here, only later).
    %We also consider the possibility of dangling ends even if both
    %terminal ntds are unbound since we might have a multiloop. 
    if ~any(listOfBondedNtds == structure{i}(1+numBonds)-(2*isParallel-1))... && any(listOfBondedNtds == structure{i}(1)-1)
            && any(realNtds == structure{i}(1+numBonds)-(2*isParallel-1)) %ismember(structure{i}(1+numBonds)-(2*isParallel-1),realNtds)        
        %this describes a 3' dangling end if it's not parallel since the 3' end of the first
            %strand is unbound and the 5' end of the second strand is bound   
            
        %first, the dangling end. 

        firstBoundNtd = sequenceInNumbers(structure{i}(1+numBonds));
        danglingNtd = sequenceInNumbers(structure{i}(1+numBonds)-(2*isParallel-1));
        if (firstBoundNtd == 4 && sequenceInNumbers(structure{i}(1))==3) ||...
                (firstBoundNtd == 3 && sequenceInNumbers(structure{i}(1))==4) %if the pair is a GU base pair,
            %then treat it as AU for the purposes of calculating dangling
            %ends (from Lu, Turner, Mathews (NAR, 2006))
            if firstBoundNtd == 3
                firstBoundNtd = 1; %i.e. change the G to an A
            end
        end
        if danglingNtd > 4
            danglingNtd = danglingNtd - 4;
        end
        if (firstBoundNtd == 0 || danglingNtd == 0)
            %if one of the sequence elements is unknown, don't do anything
        else
            if ~isParallel
                DEs{i+numStems}{1+length(DEs{i+numStems})} = {[structure{i}(1),structure{i}(1+numBonds),structure{i}(1+numBonds)-(2*isParallel-1)],...
                    dangling3Energy(firstBoundNtd,danglingNtd),dangling3Entropy(firstBoundNtd,danglingNtd)};
                DEsFE(i+numStems) = min(DEsFE(i+numStems),DEs{i+numStems}{length(DEs{i+numStems})}{2}-T*DEs{i+numStems}{length(DEs{i+numStems})}{3});
            else
                DEs{i+numStems}{1+length(DEs{i+numStems})} = {[structure{i}(1),structure{i}(1+numBonds),structure{i}(1+numBonds)-(2*isParallel-1)],...
                    dangling5Energy(firstBoundNtd,danglingNtd),dangling5Entropy(firstBoundNtd,danglingNtd)};
                DEsFE(i+numStems) = min(DEsFE(i+numStems),DEs{i+numStems}{length(DEs{i+numStems})}{2}-T*DEs{i+numStems}{length(DEs{i+numStems})}{3});
            end
        end
    end
    %now consider the flush coaxial stack
    if any(listOfBondedNtds == structure{i}(1)-1)
        isBulgeLoopOfLengthGreater1 = false;
        isBulgeLoopOfLength1 = false;
        if ~isParallel %can't have bulge loop if it's parallel
            bondedToNextIndex = find(listOfBondedNtds==structure{i}(1)-1);
            if bondedToNextIndex - length(listOfBondedNtds)/2 <= 0
                bondedToNextIndex = bondedToNextIndex + length(listOfBondedNtds)/2;
            else
                bondedToNextIndex = bondedToNextIndex - length(listOfBondedNtds)/2;
            end
            bondedToNext = listOfBondedNtds(bondedToNextIndex);
            interveningNtds = structure{i}(1+numBonds)+1:bondedToNext-1;
            if isempty(interveningNtds)
                interveningNtds = structure{i}(1+numBonds)+1:-1:bondedToNext-1;
            end
            if ~isempty(interveningNtds) && all(~any(interveningNtds == listOfBondedNtds')) %if all the intervening ntds are unbound
                isBulgeLoopOfLengthGreater1 = true;
                if length(interveningNtds) == 1
                    isBulgeLoopOfLength1 = true;
                    isBulgeLoopOfLengthGreater1 = false;
                end
            end
        end
        
        if ~isBulgeLoopOfLengthGreater1
            k = structure{i}(1)-1;
            lBondedIndex = find(listOfBondedNtds == k);
            if lBondedIndex-length(listOfBondedNtds)/2 <= 0
                l = listOfBondedNtds(lBondedIndex + length(listOfBondedNtds)/2);
            else
                l = listOfBondedNtds(lBondedIndex - length(listOfBondedNtds)/2);
            end
            firstNtd = sequenceInNumbers(k); %what is the 5' ntd?
            firstBP = sequenceInNumbers(l); %what is it bonded to?
            secondNtd = sequenceInNumbers(structure{i}(1)); %what is the 3' ntd? %this must exist by the rules we've set up which ntds can form bps.
            secondBP = sequenceInNumbers(structure{i}(1+numBonds)); %what is it bonded to? %this must exist by the rules we've set up which ntds can form bps.

            if firstNtd <=4 && firstBP <= 4 && secondNtd <=4 && secondBP <=4
                ntdType = 0; %RNARNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNARNA,bondEntropyMatrixRNARNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            elseif firstNtd <=4 && firstBP > 4 && secondNtd <=4 && secondBP > 4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            elseif firstNtd > 4 && firstBP <=4 && ~isParallel && secondNtd > 4 && secondBP <=4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,secondBP,secondNtd,firstBP,firstNtd,ntdType);
            elseif firstNtd > 4 && firstBP <=4 && isParallel && secondNtd > 4 && secondBP <=4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstBP,firstNtd,secondBP,secondNtd,ntdType);
            elseif firstNtd > 4 && firstBP > 4 && secondNtd > 4 && secondBP > 4
                ntdType = 2; %DNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixDNADNA,bondEntropyMatrixDNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            end

            FCSs{i+numStems}{1+length(FCSs{i+numStems})} = {[k,l,structure{i}(1),structure{i}(1+numBonds)],bondEnergy,bondEntropy};
            FCSsFE(i+numStems) = min(FCSsFE(i+numStems),FCSs{i+numStems}{length(FCSs{i+numStems})}{2}-T*FCSs{i+numStems}{length(FCSs{i+numStems})}{3});
        end
    end
    
    if ~any(listOfBondedNtds == structure{i}(1)-1) ...%&& any(listOfBondedNtds == structure{i}(1+numBonds)-(2*isParallel-1))...
            && any(realNtds == structure{i}(1)-1) %ismember(structure{i}(1+numBonds)-(2*isParallel-1),realNtds)   
        %this describes a 5' dangling end if it's not parallel since the 3'
        %end of the second strand is bound and the 5' end of the first strand is unbound
        firstBoundNtd = sequenceInNumbers(structure{i}(1));
        danglingNtd = sequenceInNumbers(structure{i}(1)-1);
        if (firstBoundNtd == 4 && sequenceInNumbers(structure{i}(1+numBonds))==3) ||...
                (firstBoundNtd == 3 && sequenceInNumbers(structure{i}(1+numBonds))==4) %if the pair is a GU base pair,
            %then treat it as AU for the purposes of calculating dangling
            %ends (from Lu, Turner, Mathews (NAR, 2006))
            if firstBoundNtd == 3
                firstBoundNtd = 1; %i.e. change the G to an A
            end
        end
        if danglingNtd > 4
            danglingNtd = danglingNtd - 4;
        end
        if (firstBoundNtd == 0 || danglingNtd == 0)
            %if one of the sequence elements is unknown, don't do anything
        else
            if ~isParallel
                DEs{i+numStems}{1+length(DEs{i+numStems})} = {[structure{i}(1),structure{i}(1)-1,structure{i}(1+numBonds)],...
                    dangling5Energy(firstBoundNtd,danglingNtd),dangling5Entropy(firstBoundNtd,danglingNtd)};
                DEsFE(i+numStems) = min(DEsFE(i+numStems),DEs{i+numStems}{length(DEs{i+numStems})}{2}-T*DEs{i+numStems}{length(DEs{i+numStems})}{3});
            else
                DEs{i+numStems}{1+length(DEs{i+numStems})} = {[structure{i}(1),structure{i}(1)-1,structure{i}(1+numBonds)],...
                    dangling3Energy(firstBoundNtd,danglingNtd),dangling3Entropy(firstBoundNtd,danglingNtd)};
                DEsFE(i+numStems) = min(DEsFE(i+numStems),DEs{i+numStems}{length(DEs{i+numStems})}{2}-T*DEs{i+numStems}{length(DEs{i+numStems})}{3});
            end
        end
    end
    %now consider the flush coaxial stack
    if any(listOfBondedNtds == structure{i}(1+numBonds)-(2*isParallel-1))
        isBulgeLoopOfLengthGreater1 = false;
        isBulgeLoopOfLength1 = false;
        if ~isParallel %can't have bulge loop if it's parallel
            bondedToNextIndex = find(listOfBondedNtds==structure{i}(1+numBonds)+1);
            if bondedToNextIndex - length(listOfBondedNtds)/2 <= 0
                bondedToNextIndex = bondedToNextIndex + length(listOfBondedNtds)/2;
            else
                bondedToNextIndex = bondedToNextIndex - length(listOfBondedNtds)/2;
            end
            bondedToNext = listOfBondedNtds(bondedToNextIndex);
            interveningNtds = structure{i}(1)-1:bondedToNext+1;
            if isempty(interveningNtds)
                interveningNtds = structure{i}(1)-1:-1:bondedToNext+1;
            end
            if ~isempty(interveningNtds) && all(~any(interveningNtds == listOfBondedNtds')) %if all the intervening ntds are unbound
                isBulgeLoopOfLengthGreater1 = true;
                if length(interveningNtds) == 1
                    isBulgeLoopOfLength1 = true;
                    isBulgeLoopOfLengthGreater1 = false;
                end
            end
        end
        
        if ~isBulgeLoopOfLengthGreater1
            kBondedIndex = find(listOfBondedNtds==structure{i}(1+numBonds)-(2*isParallel-1));
            if kBondedIndex - length(listOfBondedNtds)/2 <= 0
                k = listOfBondedNtds(kBondedIndex + length(listOfBondedNtds)/2);
            else
                k = listOfBondedNtds(kBondedIndex - length(listOfBondedNtds)/2);
            end
            l = structure{i}(1+numBonds)-(2*isParallel-1);
            firstNtd = sequenceInNumbers(k); %what is the 5' ntd?
            firstBP = sequenceInNumbers(l); %what is it bonded to?
            secondNtd = sequenceInNumbers(structure{i}(1)); %what is the 3' ntd? %this must exist by the rules we've set up which ntds can form bps.
            secondBP = sequenceInNumbers(structure{i}(1+numBonds));


            if firstNtd <=4 && firstBP <= 4 && secondNtd <=4 && secondBP <=4
                ntdType = 0; %RNARNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNARNA,bondEntropyMatrixRNARNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            elseif firstNtd <=4 && firstBP > 4 && secondNtd <=4 && secondBP > 4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            elseif firstNtd > 4 && firstBP <=4 && ~isParallel && secondNtd > 4 && secondBP <=4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,secondBP,secondNtd,firstBP,firstNtd,ntdType);
            elseif firstNtd > 4 && firstBP <=4 && isParallel && secondNtd > 4 && secondBP <=4
                ntdType = 1; %RNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstBP,firstNtd,secondBP,secondNtd,ntdType);
            elseif firstNtd > 4 && firstBP > 4 && secondNtd > 4 && secondBP > 4
                ntdType = 2; %DNADNA
                [bondEnergy,bondEntropy] =  deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrixDNADNA,bondEntropyMatrixDNADNA,firstNtd,firstBP,secondNtd,secondBP,ntdType);
            end

            FCSs{i+numStems}{1+length(FCSs{i+numStems})} = {[k,l,structure{i}(1),structure{i}(1+numBonds)],bondEnergy,bondEntropy};
            FCSsFE(i+numStems) = min(FCSsFE(i+numStems),FCSs{i+numStems}{length(FCSs{i+numStems})}{2}-T*FCSs{i+numStems}{length(FCSs{i+numStems})}{3});
        end
    end

    
    
    
    
    
    
    
    %if we're dealing with either end of the sequence, we might have a
    %dangling end
    if ~isParallel && any(firstNtds == structure{i}(1)) ... %ismember(structure{i}(1),firstNtds)...%structure{i}(1) == 1 
            && ~any(lastNtds == structure{i}(1+numBonds)) %~ismember(structure{i}(1+numBonds),lastNtds) %structure{i}(1+numBonds) ~= numBP 
            %if the first ntd is bound to any but the last, there is a 3' dangling end
        firstBP = sequenceInNumbers(structure{i}(1+numBonds));
        if (firstBP == 4 && sequenceInNumbers(structure{i}(1))==3) ||...
            (firstBP == 3 && sequenceInNumbers(structure{i}(1))==4) %if the pair is a GU base pair,
            %then treat it as AU for the purposes of calculating dangling
            %ends (from Lu, Turner, Mathews (NAR, 2006))
            if firstBP == 3
                firstBP = 1; %i.e. change the G to an A
            end
        end
        danglingNtd = sequenceInNumbers(1+structure{i}(1+numBonds));
        if danglingNtd > 4
            danglingNtd = danglingNtd - 4;
        end
        if (firstBP == 0 || danglingNtd == 0)
            %if one of the sequence elements is unknown, don't do anything
           %possibleNodeTerms{i} = {possibleNodeTerms{i},{[structure{i}(1+numBonds),1+structure{i}(1+numBonds),structure{i}(1)],0,0}};
        else
            DEs{i+numStems}{1+length(DEs{i+numStems})} = {[1+structure{i}(1+numBonds),structure{i}(1+numBonds),structure{i}(1)],...
                dangling3Energy(firstBP,danglingNtd),dangling3Entropy(firstBP,danglingNtd)};
            DEsFE(i+numStems) = min(DEsFE(i+numStems),DEs{i+numStems}{length(DEs{i+numStems})}{2}-T*DEs{i+numStems}{length(DEs{i+numStems})}{3});
        end
    elseif isParallel && any(firstNtds == structure{i}(1)) &&... %ismember(structure{i}(1),firstNtds) &&...%structure{i}(1) == 1 && 
            ~any(lastNtds == structure{i}(1+numBonds)) %~ismember(structure{i}(1+numBonds),lastNtds) %structure{i}(1+numBonds) ~= numBP 
            %if the first ntd is bound to any but the last, in a pseudoknot, there is a 5' dangling end
        firstBP = sequenceInNumbers(structure{i}(1+numBonds));
        if (firstBP == 4 && sequenceInNumbers(structure{i}(1))==3) ||...
            (firstBP == 3 && sequenceInNumbers(structure{i}(1))==4) %if the pair is a GU base pair,
            %then treat it as AU for the purposes of calculating dangling
            %ends (from Lu, Turner, Mathews (NAR, 2006))
            if firstBP == 3
                firstBP = 1; %i.e. change the G to an A
            end
        end
        danglingNtd = sequenceInNumbers(-1+structure{i}(1+numBonds));
        if danglingNtd > 4
            danglingNtd = danglingNtd - 4;
        end
        if (firstBP == 0 || danglingNtd == 0)
           %possibleNodeTerms{i} = {possibleNodeTerms{i},{[structure{i}(1+numBonds),-1+structure{i}(1+numBonds),structure{i}(1)],0,0}};
        else
           DEs{i+numStems}{1+length(DEs{i+numStems})} = {[structure{i}(1+numBonds),-1+structure{i}(1+numBonds),structure{i}(1)],...
               dangling5Energy(firstBP,danglingNtd),dangling5Entropy(firstBP,danglingNtd)};
           DEsFE(i+numStems) = min(DEsFE(i+numStems),DEs{i+numStems}{length(DEs{i+numStems})}{2}-T*DEs{i+numStems}{length(DEs{i+numStems})}{3});
        end
    elseif ~isParallel && any(lastNtds == structure{i}(1+numBonds))&&... ismember(structure{i}(1+numBonds),lastNtds) &&... %structure{i}(1+numBonds) == numBP && 
            ~any(firstNtds == structure{i}(1))%~ismember(structure{i}(1), firstNtds) %structure{i}(1) ~=1 
            %if the last ntd is bound to any but the first, there is a 5' dangling end
        firstBoundNtd = sequenceInNumbers(structure{i}(1));
        if (firstBoundNtd == 4 && sequenceInNumbers(structure{i}(1+numBonds))==3) ||...
            (firstBoundNtd == 3 && sequenceInNumbers(structure{i}(1+numBonds))==4) %if the pair is a GU base pair,
            %then treat it as AU for the purposes of calculating dangling
            %ends (from Lu, Turner, Mathews (NAR, 2006))
            if firstBoundNtd == 3
                firstBoundNtd = 1; %i.e. change the G to an A
            end
        end
        danglingNtd = sequenceInNumbers(structure{i}(1)-1);
        if danglingNtd > 4
            danglingNtd = danglingNtd - 4;
        end
        if (firstBoundNtd == 0 || danglingNtd == 0)
            %if one of the sequence elements is unknown, don't do anything
        else
            DEs{i+numStems}{1+length(DEs{i+numStems})} = {[structure{i}(1+numBonds),structure{i}(1)-1,structure{i}(1)],...
                dangling5Energy(firstBoundNtd,danglingNtd),dangling5Entropy(firstBoundNtd,danglingNtd)};
            DEsFE(i+numStems) = min(DEsFE(i+numStems),DEs{i+numStems}{length(DEs{i+numStems})}{2}-T*DEs{i+numStems}{length(DEs{i+numStems})}{3});
        end
    elseif isParallel && any(lastNtds == structure{i}(2*numBonds)) %ismember(structure{i}(2*numBonds),lastNtds) %structure{i}(2*numBonds) == numBP  
            %if the last ntd is bound to any but the first, in a pseudoknot, there is a 3' dangling end
        firstBoundNtd = sequenceInNumbers(structure{i}(numBonds));
        if (firstBoundNtd == 4 && sequenceInNumbers(structure{i}(1))==3) ||...
            (firstBoundNtd == 3 && sequenceInNumbers(structure{i}(1))==4) %if the pair is a GU base pair,
            %then treat it as AU for the purposes of calculating dangling
            %ends (from Lu, Turner, Mathews (NAR, 2006))
            if firstBoundNtd == 3
                firstBoundNtd = 1; %i.e. change the G to an A
            end
        end
        danglingNtd = sequenceInNumbers(structure{i}(numBonds)+1);
        if danglingNtd > 4
            danglingNtd = danglingNtd - 4;
        end
        if (firstBoundNtd == 0 || danglingNtd == 0)
            %if one of the sequence elements is unknown, don't do anything
        else
            DEs{i}{1+length(DEs{i})} = {[structure{i}(numBonds),structure{i}(numBonds)+1,structure{i}(numBonds*2)],...
                dangling3Energy(firstBoundNtd,danglingNtd),dangling3Entropy(firstBoundNtd,danglingNtd)};
            DEsFE(i) = min(DEsFE(i),DEs{i}{length(DEs{i})}{2}-T*DEs{i}{length(DEs{i})}{3});
        end
    end
    
end
%%
%Heuristic: first place all FCSs, then all TMs, then all DEs. 
% MMCSs weren't working well (their inclusion seemed to give the
%wrong results) so I moved this below TMs. Effectively, this means that
%MMCSs are never considered (since each MMCS also has a TM).
%For each, pick the best (i.e. first pick the minFE FCS, then the
%minFE TM, then the minFE DE).

ntdsUsed = [];
nodesUsed = [];
for node = 1:numNodes
    [~,i] = min(FCSsFE); %If the minimum value occurs more than once, then min returns the index corresponding to the first occurrence.
    if ~isempty(FCSs{i}) 
        canContinue = false;
        if isempty(ntdsUsed) || ~any(any(FCSs{i}{1}{1} == ntdsUsed'))
            mfeIndex = 1;
            mfe = FCSs{i}{1}{2}-T*FCSs{i}{1}{3};
            if length(FCSs{i})>1
                if FCSs{i}{2}{2}-T*FCSs{i}{2}{3} < mfe
                    mfeIndex = 2;
                end
            end
            canContinue = true;
        elseif length(FCSs{i})>1
            mfeIndex = 2;
            if ~any(any(FCSs{i}{2}{1} == ntdsUsed'))
                canContinue = true;
            end
        end
        if canContinue
            bondEnergyLocal = bondEnergyLocal + FCSs{i}{mfeIndex}{2};
            bondEntropyLocal = bondEntropyLocal + FCSs{i}{mfeIndex}{3};
            ntdsUsed = [ntdsUsed,FCSs{i}{mfeIndex}{1}];
            nodesUsed = [nodesUsed,i];
            DEs{i} = []; %so we don't consider it again
            TMs{i} = [];
            MMCSs{i} = [];
        end
    end
    FCSsFE(i) = 100; %so we don't consider the same FCS again
end


for node = 1:numNodes
    [~,i] = min(TMsFE); %If the minimum value occurs more than once, then min returns the index corresponding to the first occurrence.
    if ~any(nodesUsed == i) && ~isempty(TMs{i})
        canContinue = false; 
        %can only have one possible terminal mismatch.
        if isempty(ntdsUsed) || ~any(any(TMs{i}{1}{1} == ntdsUsed'))
            mfeIndex = 1;
            canContinue = true;
        end
        if canContinue
            bondEnergyLocal = bondEnergyLocal + TMs{i}{mfeIndex}{2};
            bondEntropyLocal = bondEntropyLocal + TMs{i}{mfeIndex}{3};
            ntdsUsed = [ntdsUsed,TMs{i}{mfeIndex}{1}];
            nodesUsed = [nodesUsed,i];
            DEs{i} = []; %so we don't consider it again
        end
    end
    TMsFE(i) = 100; %so we don't consider the same TM again
end


for node = 1:numNodes
    [~,i] = min(MMCSsFE); %If the minimum value occurs more than once, then min returns the index corresponding to the first occurrence.
    if ~any(nodesUsed == i) && ~isempty(MMCSs{i})
        canContinue = false;
        if isempty(ntdsUsed) || ~any(any(MMCSs{i}{1}{1} == ntdsUsed'))
            mfeIndex = 1;
            mfe = MMCSs{i}{1}{2}-T*MMCSs{i}{1}{3};
            if length(MMCSs{i})>1
                if MMCSs{i}{2}{2}-T*MMCSs{i}{2}{3} < mfe
                    mfeIndex = 2;
                end
            end
            canContinue = true;
        elseif length(MMCSs{i})>1
            mfeIndex = 2;
            if ~any(any(MMCSs{i}{2}{1} == ntdsUsed'))
                canContinue = true;
            end
        end
        if canContinue
            bondEnergyLocal = bondEnergyLocal + MMCSs{i}{mfeIndex}{2};
            bondEntropyLocal = bondEntropyLocal + MMCSs{i}{mfeIndex}{3};
            ntdsUsed = [ntdsUsed,MMCSs{i}{mfeIndex}{1}];
            nodesUsed = [nodesUsed,i];
            DEs{i} = []; %so we don't consider it again
            TMs{i} = [];
        end
    end
    MMCSsFE(i) = 100; %so we don't consider the same MMCS again
end



for node = 1:numNodes
    [~,i] = min(DEsFE); %If the minimum value occurs more than once, then min returns the index corresponding to the first occurrence.
    if ~any(nodesUsed == i) && ~isempty(DEs{i})
        canContinue = false; 
        %can only have one possible terminal mismatch.
        if isempty(ntdsUsed) || ~any(any(DEs{i}{1}{1} == ntdsUsed'))
            mfeIndex = 1;
            canContinue = true;
        end
        if canContinue
            bondEnergyLocal = bondEnergyLocal + DEs{i}{mfeIndex}{2};
            bondEntropyLocal = bondEntropyLocal + DEs{i}{mfeIndex}{3};
            ntdsUsed = [ntdsUsed,DEs{i}{mfeIndex}{1}];
            nodesUsed = [nodesUsed,i];
        end
    end
    DEsFE(i) = 100; %so we don't consider the same DE again
end


end





