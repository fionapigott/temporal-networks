% Fiona Pigott
% Jan 27 2014
% MATLAB v.2012b

% INPUT: 'Adoption.mat' and 'Graphs.mat'
% OUTPUT:
% 1) Calculate the Personal Network Exposure (PNEi) of a susceptible 
%    (non yet adopter) vertex i
%       OUTPUT: 'PNE' a matrix with dimension (# nodes) X (# time steps)  
%       where each entry represents the PNE of a node i at a time m.
%       Note that the first column will be identically '0'
% 2) Calculate Global Network Exposure (GNE) exposure of the entire
%    network to a social contagion in a specific month m.
%       OUTPUT: 'GNE' a vector of length (# time steps) with one 
%               value of GNE per step. Note GNE(1)=0.
%               'AdoptGNEline' line fit of Adoption as a function of GNE
%               AdoptGNEline(1) = y-intercept, (2) = slope
% 3) Calculate the Neighbor Social Influence (NSI) will estimate how an 
%    already infected neighbour h will contribute to the PNEi of each 
%    uninfected vertex i to which it is connected.
%       OUTPUT: 'NSI' a 1D cell array of dimension (# time steps). 
%       The array stores a 2D array per cell with dimension (#infected) X
%       (# nodes), this stores the interaction between a given infected  
%       h and each uninfected node i.       
%           ex: NSI{m} is a matrix of interaction between those listed
%               infected nodes h and every other node i
%               And index NSI{m}(h,i) would give the influence of the 
%               infected h (# of the node h = whoInfected{m,1}(x)) on i
%               in the case where i is infected or where i and h are not 
%               connected, NSI{m}(h,i) = 0
% 4) Calculate the probability mass func for PNE
%        OUTPUT: pmfPNE, vector with length (1001)
%        pmfPNE(x) = P(PNE < x)
% 5) Calculate the probability distribution P(K=k|Adopted) to try to 
%    characterize the difference between a group of adopters and the group 
%    as a whole
%        OUTPUT: 'PkGivenA' and 'PkGivenAAll', vectors with length (max #
%        transitions), where tally__(k) = P(k_i = k|a_i = 1). The
%        difference between PkGivenA and PkGivenAAll is that PkGivenAAll
%        accounts for all nodes, and PkGivenA only considers participating 
%        nodes. As some non-participating nodes are listed as internet 
%        adopters, it is necessary to condsider both.




% Calculations are performed as defined in Adoption_Networks.pdf (draft).
% The calculation makes use of aggregated, weighted graphs stored in 'data'
% over the time window.

% All output is saved in the file 'NetworkExposures.mat'

% tic

%% 1,2)-------------------------------------------------------------------
% Calculate PNE and GNE
% Acessing a point on the graphs: data(i,j,m) = data(row, col, time step)

weights = squeeze(sum(data,1)); % find the total weight of phone calls
                                % made by each person in each time step
                                % dim (# nodes) X (# time steps)                 
PNE = zeros(numnodes, nummat);
GNE = zeros(1,nummat);
%product = zeros(numnodes,1);

for ii = 1:numnodes
%     if isempty(whenInfected{ii,3})
%         whenInfectedEnd = nummat;
%     else
%         whenInfectedEnd = whenInfected{ii,3};
%     end
    for m = 2:nummat % for each time step
        product=dot(data(:,ii,m),Adoption(:,m-1));% dot product time step 
                                               % of data with the previous 
                                               % time step of adoption data                                  
        % divide the product for each node by the weight for each node in
        % the time step m
        % then store the column vector of length (# nodes) in PNE
        % Note: product should always be less than weights, since we are 
        % dividing dot products of a vector W*Y (where Y has entries 0/1)
        % and the dot product of W*1 (a vector of 1s).
        PNE(ii,m) = product / weights(ii,m);
    end
end
PNE(isnan(PNE))=0;
GNE = sum(PNE,1);

AdoptGNEline = polyfit(GNE,sum(Adoption,1),1);

clear m product ii
%------------------------------------------------------------------------

%% 3)---------------------------------------------------------------------

% Calcualte Neighbor Social Influence (NSI)
% Want the relationship between if after interacting with h, i adopted
% Relationship between infected -> suceptible 

NSI = cell(nummat,1);
for m = 1:nummat
    % Initialize NSI
    NSI{m} = zeros(length(whoInfected{m,1}),numnodes);
end

% weight of interaction between h (infected) -> i (suceptible)
for m = 2:nummat % for each time step
    for jj = 1:length(whoInfected{m-1}) % for each previously infected h
        h = whoInfected{m-1}(jj); % vector of infected nodes in the
                                  % previous time step
        for ii = setdiff(1:numnodes,whoInfected{m-1})
            % Present interaction between i (currently un-infected node)
            % and h (previously infected node)
            ihInteraction = data(ii,h,m);
            % Present interaction between i and all previously 
            % infected nodes
            iInfectedInteraction = sum(data(ii,whoInfected{m-1},m));
            if iInfectedInteraction ~= 0
                NSI{m}(jj,ii) = ihInteraction/iInfectedInteraction;
            end
        end
    end
end    
clear m whenInfectedEnd

%% 4)------------------------------------------------------------------
% find the probabilty mass func for the PNE
pdfPNE = ...
    accumarray(round(reshape(PNE,prod(size(PNE)),1)*1000)+1,1)...
    /(numnodes*nummat);

pmfPNE = zeros(1001,1);
for n = 1:1001
    pmfPNE(n) = sum(pdfPNE(1:n));
end

% 5)------------------------------------------------------------------
%% P(k_i = k| a_i = 1)

tallykAll = ones(nummat*10,1)*-1; % use -1 to init the vector, to 
tallyk = ones(nummat*10,1)*-1;    % easily distinguish between the 
n1 = 1;                           % init vector and times when the 
n2 = 1;                           % value should actually be zero

firstAdopt = zeros(max(transitions(1,:))+1,nummat);
for ii = 1:numnodes
    % list of all of the people who first transitioned 
    % in a given time step
    if ~isempty(whenInfected{ii,3})
        m = whenInfected{ii,3};
        firstAdopt(find(firstAdopt(:,m)==0,1,'first'),m) = ii;
    end        
end

for m = 1:nummat-1
    % people who become infected in the time step m -> m+1
    % NOT just new adopters
    hAll = setdiff(whoInfected{m+1},whoInfected{m});
    % make a list of all of the k values for people who transitioned
    tallykAll(n1:n1+length(hAll)-1) = k(hAll,m);
    n1 = n1 + length(hAll);
    % for NEW adopters in the time step m -> m+1
    hFirst = firstAdopt(firstAdopt(:,m)~=0,m);
    % make a list of all of the k values for people who transitioned
    tallyk(n2:n2+length(hFirst)-1) = k(hFirst,m);
    n2 = n2 + length(hFirst);
end

% Bin the values and divide by total adopters to get the probabilities
% of each separate value of k
tallykAll = tallykAll(tallykAll>=0);
tallyk = tallyk(tallyk>=0);
PkGivenAAll = accumarray(tallykAll+1,1)/sum(transitions(2,:));
PkGivenA = accumarray(tallyk+1,1)/sum(transitions(1,:));


% Save Network Exposure data
if randomizeTime == 1
    save('NetworkExposureRandTime.mat','PNE','GNE','AdoptGNEline','NSI','pmfPNE')
elseif randomizeInfection == 1
    save('NetworkExposureRandAdopt.mat','PNE','GNE','AdoptGNEline','NSI','pmfPNE')
else
    save('NetworkExposure.mat','PNE','GNE','AdoptGNEline','NSI','pmfPNE')
end

