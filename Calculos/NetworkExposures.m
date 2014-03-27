% Fiona Pigott
% Jan 27 2014
% MATLAB v.2012b

% INPUT: 'Adoption.mat', 'Graphs.mat',
% 'edgeDistribution.mat','TempCorrCoeff.mat'
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
% 3 is COMMENTED OUT, because I don't need it at the moment
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
% 6) Calculate the average difference between the PNE of an adopter,
%    up to the moment of adoption, and the PNE of the population in general
%    up until the same moment
%        OUTPUT: 




% Calculations are performed as defined in Adoption_Networks.pdf (draft).
% The calculation makes use of aggregated, weighted graphs stored in 'data'
% over the time window.

% All output is saved in the file 'NetworkExposures.mat'

tic

%% 1,2)-------------------------------------------------------------------
% Calculate PNE and GNE
% Acessing a point on the graphs: data(i,j,m) = data(row, col, time step)
                 
PNE = zeros(numnodes, nummat);
PNE_NeighborOverlap = PNE;
EC = zeros(numnodes, nummat);
GNE = zeros(1,nummat);
PNEunweighted = PNE;
lastm = 1;

for m = 2:nummat % for each time step
    for ii = 1:numnodes

        % dot product time step of data with the previous 
        % time step of adoption data
        product=dot(unweighted(:,ii,m),Adoption(:,m-1));

        % divide the product for each node by the weight for each node in
        % the time step m
        % then store the column vector of length (# nodes) in PNE
        % Note: product should always be less than degree, since we are 
        % dividing dot products of a vector W*Y (where Y has entries 0/1)
        % and the dot product of W*1 (a vector of 1s).
        PNE(ii,m) = product / k(ii,m);
        
        % Include temporal clustering information
        % PNE(ii,m) = C(ii,m-1)*product / s(ii,m);
        
        % unweighted PNE
        %product=dot(unweighted(:,ii,m),Adoption(:,m));
        %PNEunweighted(ii,m) = product / k(ii,m);
        
        % A.Y/(|A||Y|)
        %PNE(ii,m) = ...
        %    product/sqrt(sum(unweighted(:,ii,m-1))*sum(Adoption(:,m-1)));
        
        
%         % Prof Mauricio's suggestion for PNE calc
%         % commented out, because it takes forever to run
%         if m ~= lastm
%             [EC, ~, ~] = edge_nei_overlap_bu(unweighted(:,:,m));
%             EC(find(isnan(EC))) = 0;
%             lastm = m;
%         end
%         PNE_NeighborOverlap(ii,m) = dot(EC(:,ii),Adoption(:,m-1));
    end
end
PNE(isnan(PNE))=0;
%PNEunweighted(isnan(PNEunweighted))=0;

PNEA = zeros(nummat,1);
PNENonA = zeros(nummat,1);
PNET = PNEA;

for m = 1:nummat-1
    PNEA(m) = sum(PNE(whoInfected{m},m))/sum(Adoption(:,m));
    PNENonA(m) = sum(PNE(setdiff(1:numnodes,whoInfected{m}),m))/...
        (numnodes - sum(Adoption(:,m)));
    PNET(m) = sum(PNE(setdiff(whoInfected{m+1},whoInfected{m}),m))/transitions(2,m);
end

% one def of GNE
GNE = sum(PNE,1);
%GNEunweighted = sum(PNEunweighted,1);

%another possible def of GNE
% GNE = zeros(1,nummat);
% GNEunweihgted = GNE;
% for m = 1:nummat
%     GNE(m) = sum(sum(PNE(:,1:m),2));
%     GNEunweighted = sum(sum(PNEunweighted(:,1:m),2))
% end
    
AdoptGNEline = polyfit(sum(Adoption,1),GNE,1);
%AdoptGNElineunweighted = polyfit(GNEunweighted,sum(Adoption,1),1);

clear m product ii
%------------------------------------------------------------------------

%% 3)---------------------------------------------------------------------
NSI = 0;
% % Calcualte Neighbor Social Influence (NSI)
% % Want the relationship between if after interacting with h, i adopted
% % Relationship between infected -> suceptible 
% 
% NSI = cell(nummat,1);
% for m = 1:nummat
%     % Initialize NSI
%     NSI{m} = zeros(length(whoInfected{m,1}),numnodes);
% end
% 
% % weight of interaction between h (infected) -> i (suceptible)
% for m = 2:nummat % for each time step
%     for jj = 1:length(whoInfected{m-1}) % for each previously infected h
%         h = whoInfected{m-1}(jj); % vector of infected nodes in the
%                                   % previous time step
%         for ii = setdiff(1:numnodes,whoInfected{m-1})
%             % Present interaction between i (currently un-infected node)
%             % and h (previously infected node)
%             ihInteraction = data(ii,h,m);
%             % Present interaction between i and all previously 
%             % infected nodes
%             iInfectedInteraction = sum(data(ii,whoInfected{m-1},m));
%             if iInfectedInteraction ~= 0
%                 NSI{m}(jj,ii) = ihInteraction/iInfectedInteraction;
%             end
%         end
%     end
% end    
% clear m whenInfectedEnd iInfectedInteraction ihInteraction

%% 4)------------------------------------------------------------------
% find the probabilty mass func for the PNE
% pdfPNE = ...
%     accumarray(round(reshape(PNE,prod(size(PNE)),1)*1000)+1,1)...
%     /(numnodes*nummat);
% 
% pmfPNE = zeros(1001,1);
% for n = 1:1001
%     pmfPNE(n) = sum(pdfPNE(1:n));
% end

% 5)------------------------------------------------------------------
%% P(k_i = k| a_i = 1)

tallykAll = ones(nummat*10,1)*-1; % use -1 to init the vector, to 
tallyk = ones(nummat*10,1)*-1;    % easily distinguish between the 
n1 = 1;                           % init vector and times when the 
n2 = 1;                           % value should actually be zero

firstAdopt = zeros(max(transitions(1,:))+1,nummat);
for ii = 1:numnodes
    % list of all of the people who first transitioned 
    % in a given time step (transitioned in m-1 -> m
    if ~isempty(whenInfected{ii,3})
        m = whenInfected{ii,3};
        firstAdopt(find(firstAdopt(:,m)==0,1,'first'),m) = ii;
    end        
end

hAll = cell(nummat-1,1);
hFirst = cell(nummat-1,1);

hAll{1} = firstAdopt(find(firstAdopt(:,1)),1);
hFirst{1} = hAll{1};

% taking into account the k at the first time step during which the node
% is infected: 
for m = 1:nummat
    if m > 1
        % people who become infected in the time step m-1 -> m
        % NOT just new adopters
        hAll{m} = setdiff(whoInfected{m},whoInfected{m-1});
        % for NEW adopters in the time step m -> m+1
        hFirst{m} = firstAdopt(firstAdopt(:,m)~=0,m);
    end
    % make a list of all of the k values for people who transitioned
    tallykAll(n1:n1+length(hAll{m})-1) = k(hAll{m},m);
    n1 = n1 + length(hAll{m});
    % make a list of all of the k values for people who transitioned
    tallyk(n2:n2+length(hFirst{m})-1) = k(hFirst{m},m);
    n2 = n2 + length(hFirst{m});
end

% Bin the values and divide by total adopters to get the probabilities
% of each separate value of k
tallykAll = tallykAll(tallykAll>=0);
tallyk = tallyk(tallyk>=0);
PkGivenAAll = accumarray(tallykAll+1,1)/sum(transitions(2,:));
PkGivenA = accumarray(tallyk+1,1)/sum(transitions(1,:));

% And, in case we decide later that we care, find P(A|k=k)
PAgivenk = accumarray(round(mean(k,2))+1,sum(Adoption,2)>0)./...
    accumarray(round(mean(k,2))+1,1);


%% 6) 
diffPNEAll = zeros(nummat-1,1);
diffPNEFirst = zeros(nummat-1,1);
diffPNERand = zeros(nummat-1,1);
PNEAvgTotalExposure = zeros(nummat-1,1);

diffPNEAllU = zeros(nummat-1,1);
diffPNEFirstU = zeros(nummat-1,1);
diffPNERandU = zeros(nummat-1,1);
PNEAvgTotalExposureU = zeros(nummat-1,1);
for m = 1:nummat
    % everyone's PNE up until m
    PNEAvgTotalExposure(m) = mean(sum(PNE(:,1:m),2));
    % for all users
    diffPNEAll(m) = (mean(sum(PNE(hAll{m},1:m),2)));% - ...
        %PNEAvgTotalExposure(m));%/PNEAvgTotalExposure(m);
    % for only firs-time adopters
    diffPNEFirst(m) = (mean(sum(PNE(hFirst{m},1:m),2)));% - ...
        %PNEAvgTotalExposure(m));%/PNEAvgTotalExposure(m);
    % for a random sample of adopters the same size as hAll
    diffPNERand(m) = ...
        (mean(sum(PNE(randperm(numnodes,length(hAll{m})),1:m),2)));% - ...
        %PNEAvgTotalExposure(m));%/PNEAvgTotalExposure(m);
    
    % the unweighted version:
    % everyone's PNE up until m
    PNEAvgTotalExposureU(m) = mean(sum(PNEunweighted(:,1:m),2));
    % for all users
    diffPNEAllU(m) = (mean(sum(PNEunweighted(hAll{m},1:m),2)));
    % for only firs-time adopters
    diffPNEFirstU(m) = (mean(sum(PNEunweighted(hFirst{m},1:m),2)));
    % for a random sample of adopters the same size as hAll
    diffPNERandU(m) = ...
        (mean(sum(PNEunweighted(randperm(numnodes,length(hAll{m})),1:m),2)));
end


toc


%% Find the giant adopted component

giantAdoptedComp = zeros(1,nummat);

for m = 1:nummat
    A = repmat(Adoption(:,m)',[numnodes, 1]);
    [comps, sizes] = get_components((A & A') & unweighted(:,:,m));
    giantAdoptedComp(m) = max(sizes)/sum(Adoption(:,m));
end

%% Save Network Exposure data
clear hAll hFirst %tallyk tallkAll
if randomizeTime == 1
    save('NetworkExposureRandTime.mat','PNE','GNE','AdoptGNEline','tallyk','giantAdoptedComp')%,'NSI','pmfPNE')
elseif randomizeInfection == 1
    save('NetworkExposureRandAdopt.mat','PNE','GNE','AdoptGNEline','tallyk','giantAdoptedComp')%,'NSI','pmfPNE')
else
    save('NetworkExposure.mat','PNE','GNE','AdoptGNEline','tallyk','giantAdoptedComp')%,'NSI','pmfPNE')
end

