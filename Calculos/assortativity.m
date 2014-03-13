% Fiona Pigott
% March 13 2014
% MATALAB v2012b

% Find the assortativity coeeficient in the network over time
% for assortativity of adopters and degree assortativity

% References: [1] Mixing Patterns in Networks (M J Newman)
%             [2] Assortative Mixing in Networks (M J Newman)


% Find assortativity of adopters (i.e., how common is it for adopeters
% to link to adopters.
% r > 0 -> assortative network, to some degree

% OUTPUT: r, a vector with length (# of time steps) that gives the
%         assortativity coefficient for adoption over time.
%         rk, a vector with length (# of time steps) that gives the
%         degree assortativity coefficient over time.

% initialize values
mixingMat = zeros(2,2,nummat);
r = zeros(1,nummat);

% Pick category criteria.
% In this case, I am looking for the assortativity between the group of 
% nodes who, at some point, become adopters, or the assortativity of the
% group with the potential to adopt.
A = repmat((sum(Adoption,2)>0)',[numnodes, 1]); % potential to adopt
%A = repmat((sum(Adoption,2)>12)',[numnodes, 1]); % adoption "level"

% r is a scalar value for each month
for m = 1:nummat

    % % An alternative to the category of potential adopters, we can look
    % % at the assortativity of the group of adopters in a particular month
    % A = repmat(Adoption(:,m)',[numnodes, 1]); % actual adoption
    
    % number of links, without double counting, (the matrix is symmetric)
    All = sum(sum(triu(unweighted(:,:,m)))); 
    % links between adopters
    AtoA = sum(sum(triu(((A & A') & unweighted(:,:,m)))));
    % links that have at least one adopter
    AtoAny = sum(sum(triu(((A | A') & unweighted(:,:,m)))));
    % links that have no adopters
    NtoN = All - AtoAny;
    % links from one adopter to one non-adopter
    NtoA = AtoAny - AtoA;
    
    % define the mixing matrix, as outlined by Newman in ref [1]
    mixingMat(:,:,m) = [AtoA, NtoA/2; NtoA/2, NtoN]./All;
    % calculate the assortativity coefficient as given in ref [1]
    r(m) = (trace(mixingMat(:,:,m)) - sum(sum(mixingMat(:,:,m)^2)))...
        /(1-sum(sum(mixingMat(:,:,m)^2)));
    
    % pick a time step to evaluate error, to get a feel for the 
    % error (too costly to calculate for every time step, and it shouldn't
    % change much).
    % error eq from ref [1]
    if m == ceil(nummat/2)
        error = (1/All)*(sum(sum(mixingMat(:,:,m)^2))+...
            sum(sum(mixingMat(:,:,m)^2))^2 - ...
            sum(sum(mixingMat(:,:,m)^2*mixingMat(:,:,m)')) - ...
            sum(sum(mixingMat(:,:,m)*mixingMat(:,:,m)'^2)))...
            /(1-sum(sum(mixingMat(:,:,m)^2)));
    end   
end

clear A All AtoA NtoN NtoA

% Evaluate degree assortativity of the network
% as defined in ref [2]

% initialize a matrix to hold the scalar properties of links that 
% we are interested in.
properties = cell(1,nummat);
rk = zeros(1,nummat);

% rk is calculated for every time step
for m = 1:nummat
    
    % list the links (start, end) in the graph
    % only use upper triangular part because it's symmetric 
    [from, to] = find(triu(unweighted(:,:,m)));
    fromk = zeros(size(from));
    tok = zeros(size(to));
    % find the degree of the start and end nodes and store
    for ii = 1:length(from)
        fromk(ii) = k(from(ii),m);
        tok(ii) = k(to(ii),m);
    end
    % properties{m} = [node 1 in an edge, node 2 in an edge, ...
    %                  k of node 1, k of node 2]
    properties{m} = [from, to, fromk, tok];
    
    % 1/M, M = the total number of edges involved in the graph at t = m
    Minv = 1/length(from);
    
    % eveluate r as given in ref [2]
    rk(m) = (Minv*sum(properties{m}(:,3).*properties{m}(:,4)) - ...
        (Minv/2*sum(properties{m}(:,3)+properties{m}(:,4)))^2)/...
        (Minv/2*sum(properties{m}(:,3).^2+properties{m}(:,4).^2)+...
        (Minv/2*sum(properties{m}(:,3)+properties{m}(:,4)))^2); 
end

clear properties Minv m from fromk to tok ii
