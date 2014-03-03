% Fiona Pigott
% Jan 24 2014
% MATLAB v.2012b

% Calculate temporal correlation coefficient 'TCC' as described in
% "Graph Metrics for Temporal Networks" by Nicosia et al.
% With corrections for disconnected graphs described in TCCcorrection.pdf
% The exact method described is only useful for unweighted graphs.

% INPUT: 'unweighted', 'numnodes' and 'nummat' (all from Graphs.mat)

% OUTPUT:
%   TCC: scalar value, average temporal correlation coeff (C) over all time
%   Ctime: vector with length (# time steps) which shows the temporal 
%          correlation between each month and the next month over time.


% tic
Ctime = zeros(nummat, 1); % Initialize a vector to store average
                        % temporal overlap for each timestep
Ci = 0;          % Initialize variables 
numerator = 0;   % to store sums
denominator = 0;
amSum = 0;
aMSum = 0;
% C2 = zeros(120,1147); %debugging

links = squeeze(sum(unweighted,2));
links(find(links~=0))=1;
numlinks = sum(links);

for m = 1:nummat - 1 % loop to calculate avg temporal overlap in
                     % between two time steps
    for ii = 1:numnodes % step through each node
        for jj = 1:numnodes % evaluate the edge i <-> all other nodes
            am = unweighted(ii,jj,m);
            aM = unweighted(ii,jj,m+1);
            numerator = am*aM + numerator;
            amSum = am + amSum;
            aMSum = aM + aMSum;
        end
        denominator = sqrt(amSum * aMSum);
        if denominator > 0
 %          C2(m,ii) = numerator/denominator; %debugging
            Ci = numerator/denominator + Ci; % sum topological overlap 
                                             % for each node
        end
        numerator = 0;  % reset summing variables
        denominator = 0;
        amSum = 0;
        aMSum = 0;
    end

    Ctime(m) = Ci/(max(numlinks(m),numlinks(m+1)));
                              % store avg topological overlap over nodes i
                              % that are participating in the graph
                              % per timestep
    Ci = 0;        % reset summing variable
end
TCC = sum(Ctime)/(nummat-1);

clear ii jj m amSum aMSum denominator numerator am aM Ci links numlinks 

if randomizeTime == 1
    save('TempCorrCoeffRandTime.mat','Ctime','TCC');
else
    save('TempCorrCoeff.mat','Ctime','TCC');
end
% toc