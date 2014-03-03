% Fiona Pigott
% Jan 24 2014
% MATLAB v.2012b

% Calculate temporal correlation coefficient 'TCC' as described in
% "Graph Metrics for Temporal Networks" by Nicosia et al.
% The exact method described is only useful for unweighted graphs, 
% create an unweighted, undirected graph based on "data", then evaluate.
% OUTPUT:......

% Below (commented) are lines of code to 
%   1) Read in the .csv files and save the output, IF that is not already
%   done. Once the file Graphs.mat is saved in the "Calculos" folder,
%   this does not need to be re-run (unless data changes).
%   2) Load the saved variables into the current workspace, IF they do not
%   already appear. This will need to be run once per MATLAB session.
% Both are time consuming

% AgregateMat % Script to read in "data." If Graphs.mat is  already
              % in the folder, don't run this, it takes a long time

% % load data--Don't do this every time. If data, numnodes and
% % nummat are already in the workspace, don't run this

% load('Graphs.mat') % name of adjacency matricies -> 'data'
                        % number of nodes in the network -> 'numnodes'
                        % number of temporal snapshots -> 'nummat'
                        
% acessing a point on the graphs: data(i,j,m) = data(row, col, time step)
tic
overlap = zeros(nummat, nummat); % Initialize a vector to store average
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

for m = 1:nummat % loop to calculate avg temporal overlap in
                 % between two time steps
    for n = 1:nummat
        for ii = 1:numnodes % step through each node
            for jj = 1:numnodes % evaluate the edge i <-> all other nodes
                am = unweighted(ii,jj,m);
                aM = unweighted(ii,jj,n);
                numerator = am*aM + numerator;
                amSum = am + amSum;
                aMSum = aM + aMSum;
            end
            denominator = sqrt(amSum * aMSum);
            if denominator > 0
                Ci = numerator/denominator + Ci; % sum topological overlap 
                                                 % for each node
            end
            numerator = 0;  % reset summing variables
            denominator = 0;
            amSum = 0;
            aMSum = 0;
        end
        overlap(m,n) = Ci/(max(numlinks(m),numlinks(n)));
        % store avg topological overlap over 
                                   % nodes i per timestep
        Ci = 0;        % reset summing variable
    end
end
%%TCC = sum(Ctime)/(nummat-1);

clear ii jj m amSum aMSum denominator numerator am aM Ci %Ctime
save('Overlap.mat','overlap');
toc