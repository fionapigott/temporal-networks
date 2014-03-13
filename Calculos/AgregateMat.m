% Fiona Pigott
% Jan 23 2014
% MATLAB v. 2012b

% Script to create the matrix of temporal networks and store that
% matrix in "data" (saves data, as well as the number of nodes and times 
% as a .mat file)
% Note that the matrix has 1147x1147x120 elements:
% on my Mac laptop it ran in ~5 min.

% This can be used to aggregate other .csv files into a 3D matrix
% with one .csv per layer as long as:
    % each files is the same dimensions
    % comma delimited
    % navegate to the correct data folder
    % each file name has the same number of characters
% Note that the files will be read in in alphabetical order

% INPUT: Must be a folder with only .csv data files in chronological order
%        in this case, Datos/GephiMat

% OUTPUT:
%   data: all of the weighted adjacency matrices for the series, 
%         dim: (# nodes) X (# nodes) X (# time steps)
%   unweighted: unweighted version of 'data'
%   nummat: the number of time steps (the number of matrices)
%   numnodes: the number of nodes in the graph

% tic
% Move to the data folder
cd('..')
cd('Datos/GephiMat') %Must be a folder with only
                     %.csv data files in chronological
                     % order

% Obtain file names
filescell = dir('*.csv');
filescell = {filescell.name};
nummat = length(filescell);

files = repmat('1234A.bcd',nummat,1); % initialize an array to 
                                      % hold the strings of file names
                                      % if changing the file input, put a 
                                      % string that has the same num of 
                                      % chars as the file name
for n = 1:length(filescell)
    files(n,:) = filescell{n};
end
%clearvars filescell

test = csvread(files(1,:)); % read in test, to find the size of the matrix
sz = size(test);
data = zeros(sz(1)-1,sz(2)-1,nummat);
unweighted = data;
A = zeros(size(test));
checkDiag = eye(size(test)-1)==0;

clearvars test

for n = 1:nummat
    % load the adjacency matrix for one month
    A = csvread(files(n,:));
    if sum(sum(isnan(A)))>0
        Disp('Your data is formatted incorrectly: it contains NaN')
        break
    end
    % order the matrix in increasing order of phone numbers
    A = sortrows(sortrows(A)')';
    % remove the rows with phone numbers
    A = A(2:sz(1),2:sz(2));
    % Ensure that the diagonal entries are all "0"
    A = A .* checkDiag; 
    % Make contacts that are one-directional (A->B but not A<-B) zero
    % We only want to count contacts that occur mutually over one month
    unweightedA = A > 0;
    unweightedA = unweightedA .* unweightedA';
    % store the unweighted matrix in 'unweighted'
    unweighted(:,:,n) = unweightedA;
    % make symmetric (sum A->B & B->A) and 
    % store the weighted matrix in data
    data(:,:,n) = (A .* unweightedA) + (A .* unweightedA)';
end 

sz = sz-1;
numnodes = sz(1);
clearvars A n sz unweightedA checkDiag files

% Move back to "Calculos"
cd('..')
cd('..')
cd('Calculos')

% save 'data' as a .mat file
save('Graphs.mat','data','unweighted','numnodes','nummat')
% toc