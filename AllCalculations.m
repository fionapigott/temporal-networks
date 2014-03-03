% Fiona Pigott
% Jan 29 2014
% MATLAB v.2012b

% Perform all calculations (thus far) on given temporal network and
% adoption data set

tic
% LOAD DATA (if necessary)-----------------------------------------------
cd('Calculos')
% AgregateMat % Script to read in "data." If Graphs.mat is  already
              % in the folder, don't run this, it takes a long time

% % load data--Don't do this every time. If data, unweighted, numnodes and
% % nummat are already in the workspace, don't run this
load('Graphs.mat') % name of adjacency matricies -> 'data'
%                         % number of nodes in the network -> 'numnodes'
%                         % number of temporal snapshots -> 'nummat'

% Make fake data, if you want. Setting all options to 0
% will not generate fake data
% Choose an option by setting it to "1"
randomizeTime = 0;
randomizeInfection = 0;
% Script below loads internet adoption data from a .mat file.
% Save data as a matrix 'Adoption' with dim (# nodes) X (# time steps)
% '1' indicates that the user at that node has adopted internet
% '0' indicates that the user has not yet adopted internet.
cd('..')             
cd('Atributos')
load('Adoption.mat')
cd('..')
cd('Calculos')

% generates fake data, if randomize___ is set to 1
generateFakeData

% Calculate data about adoption (when first adoption occurs, etc)
AdoptionStats
% % Or just load the data
% load('AdoptionStats.mat')

% Calculate Network Exposure for not-yet adopters
NetworkExposures
% % Or just load the data
% load('NetworkExposure.mat')

% Calculate temporal correlation coefficient
TempCorrCoeff
% TempCorrCoeffAllCombos
% % Or just load the data
% load('TempCorrCoeff.mat')

% Calculate the probabilty distributions of edges
edgeDistribution
% % Or just load the data
% load('edgeDistribution.mat')

% Graph Everything
cd('..')
cd('Graficos')
MakeFigures
cd('..')
% ------------------------------------------------------------------------
toc
