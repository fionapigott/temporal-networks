% Fiona Pigott
% JAn 29 2014
% MATLAB v. 2012b

% Find and store the indices where the person is 'infected'
% INPUT: 'Adoption' from 'Adoption.mat'
% OUTPUT: 
%     1) whenInfected{i,1} are the indices of every time step where i is 
%           'infected,' whenInfected{i,2} are the indices of every time  
%           step where i is not infected. whenInfected{i,3} is the first 
%           infection in i, OR, if i is never infected, the empty set. 
%       fluctuations: vector with length (# nodes), stores the number
%           of times that each internet adopter transitions
%     2) 'transitions' row 1 shows number of first time adopters per time
%           interval in the first row (positive scalars). Row 2 shows the 
%           total number of people who adopted in each month. Row 3 shows 
%           the total number of people who abandoned in each month.
%     3) 'whoInfected' is a 1D cell array of length (# time steps) with 
%           lists of infected nodes at a given time step
%     4) intervalsR: the probability distribution P(interval = x|the person
%           gave up internet. vector with length (max interval without
%           internet after adopting). intervalsR(x) x = interval length
%           in time steps (months)


% tic

% 1)----------------------------------------------------------------------                               
% Find at which time step each person adopted internet, or (if applicable)
% at which they 'recovered' from internet service, and/or readopted.

[numnodesA, nummatA] = size(Adoption);

whenInfected = cell(numnodesA,3);
A = zeros(1,nummatA);
y = [];
z = [];
fluctuations = zeros(numnodesA,1);
% % debugging
% check1 = 0;
% check0 = 0;
for ii = 1:numnodesA
    A = find(diff(Adoption(ii,:))~=0);
    x = sort([1,A,A+1,nummatA]);
    for jj = 1:4:length(x)-1
        y = [y,x(jj):x(jj+1)];
        if jj+3 <= length(x)
            z = [z,x(jj+2):x(jj+3)];
        end
    end
    if Adoption(ii,1) == 1
        whenInfected{ii,1} = y;
        whenInfected{ii,2} = z;
    else
        whenInfected{ii,2} = y;
        whenInfected{ii,1} = z;
    end
    if isempty(A) %if Adoption(ii,:) is either always 1 or always 0
        if Adoption(ii,1) == 0
            whenInfected{ii,3} = [];
        else
            whenInfected{ii,3} = 1;
        end
    else
        whenInfected{ii,3} = whenInfected{ii,1}(1);
    end
    % % debugging
    % check1 = sum(Adoption(ii,infected{ii,1})) + check1;
    % check0 = sum(Adoption(ii,infected{ii,2})) + check0;
    y = [];
    z = [];
    
    % Store the number of times that each internet adopter 
    % transitioned
    fluctuations(ii) = sum(diff(Adoption(ii,:))~=0);
end

% 2)-------------------------------------------------------------------
% We are also interested to know how many new (first time) adoptions
% occur per month.
transitions = zeros(3,nummatA);
for ii = 1:numnodesA
    if ~isempty(whenInfected{ii,1})
        transitions(1,whenInfected{ii,1}(1)) =...
            transitions(1,whenInfected{ii,1}(1))+1;      
    end
end
% The people who started out as adopters don't count
transitions(1,1) = 0;

% How many adoptions total occur per month (including people who are 
% adopting after abandoning) and how many people abandon ('recover')
% per month
changes = horzcat(zeros(numnodesA,1),diff(Adoption,1,2));
transitions(2,:) = sum(changes>0);
transitions(3,:) = -sum(changes<0);


%3)----------------------------------------------------------------------
whoInfected = cell(nummatA,1);
for m = 1:nummatA
    % fill the cell matrix with the node # for each infected node
    whoInfected{m,1} = find(Adoption(:,m));
end

clear x y z ii jj A
save('AdoptionStats.mat','transitions','whenInfected',...
    'fluctuations','whoInfected')
% toc
%------------------------------------------------------------------------

%4)---------------------------------------------------------------------
% Find the amount of time people spend without internet after
% they go from 1->0

intervalsR = zeros(numnodesA,max(fluctuations));
whichMonths = zeros(numnodesA,max(fluctuations));

for ii = 1:numnodesA
    if ~isempty(whenInfected{ii,1})
        if whenInfected{ii,1}(end)~=nummatA
            infected = diff(horzcat(whenInfected{ii,1},nummatA))-1;
        else
            infected = diff(whenInfected{ii,1})-1;
        end
        whichMonths(ii,1:length(find(infected))) = ...
            whenInfected{ii,1}(find(infected));
        %intervalsR(ii,1:length(find(infected))) = infected(find(infected));
        intervalsR(ii,1) = mean(infected(find(infected)));
    end
end

% make them into prob dists
intervalsR(isnan(intervalsR)) = 0;
intervalsR = intervalsR(find(intervalsR));
intervalsR = accumarray(round(intervalsR),1);
intervalsR = intervalsR/sum(intervalsR);

whichMonths = whichMonths(find(whichMonths));
whichMonths = accumarray(whichMonths,1);
whichMonths = whichMonths/sum(whichMonths);

clear numnodesA nummatA
        





