% Fiona Pigott
% March 10 2014
% MATLAB v2012b

% Look for a relationship betweent the frequency of a node's internet
% use and their likeihood of spreading internet adoption to a neighbor

% Note: I don't actually see a correlation

% OUTPUT: bar graph of rate of adoption of neighbors for each category of
%         user

%% User Level vs contagion

% initialize a matrix to hold all of the properties of interest
userLevel = [(1:numnodes)',sum(Adoption,2);,zeros(numnodes,1),...
    zeros(numnodes,1),zeros(numnodes,1)];
userLevel = userLevel(find(userLevel(:,2)),:);
userRange = [0,0];
userRange(1) = length(userLevel(:,1));
userRange(2) = max(userLevel(:,2));
[A B] = sort(userLevel(:,2));
userLevel = userLevel(B,:);

% number of connections that each user has, average
for ii = 1:userRange(1)
    node = userLevel(ii,1);
    userLevel(ii,3) = mean(k(node,:));
end

% frst time adoptions matrix
FirstAdoptions = zeros(numnodes,nummat);
for ii = 1:numnodes
    for m = 1:nummat
        FirstAdoptions(ii,whenInfected{ii,3}) = 1;
    end
end

%number of transitions in the neighborhood of each user
for ii = 1:userRange(1)
    node = userLevel(ii,1);
    userLevel(ii,4) = nummat + 1 - whenInfected{node,3};
    for m = whenInfected{node,3}:nummat
        userLevel(ii,5) = userLevel(ii,5) + ...
            dot(unweighted(node,:,m),FirstAdoptions(:,m));
    end
end

% graph the results
bargraph = zeros(11,1);
for n = 1:11
    bargraph(n) = sum(userLevel((n-1)*26+1:n*26,5))./...
        sum(userLevel((n-1)*26+1:n*26,4));
end

figure
bar(bargraph)
xlabel({'User strength, 1-11. 1 being the users who have the fewest months';...
    'of internet use, 11 having the most total months of use'})
ylabel('Av rate (/month) of adoption of neighbors of each category of user')

clear userLever bargraph node A B userRange

%% Commented out: think about user level as a function of time spent online
% % let's only think about unp until 9/2004,
% % becuase after that things get complicated
% cd('..')
% cd('Atributos')
% load('weightedAdopt')
% cd('..')
% cd('Calculos')
% userLevel = sum(weightedAdopt(:,1:56),2);
% userLevel = [(1:numnodes)',userLevel,zeros(numnodes,1)];
% userLevel = userLevel(find(userLevel(:,2)),:);
% userRange = [0,0];
% userRange(1) = length(userLevel(:,1));
% range(2) = max(userLevel(:,2));
% [A B] = sort(userLevel(:,2));
% userLevel = userLevel(B,:);
% 
% 
% %number of transitions in the neighborhood of each user
% adoptions = diff(Adoption,1,2);
% adoptions = [zeros(numnodes,1),adoptions];
% adoptions(find(adoptions < 0)) = 0;
% for ii = 1:userRange(1)
%     for m = whenInfected{ii,3}:nummat
%         userLevel(ii,3) = userLevel(ii,3) + ...
%             dot(unweighted(userLevel(ii,1),:,m),adoptions(:,m));
%     end
% end

% let's only think about unp until 9/2004,
% becuase after that things get complicated

    
    
    
    