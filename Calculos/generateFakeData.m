% Fiona Pigott
% Jan 31 2014
% MATLAB 2012b

% Make time-randomized data set--------------------------
if randomizeTime == 1;
    randTime = randperm(nummat);
    data = data(:,:,randTime);
    unweighted = unweighted(:,:,randTime);
    Adoption = Adoption(:,randTime);
    save('randomTimeGraphs.mat','data','unweighted','nummat','numnodes')
    cd('..')
    cd('Atributos')
    save('randomTimeAdoption.mat','Adoption')
    cd('..')
    cd('Calculos')
end
%--------------------------------------------------------

% % Make random-links graph sets-----------------------------
% % a) randomize the links in each month seperately
% %    (preserving the weights/ number of calls per month)
% m = 1;
% A = reshape(triu(data(:,:,m),1),1,numnodes*numnodes);
% scramble = randperm(numnodes*numnodes);
% A = A(scramble);

% %--------------------------------------------------------

% Make random infections data set------------------------
if randomizeInfection == 1
%     numAdoptions = sum(sum(Adoption));
%     randAdopt = randperm(nummat*numnodes,numAdoptions);
%     Adoption = zeros(numnodes,nummat);
%     Adoption(randAdopt) = 1;
%     cd('..')
%     cd('Atributos')
%     save('randomTimeAdoption.mat','Adoption')
%     cd('..')
%     cd('Calculos')

% same number of infections per timestep, but randomly assigned
    numAdoptions = sum(Adoption);
    Adoption = zeros(numnodes, nummat);
    for m = 1:nummat
        randAdopt = randperm(numnodes,numAdoptions(m));
        Adoption(randAdopt,m)=1;
    end
    cd('..')
    cd('Atributos')
    save('randomTimeAdoption.mat','Adoption')
    cd('..')
    cd('Calculos')
end
%--------------------------------------------------------

%  % Make random infection data set with no 'recovery'-----
% %--------------------------------------------------------

% % Make random infection data set with similar recovery
% % rate----------------------------------------------------
% %---------------------------------------------------------

% Make a data set where contagion is definitely present----
% %--------------------------------------------------------