% Fiona Pigott
% February 24 2014
% MATLAB v2012b

%close all

cd('..')

Adoptionsvsk = figure(16);
plot(averagek,sum(Adoption),'.','MarkerSize',10)
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
xlabel('<k> (average number of nearest neighbors per timestep)')
ylabel('Total Adoptions per timestep')
title({'Relationship between total adoptions';...
    'and average number of nearest neighbors'})
print(Adoptionsvsk,'-depsc','Adoptionsvsk.eps');

AdoptionsvsGNE = figure(17);
plot(GNE,sum(Adoption),'.','MarkerSize',10)
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
xlabel('GNE per timestep')
ylabel('Total Adoptions per timestep')
title('Relationship between total adoptions and GNE')
print(AdoptionsvsGNE,'-depsc','AdoptionsvsGNE.eps');

GNEvsk = figure(18);
plot(averagek,GNE,'.','MarkerSize',10)
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
ylabel('GNE per timestep')
xlabel('<k> (average number of nearest neighbors per timestep)')
title('Relationship between GNE and <k>')
print(GNEvsk,'-depsc','GNEvsk.eps');

randomizeInfection = 1;
AllCalculations

cd('Graficos')
cd('RandomInfection')
cd('PreserveMonthlyTotals')

AdoptionsvsGNERand = figure(19);
plot(GNE,sum(Adoption),'.','MarkerSize',10)
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
xlabel('GNE per timestep')
ylabel('Total Adoptions per timestep')
title({'Relationship between total adoptions and GNE';...
    'for randomly generated adoption data'})
print(AdoptionsvsGNERand,'-depsc','AdoptionsvsGNE.eps');

AdoptionsvskRand = figure(20);
plot(averagek,GNE,'.','MarkerSize',10)
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
ylabel('GNE per timestep')
xlabel('<k> (average number of nearest neighbors per timestep)')
title({'Relationship between GNE and <k>';...
    'for randomly generated adoption data'})
print(AdoptionsvskRand,'-depsc','AdoptionsvsGNE.eps');

cd('..')
cd('..')
cd('..')





