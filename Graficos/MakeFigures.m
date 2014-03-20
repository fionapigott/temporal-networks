% Fiona Pigott
% Jan 27 2014
% MATLAB v.2012b
% 

close all

%% Create list of dates for the x-labels
dateStr = char('01-01-1998', '02-01-1998', '03-01-1998', '04-01-1998', ...
    '05-01-1998', '06-01-1998', '07-01-1998', '08-01-1998', '09-01-1998',...
    '10-01-1998', '11-01-1998', '12-01-1998', '01-01-1999', '02-01-1999',...
    '03-01-1999', '04-01-1999', '05-01-1999', '06-01-1999', '07-01-1999',...
    '08-01-1999', '09-01-1999', '10-01-1999', '11-01-1999', '12-01-1999',...
    '01-01-2000', '02-01-2000', '03-01-2000', '04-01-2000', '05-01-2000',...
    '06-01-2000', '07-01-2000', '08-01-2000', '09-01-2000', '10-01-2000',...
    '11-01-2000', '12-01-2000', '01-01-2001', '02-01-2001', ...
    '04-01-2001', '05-01-2001', '06-01-2001', '07-01-2001', '08-01-2001',...
    '09-01-2001', '10-01-2001', '11-01-2001', '12-01-2001', '01-01-2002',...
    '02-01-2002', '03-01-2002', '04-01-2002', '05-01-2002', '06-01-2002',...
    '07-01-2002', '08-01-2002', '09-01-2002', '10-01-2002', '11-01-2002',...
    '12-01-2002', '01-01-2003', '02-01-2003', '03-01-2003', '04-01-2003',...
    '05-01-2003', '06-01-2003', '07-01-2003', '08-01-2003', '09-01-2003',...
    '10-01-2003', '11-01-2003', '12-01-2003', '01-01-2004', '02-01-2004',...
    '03-01-2004', '05-01-2004', '06-01-2004', '07-01-2004',...
    '08-01-2004', '09-01-2004', '10-01-2004', '11-01-2004', '12-01-2004',...
    '01-01-2005', '02-01-2005', '03-01-2005', '04-01-2005', '05-01-2005',...
    '06-01-2005', '07-01-2005', '08-01-2005', '09-01-2005', '10-01-2005',... 
    '11-01-2005', '12-01-2005', '01-01-2006', '02-01-2006', '03-01-2006',... 
    '04-01-2006', '05-01-2006', '06-01-2006', '07-01-2006', '08-01-2006',... 
    '09-01-2006', '10-01-2006', '11-01-2006', '12-01-2006', '01-01-2007',...  
    '02-01-2007', '03-01-2007', '04-01-2007', '05-01-2007', '06-01-2007',...
    '07-01-2007', '08-01-2007', '09-01-2007', '11-01-2007',...
    '12-01-2007');
dates = datenum(dateStr);
%% Graph the Global Network Exposure for each month--------------------
%Plot Data
GNEperMonth = figure('name','GNEperMonth');
% subplot(2,2,2)
plot(dates, GNE,'-ob','MarkerSize',3)
%Set Ticks
labels = datestr(dates(1:12:120), 'yyyy');
set(gca, 'XTick', dates(1:12:120));
set(gca, 'XTickLabel', labels);
rotateXLabels( gca, 30 ) % rotateXLabel is a function downloaded from
                         % mathworks forums. refer to liscensing info
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
%Label Axes and Set Title
xlabel('Date')
ylabel('Global Network Exposure')
title('Network Exposure per Month')
if randomizeTime == 1
    title({'Network Exposure per Month';...
        'for a randomly ordered time series'})
end
if randomizeInfection == 1
    title({'Network Exposure per Month';...
        'for randomly generated infection data'})
end
%% Graph adopters per month--------------------------------------------
AdoptersperMonth = figure('name','AdoptersperMonth');
hold on
plot(dates, sum(Adoption,1),'-ob','MarkerSize',3);
plot(dates, transitions(1,:),'-om','MArkerSize',3);
plot(dates, sum(Adoption,1))
plot(dates, transitions(1,:),'m')
label = legend('Total Adopters','New Adopters');
set(label,'Location','NorthWest')
hold off
%Set Ticks
labels = datestr(dates(1:12:120), 'yyyy');
set(gca, 'XTick', dates(1:12:120));
set(gca, 'XTickLabel', labels);
rotateXLabels( gca, 30 ) % rotateXLabel is a function downloaded from
                         % mathworks forums. refer to liscensing info                         
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
%Label Axes and Set Title
xlabel('Date')
ylabel('Internet Adopters')
title('Internet Users per Month')
if randomizeTime == 1
    title({'Internet Users per Month';...
        'for a randomly ordered time series'})
end
if randomizeInfection == 1
    title({'Internet Users per Month';...
        'for randomly generated infection data'})
end
clear label
%% Graph transitions for each month
trans = figure('name','trans');
hold on
plot(dates,transitions(2,:),'-og','MarkerSize',3)
plot(dates,transitions(3,:),'-or','MarkerSize',3)
plot(dates,transitions(2,:)+transitions(3,:),'-ok','MarkerSize',3)
plot(dates,zeros(size(dates)),'k--')
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
%Set Ticks
labels = datestr(dates(1:12:120), 'yyyy');
set(gca, 'XTick', dates(1:12:120));
set(gca, 'XTickLabel', labels);
rotateXLabels( gca, 30 ) % rotateXLabel is a function downloaded from
                         % mathworks forums. refer to liscensing info 
%Label Axes and Set Title
label = legend('Adopt ("infections")','Abandon ("recovery")',...
    'Net Adoptions');
set(label,'Location','NorthWest')
ylabel('Number of people')
xlabel('Date')
title('Net Adoption Over Time')
if randomizeTime == 1
    title({'Net Adoption Over Time';...
        'for a randomly ordered time series'})
end
if randomizeInfection == 1
    title({'Net Adoption Over Time';...
        'for randomly generated infection data'})
end
clear label
%% Graph the GNE for each month vs total adopters per month-------------
GNEvsTotalAdopters = figure('name','GNEvsTotalAdopters');
hold on
plot(sum(Adoption,1),GNE,'.') 
x = 1:max(sum(Adoption,1));
plot(x,x*AdoptGNEline(1)+AdoptGNEline(2))
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
%Label Axes and Set Title
xlabel('Total Internet Adopters')
ylabel('Global Network Exposure')
title('GNE per month vs total adopters per month')
slope = num2str(roundn(AdoptGNEline(1),-2));
int = num2str(roundn(AdoptGNEline(2),-2));
text(x(end)-10,x(end)/2*AdoptGNEline(1)+AdoptGNEline(2),{'Adoption(GNE) = ',[slope,...
     '*GNE + ',int]},... 
     'HorizontalAlignment','right',...
     'FontSize',12);
if randomizeTime == 1
    title({'Total Adopters per month vs GNE per month';...
        'for a randomly ordered time series'})
end
if randomizeInfection == 1
    title({'Total Adopters per month vs GNE per month';...
        'for randomly generated infection data'})
end
%% Graph nearest neighbors for each month vs time-----------------------
NearestNeighbors = figure('name','NearestNeighbors');
hold on
plot(dates, averagek,'-ob','MarkerSize',3);
plot(dates, R, '-or','MarkerSize',3);
plot(dates,averageknn,'-og','MarkerSize',3);
axis([dates(1) dates(end) 0 max(averageknn)+5])
%Set Ticks
labels = datestr(dates(1:12:117), 'yyyy');
set(gca, 'XTick', dates(1:12:117));
set(gca, 'XTickLabel', labels);
rotateXLabels( gca, 30 ) % rotateXLabel is a function downloaded from
                         % mathworks forums. refer to liscensing info                         
set(gca,'FontSize',12);
set(findall(gcf,'type','text'),'FontSize',12);
%Label Axes and Set Title
xlabel('Date')
ylabel('Degree')
title('Degree Estimations vs Time')
label = legend('<k>','<k^2>/<k>','<k_{nn}>');
set(label,'Location','NorthEast');
clear label
if randomizeTime == 1
    title({'Degree Estimations vs Time';...
        'for a randomly ordered time series'})
end
if randomizeInfection == 1
    title({'Degree Estimations vs Time';...
        'for randomly generated infection data'})
end
hold off
%% Graph # of neighbors vs # of neighbors of neighbors-----------------
knnvsk = figure('name','knnvsk');
hold on
plot(1:max(max(k)),knnSorted(:,1),'.','MarkerSize',10);   
set(gca,'FontSize',12);
set(findall(gcf,'type','text'),'FontSize',12);
%Label Axes and Set Title
xlabel('Degree k')
ylabel('Average k_{nn} of nodes with degree k')
title('Degree of nearest neighbors as a function of degree')
%set(label);
clear label
if randomizeTime == 1
    title({'Degree vs degree of nearest neighbors';...
        'for a randomly ordered time series'})
end
if randomizeInfection == 1
    title({'Degree vs degree of nearest neighbors';...
        'for randomly generated infection data'})
end
hold off
%% Currently unused: Weighted knn-----------------
% knnW =figure('name','knnW');
% hold on
% plot(knn(:,[1 60 117]), knnweighted(:,[1 60 117]),'.');   
% set(gca,'FontSize',12);
% set(findall(gcf,'type','text'),'FontSize',12);
% %Label Axes and Set Title
% xlabel('knn')%'Degree k of the node i (k_i)')
% ylabel('knnw')
% %ylabel('Weighted average degree of the nearest neighbors of i (k_{nn})')
% title('Degree vs weighted degree of nearest neighbors')
% label = legend('Jan 1998','Jan 2003','Dec 2007');
% set(label);
% clear label
% hold off
%% P
ProbDistk = figure('name','ProbDistk');
hold on
plot((sum(N(2:end,:),2)/sum(sum(N(2:end,:),2))),'k','LineWidth',2)
plot(POnlyParticipants(:,[1 59 117]))
axis([1,100,0,.08])
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
xlabel('degree (k)')
ylabel('Fraction of nodes with degree k')
label = legend('Total All Months','Jan 1998','Jan 2003','Dec 2007');
set(label,'Location','NorthEast')
title('Probability distribution of the degree of nodes in the network')
if randomizeTime == 1
    title({'Probability distribution of the degree of nodes in the network';...
        'for a randomly ordered time series'})
end
if randomizeInfection == 1
    title({'Probability distribution of the degree of nodes in the network';...
        'for randomly generated infection data'})
end
hold off
%-----------------------------------------------------------------------
%% P(k_i = k| a_i = 1)
PkGivenAdopt = figure('name','PkGivenAdopt');
hold on
plot(0:(size(N,1)-1),(sum(N,2)/(numnodes*nummat)),'k','LineWidth',2)
plot((sum(N(2:end,:),2)/sum(sum(N(2:end,:),2))),'b','LineWidth',2)
plot((1:length(PkGivenAAll))-1,PkGivenAAll,'-og','LineWidth',1.5,'MarkerSize',3)
plot((1:length(PkGivenA))-1,PkGivenA,'-om','LineWidth',1.5,'MarkerSize',3)
axis([-1,max(max(k)),0,.3])
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
%Label Axes and Set Title
ylabel('Probability')
xlabel('k (number of nearest neighbors)')
label = legend('P(k_i = k)','P(k_i = k|k_i \neq 0)','P(k_i = k|Adopted)',...
    'P(k_i = k|Adopted for the first time)');
set(label,'Location','NorthEast')
title('Probabilty distribution of adoption')
hold off

PkGivenAdoptBox = figure('name','PkGivenAdoptBox');
kvect = reshape(k,numnodes*nummat,1);
kno0 = k(find(k~=0));
groups = cell(length(kvect)+length(kno0)+length(tallyk)+length(tallykAll),1);
groups(1:length(kvect)) = {' '};%{'Degree (all nodes)'};
groups((1+length(kvect)):(length(kno0)+length(kvect))) = ...
    {'  '};%{'Degree (participating nodes)'};
groups((length(kno0)+length(kvect)+1):(length(kno0)+length(kvect)+length(tallyk))) =...
    {'   '};%{'Degree of adopters (first time adopters)'};
groups((length(kno0)+length(kvect)+length(tallyk)+1):...
    (length(kno0)+length(kvect)+length(tallyk)+length(tallykAll))) =...
    {'    '};%{'Degree of adopters (each transition)'};
boxplot([kvect;kno0;tallyk;tallykAll],groups)

%----------------------------------------------------------------------
%% Overlap (commented out, no need to re-graph)
% %Doesn't change, and it's re-calculated unless you run "TempCorrCoeffAll"
% %so don't plot it every time
% OverlapPlot = figure('name','OverlapPlot');
% pcolor(dates,dates,overlap);
% shading flat
% axis square
% colorbar%('location','southoutside')
% title('Topological overlap between all time steps')
% %xlabel('Month (m)')
% %ylabel('Month (n)')
% %Set Ticks
% labels = datestr(dates(1:12:120), 'yyyy');
% set(gca, 'XTick', dates(1:12:120));
% set(gca, 'XTickLabel', labels);
% set(gca, 'YTick', dates(1:12:120));
% set(gca, 'YTickLabel', labels);
% rotateXLabels( gca, 30 ) % rotateXLabel is a function downloaded from
%                          % mathworks forums. refer to liscensing info 
% set(gca,'FontSize',12)
% set(findall(gcf,'type','text'),'FontSize',12)
%--------------------------------------------------------------------

%-------------------------------------------------------------------

%----------------------------------------------------------------------
%% Sum vs k
svsk = figure('name','svsk');
plot(1:length(sFunction(:,1)),sFunction(:,1),'.')
set(gca,'FontSize',12);
set(findall(gcf,'type','text'),'FontSize',12);
%Label Axes and Set Title
xlabel('Degree of i (k_i)')
ylabel('Average value of s_i')
title('Strength of the node (s_i) as a function of k')
hold off
%% Clustering coefficient
clusteringcoeff = figure('name','clusteringcoeff');
plot(1:max(max(k)),CFunction,'.')
set(gca,'FontSize',12);
set(findall(gcf,'type','text'),'FontSize',12);
%Label Axes and Set Title
xlabel('Degree of i (k_i)')
ylabel('Average value of c_i')
title('Clustering coefficient of the node (c_i) as a function of k')
%% <PNE> adopters, <PNE> all
PNEdiff = figure('name','PNEdiff');
hold on
%plot(dates(1:end),diffPNEFirstorig,'.k','MarkerSize',10)
%plot(dates,PNEAvgTotalExposureOrig,'k','LineWidth',2)
plot(dates(1:end),diffPNEFirst,'.b','MarkerSize',10)
plot(dates,PNEAvgTotalExposure,'b','LineWidth',2)
%plot(dates(1:end),diffPNEFirst,'.m')
%'Average PNE of first-time adopters at m',...
    %'Average PNE of all nodes at m',...
label = legend('Average PNE of adopters at m',...
    'Average PNE of all nodes at m');
    %'First instance of adoption in the node (first time adopter)');
set(label,'Location','NorthWest')
%Set Ticks
labels = datestr(dates(1:12:120), 'yyyy');
set(gca, 'XTick', dates(1:12:120));
set(gca, 'XTickLabel', labels);
rotateXLabels( gca, 30 ) % rotateXLabel is a function downloaded from
                         % mathworks forums. refer to liscensing info                         
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
%Label Axes and Set Title
xlabel('Date')
ylabel('<PNE>_m')%/<PNE>_m')
title({'PNE of adopters and non-adopters per month'})
hold off
%% Assortativity by degree
degreeAssort = figure('name','degreeAssort');
plot(rk,'.','MarkerSize',10)
labels = datestr(dates(1:12:120), 'yyyy');
set(gca, 'XTick', dates(1:12:120));
set(gca, 'XTickLabel', labels);
rotateXLabels( gca, 30 ) % rotateXLabel is a function downloaded from
                         % mathworks forums. refer to liscensing info  
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
%Label Axes and Set Title
xlabel('Date')
ylabel('Assortativity coefficient r')%/<PNE>_m')
title({'Assortativity of nodes by degree'})
%% Assortativity by adopter/non-adopter
adoptAssort = figure('name','adoptAssort');
hold on
plot(dates, r,'.m','MarkerSize',10)
plot(dates, rAll,'.','MarkerSize',10)
hold off
label = legend('Adopters (in the time step)',...
    'Adopters (adopt at any time)');
set(label,'Location','NorthWest')
labels = datestr(dates(1:12:120), 'yyyy');
set(gca, 'XTick', dates(1:12:120));
set(gca, 'XTickLabel', labels);
rotateXLabels( gca, 30 ) % rotateXLabel is a function downloaded from
                         % mathworks forums. refer to liscensing info  
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
%Label Axes and Set Title
xlabel('Date')
ylabel('Assortativity coefficient r')%/<PNE>_m')
title('Assortativity of adopters and non-adopters')

%% Probability of Adoption given degree k
ProbAgk = figure('name','ProbAgk');
plot(PAgivenk,'.','MarkerSize',10) 
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
xlabel('degree (k)')
ylabel('Probability of adoption')
title('P(Adoption |degree = k)')
%% % Save the graphs
%% if...

if randomizeTime == 1
    cd('RandomTimes')
end
if randomizeInfection == 1
    cd('RandomInfection')
end

print(GNEperMonth,'-depsc','GNEperMonth.eps');
%close(GNEperMonth)
print(AdoptersperMonth,'-depsc','AdoptersperMonth.eps'); 
% close(AdoptersperMonth)
print(GNEvsTotalAdopters,'-depsc','GNEvsTotalAdopters.eps');
% close(GNEvsTotalAdopters)
%savefig(GNEvsTotalAdopters,'GNEvsTotalAdopters.fig');
print(NearestNeighbors,'-depsc','NearestNeighbors.eps'); 
% close(NearestNeighbors)
print(ProbDistk,'-depsc','ProbDistk.eps'); 
% close(ProbDistk)
print(PkGivenAdopt,'-depsc','PkGivenAdopt.eps'); 
print(knnvsk,'-depsc','knnvsk.eps'); 
% close(knnvsk)
print(svsk,'-depsc','svsk.eps'); 
% close(svsk)
print(trans,'-depsc','trans.eps'); 
% close(trans)
% print(OverlapPlot,'-depsc','overlap.eps'); 
print(clusteringcoeff,'-depsc','clusteringcoeff.eps')
% close(clusteringcoeff)
print(PNEdiff,'-depsc','PNEdiff.eps')
% close(PNEdiff)
print(degreeAssort,'-depsc','degreeAssort.eps')
print(adoptAssort,'-depsc','adoptAssort.eps')
print(ProbAgk,'-depsc','ProbAgk.eps')

if randomizeTime == 1
    cd('..')
end
if randomizeInfection == 1
    cd('..')
end

clear GNEperMonth AdoptersperMonth NewAdoptersperMonth GNEvsNewAdopters
clear RperMonth kperMonth k2perMonth startDate endDate Adoptersvsk %dates
clear edgeDistribution networkExposure







