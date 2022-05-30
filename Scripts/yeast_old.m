clear
close all
load('yeast.mat')
% i = TF, j = operon, therefore i regulates j
i2 = i;
i = j;
j = i2;

n = max([i;j]);
edges = length(i);
M = zeros(n);
Mpos = zeros(n);
Mneg = zeros(n);

for z = 1:edges
    % make adjacency matrix
    M(i(z),j(z)) = 1;
end

autoreg = trace(M);

FFL = FFLcount(M);

% find in-degrees and out-degrees
for count = 1:n
    outdeg(count) = sum(M(:,count));
    indeg(count) = sum(M(count,:));
end
deg = outdeg + indeg;
% find mutual degree
mutdeg = zeros(1,n);
for i = 1:n
    for j = 1:n
        if M(i,j) == M(j,i) && M(j,i) == 1
            if i ~= j
                mutdeg(i) = mutdeg(i) + 1;
            end
        end
    end
end
Kmean = mean(outdeg);
KKminus1mean = mean(outdeg.*(outdeg-1));
Rmean = mean(indeg);
RRminus1mean = mean(indeg.*(indeg-1));
Mmean = mean(mutdeg);
RKmean = mean(outdeg.*indeg);
FFLpredicted = KKminus1mean*RKmean*RRminus1mean/(Kmean^3);

figure1 = figure('Color',[1 1 1]);
histogram(indeg,'Normalization','probability')
xlabel('In-Degree')
ylabel('Frequency')
axes1 = get(figure1,'CurrentAxes');
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontWeight','bold');

figure2 = figure('Color',[1 1 1]);
histogram(outdeg,'Normalization','probability')
xlabel('Out-Degree')
ylabel('Frequency')
axes2 = get(figure2,'CurrentAxes');
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontWeight','bold');

% make degree vectors
for degree = 0:max(outdeg)
    outdegdist(degree+1) = sum(outdeg == degree);
    kout(degree+1) = degree;
end

for degree = 0:max(indeg)
    indegdist(degree+1) = sum(indeg == degree);
    kin(degree+1) = degree;
end

for degree = 0:max(deg)
    degdist(degree+1) = sum(deg == degree);
    k(degree+1) = degree;
end
% fit degree vectors to power law to find gamma
count = 0;
for z = 2:length(k)
    if degdist(z) ~= 0
        count = count + 1;
        degdistfit(count) = degdist(z)/n;
        kfit(count) = k(z);
    end
end
[xData, yData] = prepareCurveData( kfit, degdistfit );
% Set up fittype and options.
ft = fittype( 'x^-gamma', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.959492426392903;
[fitresult, gof] = fit( xData, yData, ft, opts );
gamma = coeffvalues(fitresult);
confidenceintervals = confint(fitresult);

