clear
close all
%this program aims to implement a "switching" algorithm by taking a real
%transcriptional network and switching pairs of edges in order to randomize
%a network that maintains the same in degree and out degree as the original
% runtime is approx 40 min on original machine for N = 1000
tic
load('ecoli.mat')
% i = operon, j = TF, therefore j regulates i
% type defines type of regulation (positive = 1, negative = 2, both = 3)

n = max([i;j]); %number of nodes
edges = length(i); %number of edges
M = zeros(n);
Mpos = zeros(n);
Mneg = zeros(n);

for z = 1:edges
    % make adjacency matrix
    M(i(z),j(z)) = 1;
    % make adjacency matrix for positive regulation only
    if type(z) == 1
        Mpos(i(z),j(z)) = 1;
        % make adjacency matrix for negative regulation only
    elseif type(z) == 2
        Mneg(i(z),j(z)) = 1;
    end
end

N = 10;
%generate numerous random networks
for trial = 1:N
    % perform switches numerous times
    % Milo et al -> # of switches >= # of edges (set to 10x to be safe)
    switches = 10*edges;
    Mrand = M;
    
    for swnum = 1:switches
        %find edges of current matrix
        [row,col] = find(Mrand);
        %pick random pair of edges
        redges = randi(edges,2,1);
        %edge 1
        e1r = row(redges(1));
        e1c = col(redges(1));
        %edge 2
        e2r = row(redges(2));
        e2c = col(redges(2));
        %check that old edge pairs arent self-edge
        if e1r == e1c || e2r == e2c
            %skip current switch if either old edge pair is self edge
            continue
        end
        %check that new edge pairs dont already exist
        if Mrand(e1r,e2c) == 1 || Mrand(e2r,e1c) == 1
            %skip current switch if either new edge pair exists
            continue
        end
        %check that new edge pairs arent self-edge
        if e1r == e2c || e2r == e1c
            %skip current switch if either new edge pair is self edge
            continue
        end
        %remove old edges
        Mrand(e1r,e1c) = 0;
        Mrand(e2r,e2c) = 0;
        %add edges between new switched nodes
        Mrand(e1r,e2c) = 1;
        Mrand(e2r,e1c) = 1;
        if sum(sum(Mrand)) ~= edges
            error('Error Occured. Random simulation has different number of edges than original')
        end
    end
    
    autoreg(trial) = trace(Mrand);
    FFL(trial) = FFLcount(Mrand);
    
    % find in-degrees and out-degrees
    for count = 1:n
        outdeg(count) = sum(Mrand(:,count));
        indeg(count) = sum(Mrand(count,:));
    end
    deg = outdeg + indeg;
    % find mutual degree
    mutdeg = zeros(1,n);
    for i = 1:n
        for j = 1:n
            if Mrand(i,j) == Mrand(j,i) && Mrand(j,i) == 1
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
    FFLpredicted(trial) = KKminus1mean*RKmean*RRminus1mean/(Kmean^3);
    
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
    gamma(trial) = coeffvalues(fitresult);
    confidenceintervals = confint(fitresult);
end

gammamean = mean(gamma);
gammastd = std(gamma);
FFLmean = mean(FFL);
FFLstd = std(FFL);
FFLpredmean = mean(FFLpredicted);
FFLpredstd = std(FFLpredicted);

toc

