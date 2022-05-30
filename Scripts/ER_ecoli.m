clear
close all
t = 0;
%set number of monte carlo trials
% N > 1e4 not recommended due to run time
N = 1e4;
% initialize FFLrand vector
FFLrand = zeros(1,N);
% number of nodes
n = 423;
% number of edges
E = 578;
% number of FFLs
realFFL = 42;
tic
for z = 1:N
    % initialize matrix of all nodes with no edges
    M = zeros(n,n);
    for k = 1:E
        % pick random edge location
        RE = randi(n,1,2);
        % set matix indicies
        i = RE(1);
        j = RE(2);
        % check to make sure edge has not been placed here already
        if M(i,j) == 1
            %randomize until new edge location selected
            while M(i,j) == 1
                RE = randi(n,1,2);
                i = RE(1);
                j = RE(2);
            end
        end
        M(i,j) = 1;
    end
    FFL = 0;
    % count FFLs in this random simulation
    FFLrand(z) = FFLcount(M);
end
toc

% statistics of MC simulation
averageFFL = mean(FFLrand)
STDofFFL = std(FFLrand)
zscore = (realFFL-averageFFL)/STDofFFL
