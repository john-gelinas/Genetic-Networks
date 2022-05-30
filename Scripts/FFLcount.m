%Feed Forward Loop Counter
function FFL = FFLcount(M)
n = length(M);

FFL = 0;
% count FFLs in this simulation
for i = 1:n
    for j = 1:n
        % test if regulation exists
        if M(i,j) == 1
            % test for all other nodes i regulates
            for q = 1:n
                % test if FFL exists
                if M(i,q) == 1 && M(j,q) == 1 && M(j,i) == 0 && M(q,j) == 0 && M(q,i) == 0
                    %count FFL
                    FFL = FFL + 1;
                end
            end
        end
    end
end


end