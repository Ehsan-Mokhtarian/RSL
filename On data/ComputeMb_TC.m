function [Mb, tests] = ComputeMb_TC(D,alpha_Mb)
% Total Conditioning, see Pellet and Elisseeff

tests = 0;
n = size(D,2);
Mb = zeros(n);

for X=1:(n-1)
    for Y=(X+1):n
        S = [1:(X-1),(X+1):(Y-1),(Y+1):n];
        CI = CI_Test(X,Y,S,D,alpha_Mb);
        tests = tests+1;
        if ~CI
            Mb(X,Y)=1;
            Mb(Y,X)=1;
        end
    end
end
end

