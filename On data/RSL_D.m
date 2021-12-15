function [H,tests,SC,flag] = RSL_D(D, Mb, alpha, alpha_Mb, varargin)
% input:
% D: Data matrix with size (number of samples)*n
% Mb: n*n Markov boundary matrix. Entries can be 0 or 1
% alpha: significance level for CI tests. for example 0.01
% alpha_Mb: significance level for updating Mb. for example 2/n^2
% 
% output:
% H: The matrix of learned skeleton. Entries can be 0 or 1
% tests: the number of performed CI tests
% SC: Size of Conditioning sets. It is a (n+1)*1 vector such that i-th 
% entry indicates the number of conditioning sets with size i-1


%**************************************************************************
%********************** Start of the Algorithm ****************************
%**************************************************************************
n = size(D,2);
if nargin>4
    V = varargin{1};
    tests = varargin{2};
    SC=varargin{3};
    flag = varargin{4};
else
    V = ones(n,1);
    tests = 0;
    SC = zeros(n+1,1);
    flag = true(1,n);
end
if sum(V)==1
    H = zeros(n);
    return
else
    [X,tests,SC,flag] = FindRemovable(V,Mb,tests,SC,flag,D,alpha);
    [N,tests,SC] = FindNeighborhood(X,Mb(X,:),tests,SC,D,alpha);
    [Mb,tests,SC,flag] = UpdateMb(X,N,Mb,tests,SC,flag,D,alpha_Mb);
    V(X) = 0;
    [H,tests,SC,flag] = RSL_D(D, Mb,alpha, alpha_Mb, V, tests, SC, flag);
    H(X,N==1)=1;
    H(N==1,X)=1;
    return
end
end

%**************************************************************************
%************************* Main Functions *********************************
%**************************************************************************

function [X,tests,SC,flag] = FindRemovable(V,Mb,tests,SC,flag,D,alpha)
Mbs = sum(Mb);
[~,ind1] = sort(Mbs);
ind2 = ind1(flag(ind1)>0);
for XX=ind2
    Mb_XX = find(Mb(XX,:));
    [isR,tests,SC] = isRemovable(XX,Mb_XX,tests,SC,D,alpha);
    flag(XX) = false;
    if isR
        X = XX;
        return
    end
end
ind2 = ind1(V(ind1)>0);
X = ind2(1);
end


function [isR,tests,SC] = isRemovable(X,Mb_X,tests,SC,D,alpha)

isR = true;
mbs = length(Mb_X);
for j =1:(mbs-1)
    Y = Mb_X(j);
    for k = (j+1):mbs
        Z = Mb_X(k);
        S = [X, Mb_X([1:(j-1), (j+1):(k-1), (k+1:mbs)])];
        CI = CI_Test(Y,Z,S,D,alpha);
        ls = length(S);
        SC(ls+1)= SC(ls+1)+1;
        tests = tests + 1;
        if CI
            isR = false;
            return
        end
    end
end
end


function [N,tests,SC] = FindNeighborhood(X,Mb_X,tests,SC,D,alpha)
N = Mb_X;
Mbx = find(Mb_X);
Mbs = length(Mbx);
for i=1:(Mbs)
    Y = Mbx(i);
    for  j= [1:(i-1),(i+1):Mbs]
        Z=Mbx(j);
        S=mysetdiff(Mbx,[Y;Z]);
        ls = length(S);
        SC(ls+1)= SC(ls+1)+1;
        tests = tests+1;
        if CI_Test(X,Y,S,D,alpha)
            N(Y)=0;
            break
        end
    end
end

end


function [Mb,tests,SC,flag]= UpdateMb(X,N,Mb,tests,SC,flag,D,alpha)

Mb_X = find(Mb(X,:));
for Y = Mb_X
    Mb(X,Y) = 0;
    Mb(Y,X) = 0;
    flag(Y) = true;
end
mbs = length(Mb_X);
ns = sum(N);
if mbs > ns
    return
end
for i=1:(mbs-1)
    Y = Mb_X(i);
    for j=(i+1):mbs
        Z = Mb_X(j);
        if nnz(Mb(Y,:))>nnz(Mb(Z,:))
            S = mysetdiff(find(Mb(Z,:)),Y);
        else
            S = mysetdiff(find(Mb(Y,:)),Z);
        end
        CI = CI_Test(Y,Z,S,D,alpha);
        ls = length(S);
        SC(ls+1)= SC(ls+1)+1;
        tests = tests+1;
        if CI
            Mb(Z,Y) = 0;
            Mb(Y,Z) = 0;
            flag(Y) = true;
            flag(Z) = true;
        end
    end
end
end

%**************************************************************************
%************************* Minor Functions ********************************
%**************************************************************************

function C = mysetdiff(A,B)
% MYSETDIFF Set difference of two sets of positive integers
% (much faster than built-in setdiff)
% C = mysetdiff(A,B)
% C = A \ B = { things in A that are not in B }
% Original by Kevin Murphy, modified by Leon Peshkin

if isempty(A)
    C = [];
    return;
elseif isempty(B)
    C = A;
    return;
else % both non-empty
    bits = zeros(1, max(max(A), max(B)));
    bits(A) = 1;
    bits(B) = 0;
    C = A(logical(bits(A)));
end
end
