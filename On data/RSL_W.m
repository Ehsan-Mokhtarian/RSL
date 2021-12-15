function [H,tests,SC,flag] = RSL_W(D, V, Mb, alpha, alpha_Mb, m, varargin)

% input:
% D: Data matrix with size (number of samples)*n
% V: Vertices vector 
% Mb: n*n Markov boundary matrix. Entries can be 0 or 1
% alpha: significance level for CI tests. for example 0.01
% alpha_Mb: significance level for updating Mb. for example 2/n^2
% m: an upper bound on the clique number of the graph

% output:
% H: The matrix of learned skeleton. Entries can be 0 or 1
% tests: the number of performed CI tests
% SC: Size of Conditioning sets. It is a (n+1)*1 vector such that i-th
% entry indicates the number of conditioning sets with size i-1

%**************************************************************************
%********************** Start of the Algorithm ****************************
%**************************************************************************

if nargin > 6
    tests = varargin{1};
    SC = varargin{2};
    flag = varargin{3};
else
    n = length(V);
    tests = 0;
    SC = zeros(n+1,1);
    flag = true(1,n);
end
if sum(V)==1
    n = length(V);
    H = zeros(n);
    return
else
    [X,tests,SC,flag] = FindRemovable_G(V,Mb,m,tests,SC,flag,D,alpha);
    [N,tests,SC] = FindNeighborhood(X,Mb(X,:)==1,m,tests,SC,D,alpha);
    [Mb,tests,SC,flag] = UpdateMb(X,N,Mb,tests,SC,flag,D,alpha_Mb);
    V(X) = 0;
    [H,tests,SC,flag] = RSL_W(D,V,Mb,alpha,alpha_Mb,m,tests,SC,flag);
    H(X,N)=1;
    H(N,X)=1;
    return
end
end

% ************************** Major Functions ******************************

function [X,tests,SC,flag] = FindRemovable_G(V,Mb,m,tests,SC,flag,D,alpha)
Mbs = sum(Mb);
[~,ind1] = sort(Mbs);
ind2 = ind1(flag(ind1)>0);
for XX=ind2
    [isR,tests,SC] = isRemovable_G(XX,Mb,m,tests,SC,D,alpha);
    flag(XX) = false;
    if isR
        X = XX;
        return
    end
end
ind1 = ind1(V(ind1)>0);
X = ind1(1);
end

function  [isR,tests,SC] = isRemovable_G(X,Mb,m,tests,SC,D,alpha)
isR = true;
Mb_X = find(Mb(X,:));
mbs = length(Mb_X);
m = min(m,mbs);
if mbs<=1
    return
end

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
for s = 1:m-2
    for j = 1:mbs-1
        Y = Mb_X(j);
        for k = j+1:mbs
            Z = Mb_X(k);
            subsets = subsets1(Mb_X([1:(j-1) (j+1):(k-1) (k+1):mbs]), mbs-s-2);
            for w=1:length(subsets)
                S = subsets{w};
                tests = tests+1;
                ls = length(S)+1;
                SC(ls+1)= SC(ls+1)+1;
                CI=CI_Test(Y,Z,[S X],D,0.7);
                if CI
                    isR = false;
                    return
                end
            end
        end
    end
    for j = 1:mbs
        Y = Mb_X(j);
        subsets = subsets1(Mb_X([1:(j-1), (j+1):mbs]),mbs-s-1);
        for w=1:length(subsets)
            S = subsets{w};
            ls = length(S);
            SC(ls+1)= SC(ls+1)+1;
            tests = tests+1;
            CI = CI_Test(Y,X,S,D, 0.7);
            if CI
                isR = false;
                return
            end
        end
    end
end
end

function [N,tests,SC] = FindNeighborhood(X,Mb_X,m,tests,SC,D,alpha)
N = Mb_X;
Mb_X = find(Mb_X);
mbs = length(Mb_X);
if mbs < m
    return;
end
if mbs<=1
    return
end
for j=1:mbs
    Y = Mb_X(j);
    sets = subsets1(Mb_X([1:(j-1), (j+1):mbs]), mbs-m);
    for i = 1:length(sets)
        S = sets{i};
        CI=CI_Test(Y,X,S,D,alpha);
        ls = length(S);
        SC(ls+1)=SC(ls+1)+1;
        tests=tests + 1;
        if CI
            N(Y) = false;
        end
    end
end
end

function [Mb,tests,SC,flag]= UpdateMb(X,N,Mb,tests,SC,flag,D,alpha_Mb)

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
        CI = CI_Test(Y,Z,S,D,alpha_Mb);
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

% ************************ Minor Functions ********************************
function C = mysetdiff(A,B)
% MYSETDIFF Set difference of two sets of positive integers (much faster than built-in setdiff)
% C = mysetdiff(A,B)
% C = A \ B = { things in A that are not in B }
%
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
function sub_s=subsets1(s,k)
% SUBSETS1 creates sub-sets of a specific from a given set
% SS = subsets1(S, k)
%
% S is the given set
% k is the required sub-sets size
%
% Example:
%
% >> ss=subsets1([1:4],3);
% >> ss{:}
% ans =
%      1     2     3
% ans =
%      1     2     4
% ans =
%      1     3     4
% ans =
%      2     3     4
%
% Written by Raanan Yehezkel, 2004

if k<0 % special case
    error('subset size must be positive');
elseif k==0 % special case
    sub_s={[]};
else
    l=length(s);
    ss={};
    if l>=k
        if k==1 % Exit condition
            for I=1:l
                ss{I}=s(I);
            end
        else
            for I=1:l
                ss1=subsets1(s([(I+1):l]),k-1);
                for J=1:length(ss1)
                    ss{end+1}=[s(I),ss1{J}];
                end
            end
        end
    end
    sub_s=ss;
end
end
