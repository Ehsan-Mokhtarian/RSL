function [H,tests,SC,flag] = RSL_D_oracle(G, Mb, varargin)
% input:
% G: The true DAG. It is only used in dsep function to perform oracle CI
% tests
% Mb: n*n Markov boundary matrix. Entries can be 0 or 1
%
% output:
% H: The matrix of learned skeleton. Entries can be 0 or 1
% tests: the number of performed CI tests
% SC: Size of Conditioning sets. It is a (n+1)*1 vector such that i-th 
% entry indicates the number of conditioning sets with size i-1


%**************************************************************************
%********************** Start of the Algorithm ****************************
%**************************************************************************

n = size(Mb,1);
if nargin > 2
    V = varargin{1};
    tests = varargin{2};
    SC = varargin{3};
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
    [X,tests,SC,flag] = FindRemovable(V,G,Mb,tests,SC,flag);
    X
    [N,tests,SC] = FindNeighborhood(X,Mb(:,X),G,tests,SC);
    [Mb,tests,SC,flag]= UpdateMb(X,N,Mb,G,tests,SC,flag);
    V(X)=0;
    [H,tests,SC,flag] = RSL_D_oracle(G,Mb,V,tests,SC,flag);
    H(X,N==1)=1;
    H(N==1,X)=1;
    return
end
end

%**************************************************************************
%************************* Main Functions *********************************
%**************************************************************************

function [X,tests,SC,flag] = FindRemovable(V,G,Mb,tests,SC,flag)
Mbs = sum(Mb);
[~,ind] = sort(Mbs);
ind1 = ind(flag(ind)>0);
for XX=ind1
    Mb_XX = find(Mb(XX,:));
    [isR,tests,SC] = isRemovable(XX,Mb_XX,G,tests,SC);
    flag(XX) = false;
    if isR
        X = XX;
        return
    end
end
ind2 = ind(V(ind)>0);
X = ind2(1);
end

function  [isR,tests,SC] = isRemovable(X,Mb_X,G,tests,SC)

isR = true;
mbs = length(Mb_X);

for j =1:(mbs-1)
    Y = Mb_X(j);
    for k = (j+1):mbs
        Z = Mb_X(k);
        S = [X, Mb_X([1:(j-1), (j+1):(k-1), (k+1:mbs)])];
        CI = dsep(G,Y,Z,S);
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

function [N,tests,SC] = FindNeighborhood(X,Mb_X,G,tests,SC)

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
        if dsep(G, X, Y, S)
            N(Y)=0;
            break
        end
    end
end

end

function [Mb,tests,SC,flag]= UpdateMb(X,N,Mb,G,tests,SC,flag)

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
        CI = dsep(G,Y,Z,S);
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
