function isD = isDiamondfree(G)

isD = true;
n=size(G,1);
H = G+G';
for X=1:n
    Pa = find(G(:,X));
    Pas = length(Pa);
    if Pas >= 3
        for i=1:Pas
            Y = Pa(i);
            N_Y = find(H(:,Y));
            S = myintersect(Pa,N_Y);
            ls = length(S);
            if nnz((H(S,S)+eye(ls))==0)>0
                isD = false;
                return
            end
        end
    end
end
end

function C = myintersect(A,B)
% MYINTERSECT Intersection of two sets of positive integers
% (much faster than built-in intersect)
% C = myintersect(A,B)

A = A(:)'; B = B(:)';

if isempty(A)
    ma = 0;
else
    ma = max(A);
end

if isempty(B)
    mb = 0;
else
    mb = max(B);
end

if ma==0 || mb==0
    C = [];
else
    %bits = sparse(1, max(ma,mb));
    bits = zeros(1, max(ma,mb));
    bits(A) = 1;
    C = B(logical(bits(B)));
end
end