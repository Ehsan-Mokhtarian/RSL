function [A,D] = Generate_Data(graph_type,number_of_samples,n,p)

if nargin<3
    A = Generate_Graph(graph_type);
    n = size(A,1);
else
    A = Generate_Graph(graph_type,n,p);
end
N = randn(number_of_samples,n)*diag(0.7 + 0.52*rand(1,n));
AA = (1+0.5*rand(n)) .* ((-1).^(rand(n)>0.5));
AA = A.*AA;
D = N/(eye(n)-AA);
end
