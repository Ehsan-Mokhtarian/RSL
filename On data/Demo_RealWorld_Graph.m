close all
clear
clc

addpath('structures/')
addpath('functions/')

%% Load graph and generate data
% See more graphs in 'structures' folder
graph_name = 'insurance';
number_of_samples = 10000;
[G,D] = Generate_Data(graph_name,number_of_samples);
n = size(G,1); % number of vertives
V = ones(n,1);
m = clique_number(G); % clique number of the graph
if isDiamondfree(G)
    disp("G is diamond-free")
else
    disp("G is not diamond-free")
end
fprintf("Clique number: %d \n\n",m);

%% Markov boundary discovery
alpha_Mb = 2/n^2;
Mb = ComputeMb_TC(D, alpha_Mb); % Learning  Markov boundaries using TC algorithm

%% Calling RSL functions
alpha = 0.01;
[G_RSLD, tests_RSLD, SC_RSLD] = RSL_D(D, Mb, alpha, alpha_Mb);
[G_RSLW, tests_RSLW, SC_RSLW] = RSL_W(D, V, Mb, alpha, alpha_Mb, m);

%% Accuracy of the learned graph
report_accuracy("RSL_D",G,G_RSLD,tests_RSLD);
report_accuracy("RSL_W",G,G_RSLW,tests_RSLW);