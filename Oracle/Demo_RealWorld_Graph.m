close all
clear
clc

addpath('structures/')
addpath('functions/')

%% Load graph and generate data
% See more graphs in 'structures' folder
graph_name = 'diabetes';
G = Generate_Graph(graph_name);
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
Mb = ComputeMb_oracle(G);

%% Calling RSL functions
[G_RSLD, tests_RSLD, SC_RSLD] = RSL_D_oracle(G, Mb);
[G_RSLW, tests_RSLW, SC_RSLW] = RSL_W_oracle(G,V,Mb,m);

%% Accuracy of the learned graph
report_accuracy("RSL_D",G,G_RSLD,tests_RSLD);
report_accuracy("RSL_W",G,G_RSLW,tests_RSLW);