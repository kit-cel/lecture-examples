% use files from previous chapter
addpath('../ch5_LDPC_Finite_Length');

% construct an irregular matrix
% degree distribtion
% variables node of degrees 2 3 and 6. nodes_dv contains the factors
% Lambda_i, i.e., we have 6 variable nodes of degree 2, 12 variable nodes
% of degree 3 and 12 variable nodes of degree 6
dv = [2 3 6];
nodes_dv = [6 12 12];

% 15 check nodes of degree 8
dc = 8;
nodes_dc = [15];

% normalize degree distributions
nodes_dv_norm = nodes_dv / sum(nodes_dv);
nodes_dc_norm = nodes_dc / sum(nodes_dc);

% get polynomials lambda and rho
lambda = node_to_edge(dv, nodes_dv_norm);
rho = node_to_edge(dc, nodes_dc_norm);


H = getIrregularH(sum(nodes_dv), lambda, rho);
spy(H);

% generate random information word
u = randi(2,1,15)-1;

% uses a slightly different representation using an upper triangular matrix
% T than in the lecture. Conceptually identical though
x = Encode_LDPC(H, u);

% verify that x is a codeword
mod(H*x,2)