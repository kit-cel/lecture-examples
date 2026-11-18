% This code is provided as supplementary material of the lecture Channel Coding - Graph Based Codes (CC-GBC)
%
% This code illustrates
%
%     Rank of random parity-check matrix


% size of matrix
m = 10;
n = 20;

% construct random matrix
H = randi(2, m, n) - 1;

% show rank and maximum possible rank
fprintf('Rank of matrix is %d (maximum rank: %d)\n', rank(gf(H)), min([m,n]));
