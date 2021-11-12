% size of matrix
m = 10;
n = 20;

% construct random matrix
H = randi(2, m, n) - 1;

% show rank and maximum possible rank
fprintf('Rank of matrix is %d (maximum rank: %d)\n', rank(gf(H)), min([m,n]));
