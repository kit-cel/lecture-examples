% This code is provided as supplementary material of the lecture Channel Coding: Algebraic Methods
% given by Laurent Schmalen at Karlsruhe Institute of Technology (KIT)
%
% This code illustrates:
*
% * Erasure recovery of GRS codes using simple Gaussian elimination
%
% This code (and other examples) is available on https://github.com/kit-cel/lecture-examples
%

%% construct GF(2^4) using primitive polynomial Z^4 + Z^3 + 1
% primitive element alpha
alpha = gf(2,4,'D^4+D^3+1');

m = alpha.m;
q = 2^m;

% array with all q elements of the Galois field
gf_elements = gf(0:(q-1), 4, alpha.prim_poly);




% correctable erasures of the RS code
t = 6;



% Construct primitive GRS code (i.e., q = n-1)
n = q-1;
k = n - t;  % k is immediately given as the GRS code is MDS (k = n-d+1 and d = t+1)

% use all q-1 disting elements of GF(q) as code locators
gamma = alpha.^[0:(n-1)];


% use random column multipliers
v = alpha.^[randi(n,1,n)-1];


% construct parity-check matrix of the code
H = [];
for j = [0:(n-k-1)]
    H = [H; gamma.^j];
end
H = H * diag(v);




% column multipliers of the dual code are given by gamma/v, as our GRS code
% is primitive
v_prime = gamma ./ v;

% construct generator matrix (parity-check matrix of the dual code)
G = [];
for j = [0:(k-1)]
    G = [G; gamma.^j];
end
G = G * diag(v_prime);


% check if G*H' is actually zero
if any(G*H' ~= 0)
    error('Something went wrong during code construction');
end

%% encode a randomly generated codeword
u = gf_elements(randi(q,1,k));

% encode
x = u*G;

% generate erasures patterns
% number of erasures
E = 6;
if E > t
    warning('number of errors exceeds the correction capabilities of the code');
end

% set of erased positions
E = randperm(n,E);

% generate received codeword
y = x;

% erase positions (set to zero in received codeword)
y(E) = gf(0, m, alpha.prim_poly);



%% carry out decoding
% setup system of equation

% sets containing erased and non-erased position
idx_erased = sort(E);
idx_nonerased = setdiff(1:n, E);

% setup system of equations
A = [];
for i = 1:numel(idx_erased)
    A = [A; v(idx_erased).*gamma(idx_erased).^(i-1)];
end

% right hand side
b = [];
for i = 1:numel(idx_erased)
    b = [b; (v(idx_nonerased).*gamma(idx_nonerased).^(i-1)) * y(idx_nonerased)'];
end

% solve system of equations
erased_values = A \ b;

y(idx_erased) = erased_values;


if all(y == x)
    fprintf('Decoding successful!\n');
end

