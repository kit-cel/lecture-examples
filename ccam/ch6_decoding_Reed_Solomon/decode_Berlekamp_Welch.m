%% This code is provided as supplementary material of the lecture Channel Coding: Algebraic Methods
% given by Laurent Schmalen at Karlsruhe Institute of Technology (KIT)
%
% This code illustrates:
%
% * Decoding of GRS codes using the Berlekamp-Welch decoding algorithm
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




% correctable errors of the RS code
t = 3;



% Construct primitive GRS code (i.e., q = n-1)
n = q-1;
k = n - 2*t;  % k is immediately given as the GRS code is MDS (k = n-d+1 and d = 2*t+1)

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

% generate error patterns
% number of errors
E = 3;
if E > t
    warning('number of errors exceeds the correction capabilities of the code');
end

% set of error positions
J = randperm(n,E);

% construct error vector (random elements)
e = gf(zeros(1,n),m,alpha.prim_poly);
e(J) = gf_elements(randi(q,1,E));

% generate noisy codeword
y = x + e;

%% carry out decoding
% Initial step, remove influence of column multiplier
z = y ./ v_prime;

% First step, setup system of equations 
EQ_system = gf([], m, alpha.prim_poly);
for i = 1:n
    EQ_system(i,:) = [z(i)*gamma(i).^[0:(t-1)], gamma(i).^[0:(t+k-1)]];
end

if rank(EQ_system) == n
    % solve system of equations directly as it is full rank
    EQ_sol = EQ_system \ (z(:) .* gamma(:).^t);
else
    % solve using Gaussian elimination, assume remaining entries of the
    % polynomials are zero
    [A, idx] = gfrref([EQ_system, (z(:).*gamma(:).^t)]);
    EQ_sol = gf(zeros(n,1), m, alpha.prim_poly);
    EQ_sol(idx) = A(1:numel(idx),end);
end

% extract polynomials E and Q
E_poly = [1,fliplr(EQ_sol(1:t)')];
Q_poly = fliplr(EQ_sol((t+1):end)');

% remove leading non-zero elements for E, if we have less than t errors
first_nonzero = find(E_poly.x > 0, 1);
E_poly = E_poly(first_nonzero:end);

% carry out polynomial devision
[info_poly_hat, ~] = deconv(Q_poly, E_poly);

% recover codeword and information word
x_hat = gf(zeros(1,n), m, alpha.prim_poly);
for i = 1:n
    x_hat(i) = polyval(info_poly_hat, gamma(i)) * v_prime(i);
end

u_hat = gf(zeros(1,k), m, alpha.prim_poly);
for i = 1:k
    u_hat(i) = info_poly_hat(end-i+1);
end


if all(x_hat == x)
    fprintf('Decoding successful!\n');
else
    fprintf('ERROR!\n');
end

