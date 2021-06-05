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
% First step: compute syndrome
S = y * H';

if all(S==0)
    fprintf('No errors, aborting');
    return;
end

% generate matrix Stilde (try all possible tau, starting from the largest)
for tau = t:-1:1
    
    Stilde = [];
    for j = 1:tau
        Stilde = [Stilde; S([tau:-1:1]+(j-1))];
    end
    sv = -S([(tau+1):2*tau])';

    % if matrix is invertible, invert
    if rank(Stilde) == tau
        Lambda = inv(Stilde)*sv;
        break;
    end
end

% construct error locator polynomial
Lambda_poly = [fliplr(Lambda'), 1];

error_pos = [];
% find roots of error locator polynomial, by trying out all possible nonzero GF
% elements
for j=1:n
    if polyval(Lambda_poly, gf(j,m,alpha.prim_poly)) == 0
        error_locator = 1/gf(j,m,alpha.prim_poly);
        % find code locator that corresponds to error locator
        error_pos(end+1) = find(gamma == error_locator,1);
    end
end




% find error values using Forney's method
% first, compute the error evaluator polynomial
Tm = gf(zeros(tau, tau+1), m, alpha.prim_poly);
for j=1:tau
    Tm(j, j:-1:1) = S(1:j);
end
Omega = Tm*fliplr(Lambda_poly)';

% generate polynomial
Omega_poly = fliplr(Omega');

% compute error values using Forney's algorithm
% calculate the derivative of the error locator polynomial. Attention must
% be paid that the multiplication with the exponent is actually carried out
% as additions in the Galois field
derLambda_poly = Lambda_poly(1:end-1) .* mod([(prod(size(Lambda_poly))-1):-1:1],2);

e_hatF = gf(zeros(1,n), m, alpha.prim_poly);
for j=1:numel(error_pos)
    polyval(derLambda_poly, 1/gamma(error_pos(j)))
    e_hatF(error_pos(j)) = - gamma(error_pos(j)) * polyval(Omega_poly, 1/gamma(error_pos(j))) / polyval(derLambda_poly, 1/gamma(error_pos(j))) / v(error_pos(j));
end

% reconstruct codeword
x_hatF = y - e_hatF;

if all(x_hatF == x)
    fprintf('Decoding using Forney''s algorithm successful!\n');
end







%% find error values (direct method, using matrix inversion, may run into issues as matrix is not invertible)
Em = [];
for j=1:prod(size(S))
    Em = [Em; v(error_pos) .* gamma(error_pos).^(j-1)];
end

% Solve system of equations using Gaussian elimination
error_val = Em \ S';

% estimate error vector
e_hat = gf(zeros(1,n), m, alpha.prim_poly);
e_hat(error_pos) = error_val;


% reconstruct codeword
x_hat = y - e_hat;

if all(x_hat == x)
    fprintf('Decoding successful!\n');
end

