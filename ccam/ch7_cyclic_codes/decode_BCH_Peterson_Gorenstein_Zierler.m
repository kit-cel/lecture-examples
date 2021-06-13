%% construct GF(2^4) using primitive polynomial Z^4 + Z^3 + 1
% primitive element alpha
alpha = gf(2,4,'D^4+D^3+1');

m = alpha.m;
q = 2^m;

% array with all q elements of the Galois field
gf_elements = gf(0:(q-1), 4, alpha.prim_poly);




% correctable errors of the BCH code
t = 3;

% length of code
n = q-1;


% obtain table of minimal polynomials
Mpoly_list = [];
for j=1:(2*t)
    M = minpol(alpha^j);
    % store in list and convert to GF(2)
    Mpoly_list(end+1,:) = [zeros(1,n-prod(size(M))), double(M.x)];
end
% lcm operation corresponds to unique rows of the list of minimal
% polynomials
Mpoly_list = unique(Mpoly_list, 'rows');

% compute generator polynomial by multiplying minimal polynomials
% (corresponds tp convolution). Attention, operations are in GF(2)
gen_poly = [1];
for j=1:size(Mpoly_list)
    gen_poly = mod(conv(gen_poly, Mpoly_list(j,:)), 2);
end
% remove leading zeros and flip, such that the g_0 is the first entry
% gen_poly(1)
gen_poly = fliplr(gen_poly(find(gen_poly > 0, 1, 'first'):end));

% dimension of BCH code is obtained from degree of polynomial
k = n - (numel(gen_poly) - 1);  % number of elements in polynomial vector is degree+1




%% encode a randomly generated codeword
u = randi(2,1,k) - 1;

% systematic encoding, compute remainder after polynomial division of Z^(n-k) *
% u(Z) by g(Z)
% ATTENTION: gfdeconv uses a different representation of polynomials
[q,r] = gfdeconv([zeros(1,n-k), u], gen_poly);

x = [r,zeros(1,n-k-numel(r)), u];


% cross check if x is actually a codeword
% evaluating the polynomial x at each alpha^j must be zero
for j=1:2*t
    temp = gf(0, m, alpha.prim_poly);
    for l = 1:n
        temp = temp + x(l) * alpha^(j * (l-1));
    end
    if temp ~= 0
        error('x is not a codeword of the BCH code');
    end
end


% generate error patterns
% number of errors
E = 3;
if E > t
    warning('number of errors exceeds the correction capabilities of the code');
end

% set of error positions
J = randperm(n,E);

% construct error vector (random elements)
e = zeros(1,n);
e(J) = ones(1,E);

% generate noisy codeword
y = mod(x + e, 2);



%% carry out decoding
% First step: compute syndrome
S = gf(zeros(1,2*t), m, alpha.prim_poly);
for j=1:2*t
    temp = gf(0, m, alpha.prim_poly);
    for l = 1:n
        temp = temp + y(l) * alpha^(j * (l-1));
    end
    S(j) = temp;
end

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
        
        % position of error is immediately obtained from error locator,
        % attention, we need to flip the positions
        error_pos(end+1) = 1+log(error_locator);
    end
end

% as the code is binary, we do not need the error values

% reconstruct codeword
e_hat = zeros(1,n);
e_hat(error_pos) = 1;

x_hat = mod(y + e_hat, 2);

if all(x_hat == x)
    fprintf('Decoding successful!\n');
end



