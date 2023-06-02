% This code is provided as supplementary material of the lecture Channel Coding: Algebraic Methods
% given by Laurent Schmalen at Karlsruhe Institute of Technology (KIT)
%
% This code illustrates:
%
% * Construction, encoding and decoding of a two-error correcting code
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

%% construct 2-error correcting code
% generate parity-check matrix
H_gf = [alpha.^[0:14]; alpha.^([0:14]*3)];

% create binary parity-check matrix
H_bin = zeros(m*size(H_gf,1), size(H_gf,2));
for i=1:size(H_gf,2)
    temp = H_gf(1,i);
    H_bin(1:m,i) = de2bi(temp.x, m);
    temp = H_gf(2,i);
    H_bin(m+[1:m],i) = de2bi(temp.x, m);
end
% parity-check matrix part related to information bits
n = size(H_bin,2);
k = n - size(H_bin,1);
H_bin_u = gf(H_bin(:,1:(n-size(H_bin,1))), 1);  % binary matrix
H_bin_p = H_bin(:,(n-size(H_bin,1)+1):end);
H_bin_p_inv = inv(gf(H_bin_p,1));

% from H, we can calculate the generator matrix in the following way
% We assume a systematic code
% Then H_u corresponds to the first k columns of H and H_p to the last
% columns of H
% x is a codeword if H*x' = 0
% If x = [u p], that means that H_u*u + H_p*p = 0 must hold
% Hence we get p = -inv(H_p)*H_u*u
% By assuming u to be the vectors with a single 1, we get the generator
% matrix

% generate generator matrix
G = [eye(k), zeros(k,n-k)];

for i = 1:k
    u = zeros(1,k);
    u(i) = 1;
    p = H_bin_p_inv * H_bin_u * u';
    G(i,(k+1):end) = p.x;
end

%% encode information 
% random information vector
u = randi(2,1,k)-1;

% codeword
x = mod(u*G,2);
fprintf('Transmitted codeword: ['); fprintf('%d ', x(1:end-1)); fprintf('%d]\n',x(end));

%% add errors
errors = 3;

error_pos = randperm(n,errors);
x(error_pos) = mod(x(error_pos) + 1,2);
fprintf('Erroneous codeword:   ['); fprintf('%d ', x(1:end-1)); fprintf('%d]\n',x(end));

%% decoding
% Compute syndrome
S = H_gf*x(:);

% check if only one error
if S(1) == 0 && S(2) == 0
    fprintf('No errors detected\n');
else
    if S(1).^3 == S(2)
        error_locations = log(S(1)) + 1;
        fprintf('One error detected at position %d\n', error_locations)    
    else
        error_locations = [];
        % find roots of the polynomial x^2 + s1*x + s2/s1 + s1^2 = 0
        for i = 1:n
            tx = alpha^(i-1);
            if tx.^2 + S(1)*tx + S(2)/S(1) + S(1)^2 == 0
                error_locations(end+1) = i;
            end
        end
        if numel(error_locations) == 2
            fprintf('Two errors detected at positions %d and %d\n', error_locations(1), error_locations(2));
        else
            fprintf('Could not detect errors\n');
        end
    end
    x(error_locations) = mod(x(error_locations)+1,2);
end
fprintf('Corrected codeword:   ['); fprintf('%d ', x(1:end-1)); fprintf('%d]\n',x(end));

% 