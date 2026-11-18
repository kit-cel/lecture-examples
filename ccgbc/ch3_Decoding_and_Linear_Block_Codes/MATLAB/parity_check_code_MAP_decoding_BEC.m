% This code is provided as supplementary material of the lecture Channel Coding - Graph Based Codes (CC-GBC)
%
% This code illustrates
%
%     MAP/ML decoding of a random parity-check code on the BEC

%% Code generation
% size of the matrix
m = 12;
n = 24;

% erasure probablity of BEC
epsilon = 0.4;

while(1)
    % construct parity-check matrix
    H = randi(2, m, n) - 1;

    % for encoding, split matrix into two parts
    % We have a systematic code x = (u p), and we split the parity-check matrix
    % into two parts H = (H_u H_p)
    % The condition H*x^T = 0 reads then H_u*u^T + H_p*p^T = 0
    % which we can solve for p^T = inv(H_p)*H_u*u^T
    % This means that H_p must be invertible, which we check first
    
    Hu = H(1:m, 1:m);
    Hp = H(1:m, (m+1):end);
    
    % abort construction if H_p is of full rank
    if rank(gf(Hp)) == m
        break;
    end
end

%% Encoding and Decoding using the given code
% encode
k = n-m;

% random information word
u = randi(2,1,k)-1;

% parity bits
p = inv(gf(Hp)) * gf(Hu) * gf(u)';
p = double(p.x)';

% codeword
x = [u, p];

% verify if codeword
if ~all(mod(H*x', 2) == 0)
    warning('Not a codeword!');
end


% transmission over BEC (-1 marks erasure)
y = x;
y(rand(size(x)) < epsilon) = -1;

% erasure index set
E = find(y == -1);
Ebar = setdiff(1:n, E);

% get sub-matrices
HE = H(:,E);
HEbar = H(:,Ebar);

% calculate s
s = mod(HEbar * y(Ebar)', 2);

% can we solve the system of equations?
r = rank(gf(HE));

if r == numel(E)  % full rank? 
    % if full rank, solve directly using Gaussian elimination
    xE = gflineq(HE, s);
    
    % fill missing positions in y
    y(E) = xE;
    
    fprintf('Decoding successful: Transmitted and decoded word:\n');
    x
    y
else
    fprintf('H_E is not of full rank  ... trying to generate the set of all compatible solutions. This will take some time for large n\n');
    xEposs = de2bi(0:(2^numel(E)-1),numel(E),'left-msb');
    idx_compat = find(sum(mod(HE*xEposs'+repmat(s,1,size(xEposs,1)),2),1) == 0);
    
    % all compatible solutions
    X_compat = xEposs(idx_compat,:);
    
    % as there are more than one compatible solution, carry out bitwise MAP
    % decoding
    % find positions that are equal in the set of compatible codewords
    idx_bitwise = all(~diff(X_compat),1);
    
    % replace erased positions with recovered bits
    xh_MAP = y;
    xh_MAP(E(idx_bitwise)) = X_compat(1,idx_bitwise);
    
    fprintf('Transmitted:\n');
    x
    fprintf('Received:\n');
    y
    fprintf('Bitwise MAP decoding result:\n');
    xh_MAP    
end
fprintf('----------------------------------------------------------\n');