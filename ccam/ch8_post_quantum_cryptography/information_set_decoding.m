% parameters of a code
n = 20;
k = 10;
t = 3;


% generate a random parity-check matrix
G = randi(2,k,n)-1;

% generate a random information word
u = randi(2,1,k)-1;

% encode
x = mod(u*G,2);

% add t errors
e = zeros(1,n);
e(randperm(n,t)) = 1;

y = mod(x + e,2);

% information set decoding
trials = 0;

% try until we have found the information
while true
    trials = trials + 1;
    
    % select k random columns
    idx = randperm(n,k);
    G_prime = G(:,idx);
    
    % check if G' is invertible
    % only continue if G' is full rank
    if gfrank(G_prime) ~= k        
        continue
    end
    
    % extract code bits corresponding to G'
    y_prime = y(idx);
    
    % recover information part1
    G_prime_inv = inv(gf(G_prime));
    G_prime_inv = double(G_prime_inv.x);
    
    u_hat = mod(y_prime * G_prime_inv,2);
    
    if all(u == u_hat)
        break;
    end
end     


fprintf('Required number of trials until information was guessed: %d\n',trials);
