function RS_example()
%RS_EXAMPLE Example run. Create an RS codes, encode a message, introduce
%errors, decode
%   Detailed explanation goes here
    
    % Basic Parameters
    q = 2^4; % Field Size
    n = q - 1; % Blocklength
    t = 3; % Correctable errors
    k = n - 2*t; % Information length
    
    % Sequence of zeros
    alphas = 1:2*t; 

    % Construct generator polynomial
    g = gf(1,log2(q));
    a = gf(2,log2(q)); % 2 is always a primitive element in MATLABs gf implementation
    for i = 1:length(alphas)
        % Convolution is polynomial multiplication
        g = conv(g,[1 -a^alphas(i)]);
    end
    
    % Create random message and encode
    m = randi([0,q-1],1,k);
    c = RSenc(m,g,q);
    
    % Corrupt randomly 
    n_e = 3; % No. of Errors if n_e <= t the erros should be corrected
    E_pos = randperm(n,n_e);

    r = c;
    r(E_pos) = randi([0,q-1],1,n_e); % Set values at E_pos to trandom symbols

    % Decode
    c_hat = RSdec(r,q,t);

    % Check for Errors
    errs = sum(c_hat~=c);
    fprintf('%d wrong symbols \n',errs);
end