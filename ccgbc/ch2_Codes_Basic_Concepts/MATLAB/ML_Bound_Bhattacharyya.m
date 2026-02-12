% This code is provided as supplementary material of the lecture Channel Coding - Graph Based Codes (CC-GBC)
%
% This code illustrates
%
%     Calculating the Bhattacharryya bound on compare it with ML decoding


% specify code by its generator matrix
G = [1 0 0 0 1 1 0;
     0 1 0 0 1 0 1;
     0 0 1 0 0 1 1];

% get parameters k and n of the code
k = size(G,1);
n = size(G,2);

% generate table of all codewords
% GF(2) operations are obtained by calculating modulo-2
u = de2bi(0:(2^k-1),k,'left-msb');  % all possible input patterns
code = mod(u * G, 2);  % generate codewords


% compute weights
w = sum(code,2);
dmin = min(w(2:end));  % exclude all-zero codeword

% weight enumerator polynomial
A = [fliplr(histc(w, 1:n)'), 0];



% channel
esno_dB = 5; % specify Es/N0 in dB

esno_lin = 10^(esno_dB/10);


% Bhattacharyya parameter
B = exp(-esno_lin);

% compute bound on ML performance
P_ML_Bound = polyval(A, B);

fprintf('According to the Bhattacharyya parameter, the ML error probability is upper bounded by %1.4g\n',P_ML_Bound);



% calculate Monte-Carlo-Estimate of ML performance
% transmit N codewords
N = 1000000;

% standard deviation of noise distribution
sigma_n = sqrt(1/2/esno_lin);

errors = 0;
for j = 1:N
    % select random codeword
    idx = randi(size(code,1));
    c = code(idx, :);
    
    % after channel
    y = (1-2*c) + sigma_n*randn(size(c));
    
    % calculate correlation with all codewords
    correlation = sum((1-2*code).*repmat(y, size(code,1), 1), 2);
    
    [~,ML_decision_idx] = max(correlation);
    
    errors = errors + (ML_decision_idx ~= idx);
end
P_ML = errors/N;
fprintf('After ML decoding, the error rate is approximately %1.4g\n', P_ML);
