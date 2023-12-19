function simulate_LDPC_BEC
% parameters of regular LDPC code
dv = 3;
dc = 6;

% specify epsilon (erasure probability) at which simulation takes place
epsilon = 0.2;

% number of frames to simulate
frames = 100;


% decoding iterations
iterations = 50;

% generate parity-check matrix of regular LDPC code
H = generate_Gallager(dv, dc, 48);

n = size(H,2);

% simulate all-zero codeword
x = zeros(1,n);

errors = 0;
for frame = 1:frames
    % erasure channel, first map to bipolar and map to very large value as
    % approximation to infinite LLR
    y = (1 - 2*x) * 9999;           
    y(rand(size(x)) < epsilon) = 0;  % erasures (LLR of zero)
    
    xh = decode_LDPC_BEC_nosparse(y, full(H), iterations);
    
    errors = errors + isempty(xh);
end
FER = errors / frames;    % divide by two, as we may correctly guess the residual erasures
fprintf('epsilon = %1.2f: FER = %1.4g\n',epsilon, FER);
end



% LDPC decoder using the full (simplified) update rule, inner loop
% completely vectorized but still relatively slow
% based on a non-sparse parity-check matrix as input
function xh = decode_LDPC_BEC_nosparse(L, H, iterations)
    n = size(H,2);
    m = size(H,1);
    
    [row_i, col_i] = find(H);  % get row and column indices
    
    % initialize variable to check node messages with channel output
    VtoC = repmat(L(:)', m, 1) .* H;
    
    % main iterations
    for it = 1:iterations
        % compute check to variable sum(CtoV,1)node messages
        VtoC_sign = mysign(VtoC);
        VtoC_abs = abs(VtoC);
        
        phiVtoC = phifun(VtoC_abs/2 + (1-H).*9e9);  % mask out zero entries
        phiVtoC_sum = sum(phiVtoC,2);
        
        % multiply signs
        [tri,~,values] = find(VtoC_sign);
        totalsign_VtoC = accumarray(tri,values,[],@prod);
       
        CtoV_abs =  phifun((repmat(phiVtoC_sum, 1, n).*H - phiVtoC)/2 + (1-H).*9e9);
        CtoV_sign = repmat(totalsign_VtoC, 1, n) .* VtoC_sign;
        CtoV = CtoV_sign .* CtoV_abs;        
       
        
        % compute variable to check node messages, pretty simple
        CtoV_sum = sum(CtoV,1);
        VtoC = (repmat(CtoV_sum + L(:)', m, 1) - CtoV) .* H;        

        % stopping criterion, all parity checks are fulfilled
        L_total = CtoV_sum + L;
        
        % binary decision
        if any(abs(L_total) < 1)
            % erasures left
            xh = [];
        else
            xh = L_total < 0;
        end       
    end
    
end




function y = phifun(x)
    y = log(coth(x));
    y(isinf(y)) = 9e9;
end




% generate a parity-check matrix according to Gallager's method
% do not care about 4-cycles
function H = generate_Gallager(dv, dc, n)
    if mod(n,dc) ~= 0
        error('n must be a multiple of check node degree dc');
    end
    rows = floor(n / dc);
    % column indices
    jj = 1:n;
    ii = reshape(repmat([1:rows],dc,1), 1, []);
    Ho = sparse(ii,jj,ones(size(jj)),rows,n);
    H = Ho;
    for k=1:(dv-1)
        H = [H; Ho(:,randperm(n))];
    end
end    

function y = mysign(x)
    y = ones(size(x));
    y(x < 0) = -1;
end
