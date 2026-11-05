function simulate_LDPC
% parameters of regular LDPC code
dv = 3;
dc = 6;

% specify Es/N0 at which simulation takes place
esno_dB = 1;

% number of frames to simulate
frames = 1000;


% decodsing iterations
iterations = 5;

% compute noise standard deviation
sigma = sqrt(0.5 * 10.^(-esno_dB/10));

% channel parameter for LLR calculation
Lc = 4*10^(esno_dB/10);




% generate parity-check matrix of regular LDPC code
H = generate_Gallager(dv, dc, 1200);

n = size(H,2);

% simulate all-zero codeword
x = ones(1,n);

errors = 0;
for frame = 1:frames
    y = x + sigma*randn(1,n);

    % calculate LLRs
    L = Lc*y;

    xh = decode_LDPC(L, H, iterations);
    
    % alternative using the full matrix (faster for very small matrices,
    % slower for large matrices)
    %xh = decode_LDPC_nosparse(L, full(H), iterations);
   
    errors = errors + sum(xh ~= 0);
end
BER = errors / frame / n;
fprintf('Es/N0 = %1.2f: BER = %1.4g\n',esno_dB, BER);
end
				
% LDPC decoder using the full (simplified) update rule, inner loop
% completely vectorized but still relatively slow
% based on a sparse matrix as input
function xh = decode_LDPC(L, H, iterations)
    [row_i, col_i] = find(H);  % get row and column indices
    
    % initialize variable to check node messages with channel output
    VtoC = sparse(row_i, col_i, L(col_i));

    % main iterations
    for it = 1:iterations
        % compute check to variable sum(CtoV,1)node messages
        VtoC_sign = spfun(@sign, VtoC);
        VtoC_abs = spfun(@abs, VtoC);
        
        phiVtoC = spfun(@(x)(log(coth(x/2))), VtoC_abs);
        phiVtoC_sum = sum(phiVtoC,2);
        
        [tri,~,values] = find(VtoC_sign);
        totalsign_VtoC = accumarray(tri,values,[],@prod);
       
        CtoV_abs = spfun(@(x)(log(coth(x/2))), sparse(row_i, col_i,phiVtoC_sum(row_i)) - phiVtoC);
        CtoV_sign = sparse(row_i, col_i, totalsign_VtoC(row_i)) .* VtoC_sign;
        CtoV = CtoV_sign .* CtoV_abs;        
       
        
        % compute variable to check node messages, pretty simple
        CtoV_sum = sum(CtoV,1);
        VtoC = sparse(row_i, col_i, L(col_i) + CtoV_sum(col_i)) - CtoV;        

        % stopping criterion, all parity checks are fulfilled
        L_total = CtoV_sum + L;

        % binary decision
        xh = L_total < 0;
        if all(mod(H*xh(:),2) == 0)
            % all parity-checks fulfilled?
            break;
        end       
    end
    
end




% LDPC decoder using the full (simplified) update rule, inner loop
% completely vectorized but still relatively slow
% based on a non-sparse parity-check matrix as input
function xh = decode_LDPC_nosparse(L, H, iterations)
    n = size(H,2);
    m = size(H,1);
    
    [row_i, col_i] = find(H);  % get row and column indices
    
    % initialize variable to check node messages with channel output
    VtoC = repmat(L(:)', m, 1) .* H;
    
    % main iterations
    for it = 1:iterations
        % compute check to variable sum(CtoV,1)node messages
        VtoC_sign = sign(VtoC);
        VtoC_abs = abs(VtoC);
        
        phiVtoC = log(coth(VtoC_abs/2 + (1-H).*9e9));  % mask out zero entries
        phiVtoC_sum = sum(phiVtoC,2);
        
        % multiply signs
        [tri,~,values] = find(VtoC_sign);
        totalsign_VtoC = accumarray(tri,values,[],@prod);
       
        CtoV_abs =  log(coth((repmat(phiVtoC_sum, 1, n).*H - phiVtoC)/2 + (1-H).*9e9));
        CtoV_sign = repmat(totalsign_VtoC, 1, n) .* VtoC_sign;
        CtoV = CtoV_sign .* CtoV_abs;        
       
        
        % compute variable to check node messages, pretty simple
        CtoV_sum = sum(CtoV,1);
        VtoC = (repmat(CtoV_sum + L(:)', m, 1) - CtoV) .* H;        

        % stopping criterion, all parity checks are fulfilled
        L_total = CtoV_sum + L;

        % binary decision
        xh = L_total < 0;
        if all(mod(H*xh(:),2) == 0)
            % all parity-checks fulfilled?
            break;
        end       
    end
    
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


