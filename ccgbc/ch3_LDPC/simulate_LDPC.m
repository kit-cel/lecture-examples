function simulate_LDPC
% parameters of regular LDPC code
dv = 3;
dc = 6;

% specify Es/N0 at which simulation takes place
esno_dB = 2;

% number of frames to simulate
frames = 100;


% decodsing iterations
iterations = 5;

% compute noise standard deviation
sigma = sqrt(0.5 * 10.^(-esno_dB/10));

% channel parameter for LLR calculation
Lc = 4*10^(esno_dB/10);




% generate parity-check matrix of regular LDPC code
H = generate_Gallager(dv, dc, 30000);

n = size(H,2);

% simulate all-zero codeword
x = ones(1,n);

errors = 0;
errors_ref = 0;
for frame = 1:frames
    y = x + sigma*randn(1,n);

    % calculate LLRs
    L = Lc*y;

    xh = decode_LDPC(L, H, iterations);
    
   
    errors = errors + sum(xh ~= 0);
    
    % abort if positive 95% confidence interval of estimated BER within
    % 1/100th of estimated BER (disabled, can be uncommented for faster
    % speed when simulating good channels
    %BER = errors / frame / n;    
    %sigma_BER = sqrt(BER*(1-BER)/frame/n);
    %if BER > 0 && 2*sigma_BER < BER/100
    %    break;
    %end
end
BER = errors / frame / n;
fprintf('Es/N0 = %1.2f: BER = %1.4g\n',esno_dB, BER);
end
				
% LDPC decoder using the full (simplified) update rule, inner loop
% completely vectorized but still relatively slow
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


