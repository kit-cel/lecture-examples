function simulate_LDPC_BEC
% parameters of regular LDPC code
dv = 3;
dc = 6;

% size of Code
n = 200;

% specify epsilon (erasure probability) at which simulation takes place
epsilon = 0.2;

% number of frames to simulate
frames = 5000;


% decoding iterations
iterations = 200;


H = getRegularH(n, dv, dc);

% simulate all-zero codeword
x = zeros(1,n);

errors = 0;
h = waitbar(0,'Simulating ...','CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
setappdata(h, 'canceling', 0);
for frame = 1:frames
    if mod(frame,20) == 0
        % generate parity-check matrix of regular LDPC code on a regular
        % basis so that we average over the parity-check matrices
        H = getRegularH(n, dv, dc);
    end
    % erasure channel, first map to bipolar and map to very large value as
    % approximation to infinite LLR
    y = 1 - 2*x;        
    y(rand(size(x)) < epsilon) = 2;  % 2 denotes erasure
    
    xh = decode_LDPC_BEC(y, H, iterations);
    
    errors = errors + isempty(xh);
    
    if getappdata(h,'canceling')
        break
    end
    if mod(frame,20) == 0
        waitbar(frame/frames, h, sprintf('Simulating ... current estimate P_E = %1.2e', errors/frame));
    end
end
delete(h);
FER = errors / frame    % divide by two, as we may correctly guess the residual erasures
fprintf('epsilon = %1.2f: FER = %1.4g\n',epsilon, FER);
end


% variable node 
function y = MP_BEC_VN(x, yc)
    y = yc*ones(size(x));
    if yc == 2 % channel is erased
        if sum(x == 2) == numel(x)-1  % all erased but one message
            y(x == 2) = x(find(x < 2,1));
        elseif sum(x==2) < numel(x)-1 % at least two non-erased messages
            y = x(find(x < 2,1))*ones(size(x));
        end
    end
end

% check node update function
function y = MP_BEC_CN(x)
    y = 2*ones(size(x));

    % Case 1, no erasure
    if sum(x == 2) == 0
        y = x.*prod(x);
    % Case 2, one single erasures
    elseif sum(x == 2) == 1
        y(x==2) = prod(x(x ~= 2));
    end
end
				
% LDPC decoder for the BEC, inner loop
% completely vectorized but still relatively slow
% based on a sparse matrix as input
function xh = decode_LDPC_BEC(y, H, iterations)
    [row_i, col_i] = find(H);  % get row and column indices
    
    % initialize variable to check node messages with channel output
    VtoC = sparse(row_i, col_i, y(col_i));

    % initialize check to variable node messages with erasures
    CtoV = sparse(row_i, col_i, 2*ones(size(col_i)));

    xh = [];
    % main iterations   
    for it = 1:iterations

        % iterate through all check nodes        
        for j = 1:size(H,1)
            CtoV(j,:) = spfun(@MP_BEC_CN, VtoC(j,:));
        end

        % iterate through all variable nodes
        for j = 1:size(H,2)
            VtoC(:,j) = spfun(@(x)(MP_BEC_VN(x, y(j))), CtoV(:,j));
        end



        % binary decision
        x_dec = y;
        min_col = accumarray(col_i, nonzeros(CtoV), [], @min);
        x_dec(y == 2) = min_col(y == 2);        
        if all(x_dec < 2)
            xh = (1 - x_dec)/2;            
            if all(mod(H*xh(:),2) == 0)
                % all parity-checks fulfilled?
                break;
            else
                xh = [];
            end
        end       
    end
end


