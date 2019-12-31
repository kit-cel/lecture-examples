function LDPC_Decode_Image
%codefile = 'Code_Rate12_regular.alist';
codefile = 'Code_Rate0.63_irregular.alist';     % alternative file
H = read_alist_file(codefile);

% generate encoder object (takes some time to initialize)
henc = comm.LDPCEncoder(H);


% specify Es/N0 at which simulation takes place
esno_dB = 1;

% channel output can be 'soft' or 'hard'
channel_output = 'soft';


% load image, size 128x128
cat = imread('LogoCEL.png');

% transform to black and white image
cat(cat <= 127) = 0;
cat(cat > 127) = 1;
cat = 1-cat;

% encode
x = step(henc,double(reshape(cat',[],1)));

% now transmit using bpsk
bpsk = 1-2*double(x);

% noise standard deviation
sigman = sqrt(0.5*10.^(-esno_dB/10));

noisy = bpsk + sigman*randn(size(bpsk));

if strcmpi(channel_output,'soft')
    Lc = 4*10^(esno_dB/10);
    LLR = Lc * noisy;    
elseif strcmpi(channel_output,'hard')
    delta = qfunc(1/sigman);
    Lc = log((1-delta)/delta);
    LLR = Lc * sign(noisy);
else
    error('uncorecognized channel output parameter');
end


figure(1);
subplot(2,5,1);
imagesc(reshape(bpsk, 128,[])');
colormap(gray);
set(gca,'xtick',[]); set(gca,'ytick',[]);
axis equal;
xlim([1,128]);
ylim([1,numel(bpsk)/128])
title('Encoded');

subplot(2,5,2);
imagesc(reshape(LLR, 128,[])');
colormap(gray);
set(gca,'xtick',[]); set(gca,'ytick',[]);
axis equal;
xlim([1,128]);
ylim([1,numel(bpsk)/128])
title('Channel output');

iti = 3;
for iter = [1 2 3 4 5 10 15 20]
    subplot(2,5,iti);
    % generate decoder object, built-in MATLAB decoder
    %hdec = comm.LDPCDecoder('ParityCheckMatrix',H,'OutputValue','Whole codeword','MaximumIterationCount',iter);
    
    % use own decoder
    decoded = decode_LDPC(LLR(:)', H, iter);
    
   
    imagesc(reshape(1-decoded, 128,[])');
    colormap(gray);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    axis equal;
    xlim([1,128]);
    ylim([1,numel(decoded)/128])
    iti = iti+1;
    title(sprintf('after %d iter.',iter));
end
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
        if all(L_total > 0)
            break;
        end       
    end
    
    % binary decision    
    xh = L_total < 0;
end
