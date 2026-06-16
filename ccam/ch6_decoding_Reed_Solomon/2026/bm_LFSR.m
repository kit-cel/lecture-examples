function connectionPoly = bm_LFSR(sequence,fieldSize)
%BM_LFSR Given sequence, compute an LSFR that explains sequence using the
% Berlekamp-Massey algorithm
arguments (Input)
    sequence (1,:) {mustBeA(sequence,'gf')} = gf([0,8,3,11],4)
    fieldSize = 4
end

arguments (Output)
    connectionPoly (1,:) {mustBeA(connectionPoly,'gf')}
end
    
    S       = sequence;
    zeroPad = @(x,size) zeroPadGF(x,size,fieldSize); % Useful function

    % Init
    dm = gf(1,fieldSize);  
    c = gf(1,fieldSize);
    p = c;
    L = 0;
    l = 1;

    n = length(sequence);
    for kk = 1:n
      d = S(kk-L:kk)*flip(c(1:L+1))';
      if d==0
          l = l+1;
          continue
      else
          shifted = [gf(zeros(1,l),fieldSize),p];
          if 2*L>=kk
              % No size increase
              l = l+1;
              shifted = zeroPad(shifted,length(c));
              c = c - d/dm*shifted;
          else
              % Size increase
              p = c;
              c = zeroPad(c,kk-L+1);
              l = 1;
              L = length(c) - 1;  
              shifted = zeroPad(shifted,length(c));
              c = c - d/dm*shifted;
              dm = d;
          end
          
      end
    end

    connectionPoly = c;
end

function out = zeroPadGF(x,size,fieldSize)
    out = gf(padarray(x.x,[0,size-length(x)],'post'),fieldSize);
end