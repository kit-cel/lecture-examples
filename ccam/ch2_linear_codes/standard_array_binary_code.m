% Compute the standard array of a binary code

G = [1 0 1 1 0; 
     0 1 0 1 1];
 
 k = size(G,1);
 n = size(G,2);
 
 all_patterns = de2bi(0:2^n-1,n,'left-msb');
 
 % sort patterns that lowest weights start first
 [~,idx] = sort(sum(all_patterns,2));
 all_patterns = all_patterns(idx,:);
 
 used = false(2^n,1);
 
 % first row, all the codewords
 inputs = de2bi(0:2^k-1,k,'left-msb');
 for ki = 1:size(inputs,1)
     % generate codeword    
     x = mod(inputs(ki,:)*G,2);
     
     % mark pattern as used
     idx = find(ismember(all_patterns, x, 'rows')==true, 1, 'first');
     used(idx) = true;
     
     fprintf('%s',strrep(num2str(x),' ',''));
     if ki == 1
         fprintf(' | ');
     elseif ki < size(inputs,1)
         fprintf('   ');
     else
         fprintf('\n');
     end
 end
 fprintf('-----------------------------\n');
 
 
 % complete next rows
 while any(used == false)  % as long as we still have unused patterns
     % find first unused pattern (recall that they are sorted according to
     % their weight)
     free_pattern_idx = find(used==false,1,'first');
     
     used(free_pattern_idx) = true;
     fprintf('%s | ',strrep(num2str(all_patterns(free_pattern_idx,:)),' ',''));
     
     % add error pattern to every codeword
     for ki=2:size(inputs,1)
         x = mod(inputs(ki,:)*G + all_patterns(free_pattern_idx,:), 2);
     
         % mark pattern as used
         idx = find(ismember(all_patterns, x, 'rows')==true, 1, 'first');
         used(idx) = true; 
        
         fprintf('%s',strrep(num2str(x),' ',''));
         if ki < size(inputs,1)
             fprintf('   ');
         else
            fprintf('\n');
         end
     end
 end
     