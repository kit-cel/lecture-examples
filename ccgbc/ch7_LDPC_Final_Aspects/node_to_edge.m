function poly = node_to_edge(dv, a)
% poly = node_to_edge(dv, a)
%
% converts degree distribution from node perspective to distribution from
% edge perspective, as required in density evolution for example
max_dv = max(dv);
L_poly = zeros(max_dv+1,1);
L_poly(max_dv-dv+1) = a;

poly = polyder(L_poly);
poly = poly / polyval(poly,1);

end