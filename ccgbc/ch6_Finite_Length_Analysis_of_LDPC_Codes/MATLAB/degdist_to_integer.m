function [new_dv, num_dv, new_dc, num_dc] = degdist_to_integer(dv, a_dv, dc, a_dc, N)
% function [new_dv, num_dv, new_dc, num_dc] = degdist_to_integer(dv, a_dv, dc, a_dc, N)    
%
    rate = get_rate_LDPC(dv, a_dv, dc, a_dc);
    M = round(N*(1-rate));
    
    
    degprob = optimproblem;
    
    
    num_dc_var_slack = optimvar('num_dc_var_slack', max(dc),'LowerBound',0);     
    num_dv_var_slack = optimvar('num_dv_var_slack', max(dv),'LowerBound',0);     
    
    cost = sum(num_dv_var_slack) + sum(num_dc_var_slack);
    degprob.Objective = cost;
    
    num_dc_var = optimvar('num_dc_var', max(dc), 'Type','integer','LowerBound',0,'UpperBound',M);     
    num_dv_var = optimvar('num_dv_var', max(dv), 'Type','integer','LowerBound',0,'UpperBound',N);     
    
    a_dv_temp = zeros(1,max(dv));
    a_dv_temp(dv) = a_dv;
    dv_temp = 1:max(dv);    
    
    a_dc_temp = zeros(1,max(dc));
    a_dc_temp(dc) = a_dc;
    dc_temp = 1:max(dc);
        
    degprob.Constraints.sumvc = sum(num_dv_var) == N;
    degprob.Constraints.sumdc = sum(num_dc_var) == M;
    degprob.Constraints.edges = dv_temp*num_dv_var - dc_temp*num_dc_var == 0;
    degprob.Constraints.nodv2inc = num_dv_var(2) <= a_dv_temp(2)*N;
    degprob.Constraints.nodv1 = num_dv_var(1) == 0;
    degprob.Constraints.nodc1 = num_dc_var(1) == 0;
    degprob.Constraints.nodc2 = num_dc_var(2) == 0;
    degprob.Constraints.slackv1 = num_dv_var/N - a_dv_temp(:) <= num_dv_var_slack;
    degprob.Constraints.slackv2 = num_dv_var/N - a_dv_temp(:) >= -num_dv_var_slack;
    degprob.Constraints.slackc1 = num_dc_var/M - a_dc_temp(:) <= num_dc_var_slack;
    degprob.Constraints.slackc2 = num_dc_var/M - a_dc_temp(:) >= -num_dc_var_slack;

    [sol,~] = solve(degprob);
 
     if isempty(sol.num_dv_var)
         error('Could not solve optimization problem for finding integer degree distribution');
     end
         
    num_dv = round(double(sol.num_dv_var));
    new_dv = dv_temp;
    
    new_dv = new_dv(num_dv ~= 0);
    num_dv = num_dv(num_dv ~= 0);
    
    num_dc = round(double(sol.num_dc_var));
    new_dc = dc_temp;

    new_dc = new_dc(num_dc ~= 0);
    num_dc = num_dc(num_dc ~= 0);
end

function r = get_rate_LDPC(dv, a_dv, dc, a_dc)  % degree distribution from node perspective
% function r = get_rate_LDPC(dv, a_dv, dc, a_dc)  % degree distribution from node perspective
    vn_poly = node_to_edge(dv, a_dv);
    cn_poly = node_to_edge(dc, a_dc);

    r = 1 - polyval(polyint(cn_poly),1) ./ polyval(polyint(vn_poly),1);
end

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