%% Simulate MST data
%
%

%%

function [d, p] = SimMST(model, param, Ntrial)

p = zeros(3,3);

switch model
    case 'SDT1'
        dp = [0, param(1), param(1) + param(2)];
        s = [param(3), param(3) + param(4)];
        r0 = 0;
        for t = 1:3 % loop trial type
            % sample evidence
            p(t, 1) = normcdf(s(1), dp(t), 1);
            p(t, 3) = 1 - normcdf(s(2), dp(t), 1);
            p(t, 2) = 1 - p(t, 1) - p(t, 3);
        end
        p(1,1) = r0 + (1-r0)*p(1,1);
        p(1,2) = (1-r0)*p(1,2);
        p(1,3) = (1-r0)*p(1,3);
        
        % simulate multinomial dist
        d = mnrnd(Ntrial, p);
        
    case 'SDT2'
        dp1=param(1); % recognition d-prime
        dp2 = param(2); % discrimination d-prime
        s1 = param(3); % recognition criterion
        s2 = param(4); % discrimination criterion
        dp3 = param(5); % novelty detection d-prime
        r0 = 0;
        
        p(1, 1) = r0 + (1 - r0)*normcdf(s1, 0, 1); % p(new | new)
        p(2, 1) = normcdf(s1, dp1, 1);  % p(new | similar)
        p(3, 1) = p(2, 1); % p(new | old)
        
        p(1, 2) = (1 - r0)*(1 - normcdf(s1, 0, 1)) .* (1 - normcdf(s2, dp3, 1));
        p(2, 2) = (1 - normcdf(s1, dp1, 1)) .* (1 - normcdf(s2, dp2, 1));
        p(3, 2) = (1 - normcdf(s1, dp1, 1)) .* (1 - normcdf(s2, 0, 1));
        
        p(1, 3) = (1 - r0)*(1 - normcdf(s1, 0, 1)) .* normcdf(s2, dp3, 1);
        p(2, 3) = (1 - normcdf(s1, dp1, 1)) .* normcdf(s2, dp2, 1);
        p(3, 3) = (1 - normcdf(s1, dp1, 1)) .* normcdf(s2, 0, 1);
        
        % simulate multinomial dist
        d = mnrnd(Ntrial, p);
        
    case 'SDT1-U'
        dp = [0, param(1), param(1) + param(2)];
        s = [param(3), param(3) + param(4)];
        r0 = 0;
        for t = 1:3 % loop trial type
            p(t, 1) = unifcdf(s(1), dp(t) - 1, dp(t) + 1);
            p(t, 3) = 1 - unifcdf(s(2), dp(t) - 1, dp(t) + 1);
            p(t, 2) = 1 - p(t, 1) - p(t, 3);
        end
        p(1,1) = r0 + (1-r0)*p(1,1);
        p(1,2) = (1-r0)*p(1,2);
        p(1,3) = (1-r0)*p(1,3);
        
        % simulate multinomial dist
        d = mnrnd(Ntrial, p);
        
    case 'SDT2-U'
        dp1=param(1); % recognition d-prime
        dp2 = param(2); % discrimination d-prime
        s1 = param(3); % recognition criterion
        s2 = param(4); % discrimination criterion
        dp3 = param(5); % novelty detection d-prime
        r0 = 0;
        
        p(1, 1) = r0 + (1 - r0)*unifcdf(s1, 0, 1); % p(new | new)
        p(2, 1) = unifcdf(s1, dp1, 1);  % p(new | similar)
        p(3, 1) = p(2, 1); % p(new | old)
        
        p(1, 2) = (1 - r0)*(1 - unifcdf(s1, -1, 1)) .* (1 - unifcdf(s2, dp3-1, dp3+1));
        p(2, 2) = (1 - unifcdf(s1, dp1-1, dp1+1)) .* (1 - unifcdf(s2, dp2-1, dp2+1));
        p(3, 2) = (1 - unifcdf(s1, dp1-1, dp1+1)) .* (1 - unifcdf(s2, -1, 1));
        
        p(1, 3) = (1 - r0)*(1 - unifcdf(s1, -1, 1)) .* unifcdf(s2, dp3-1, dp3+1);
        p(2, 3) = (1 - unifcdf(s1, dp1-1, dp1+1)) .* unifcdf(s2, dp2-1, dp2+1);
        p(3, 3) = (1 - unifcdf(s1, dp1-1, dp1+1)) .* unifcdf(s2, -1, 1);
        
        % simulate multinomial dist
        d = mnrnd(Ntrial, p);
        
    case 'HT'
        % Specify parameters
        r=param(1); % p(remember)
        g_old = param(2); % p(guess old)
        g_sim = param(3); % p(guess similar)
        g_new = 1 - g_old - g_sim;
        sigma_l = param(4);
        r0 = param(5);
        r1 = param(6);
        
        % new
        p(1, 1) = r0 + (1 - r0)*g_new; % p(new | new)
        p(2, 1) = (1 - r1)*g_new;  % p(new | similar)
        p(3, 1) = (1 - r)*g_new; % p(new | old)
        % similar
        p(1, 2) = (1- r0)*g_sim;
        p(2, 2) = r1*sigma_l + (1 - r1)*g_sim;
        p(3, 2) = (1 - r)*g_sim;
        % old
        p(1, 3) = (1 - r0)*g_old;
        p(2, 3) = r1*(1 - sigma_l) + (1 - r1)*g_old;
        p(3, 3) = r + (1 - r)*g_old;
        
        % simulate multinomial dist
        d = mnrnd(Ntrial, p);
        
end

end