%% MST: Two-step Signal Detection
%
%

function Output = MST_SD_S2(param, Data, Input)
% Specify parameters
dp1=param(1); % recognition d-prime
dp2 = param(2); % discrimination d-prime
s1 = param(3); % recognition criterion
s2 = param(4); % discrimination criterion
Nparam = 4;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
if any(strcmp(Input.Variants,'Newness Threshold'))
    Nparam=Nparam+1;
    r0 = param(Nparam); % Threshold
else
    r0 = 0;
end
if any(strcmp(Input.Variants,'Newness Separation'))
    Nparam=Nparam+1;
    dp3 = param(Nparam); % Threshold
else
    dp3 = 0;
end

RespMat = Data.Nresp;
% N_total = sum(Data.Nresp, 2);
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

p_now = zeros(3,3);

p_now(1, 1) = r0 + (1 - r0)*normcdf(s1, 0, 1); % p(new | new)
p_now(2, 1) = normcdf(s1, dp1, 1);  % p(new | similar)
p_now(3, 1) = p_now(2, 1); % p(new | old)

p_now(1, 2) = (1 - r0)*(1 - normcdf(s1, 0, 1)) .* (1 - normcdf(s2, dp3, 1));
p_now(2, 2) = (1 - normcdf(s1, dp1, 1)) .* (1 - normcdf(s2, dp2, 1));
p_now(3, 2) = (1 - normcdf(s1, dp1, 1)) .* (1 - normcdf(s2, 0, 1));

p_now(1, 3) = (1 - r0)*(1 - normcdf(s1, 0, 1)) .* normcdf(s2, dp3, 1);
p_now(2, 3) = (1 - normcdf(s1, dp1, 1)) .* normcdf(s2, dp2, 1);
p_now(3, 3) = (1 - normcdf(s1, dp1, 1)) .* normcdf(s2, 0, 1);

% for t = 1:3 % loop trial type
%     p_now(t, 1) = normcdf(s(1), dp(t), 1);
%     p_now(t, 3) = 1 - normcdf(s(2), dp(t), 1);
%     p_now(t, 2) = 1 - p_now(t, 1) - p_now(t, 3);
% end
%  LLH_total = -log(mnpdf(RespMat, p_now));
%  LLH = sum(LLH_total); 
 LLH_total = -log(p_now).*RespMat;
 LLH = sum(sum(LLH_total));   

% Posterior
LP=-log(Prior)+LLH; % likelihood*prior

if LP==Inf || isnan(LP)
    LP=realmax('double'); % Output should be a real value
    LLH=LP;
end

% PPD
p_LH=zeros(sum(sum(RespMat)),1);
N_LH=0;

% for i=1:Nc
%     p_LH(N_LH+1:N_LH+Ntarget+Nlure,1)=[p_hit(i)*ones(Nhit(i),1); (1-p_hit(i))*ones(Ntarget-Nhit(i),1); p_fa(i)*ones(Nfa(i),1); (1-p_fa(i))*ones(Nlure-Nfa(i),1)];
%     N_LH=N_LH+Ntarget+Nlure;
% end

% Decide output
if strcmp(Input.Output,'LP')
    Output=LP;
elseif strcmp(Input.Output,'LLH')
    Output=LLH;
elseif strcmp(Input.Output,'Prior')
    Output=Prior;
elseif strcmp(Input.Output,'LPPD')
    Output=log(p_LH);
elseif strcmp(Input.Output,'All')
    Output.LP=LP;
    Output.LLH=LLH;
    Output.Prior=Prior;
    Output.LPPD=log(p_LH);
end

end

% Define prior
function p=prior(param, Input)
% 
% % Specify parameters
% K=param(1); % Capacity
% % weibull prior for capacity
% p0(1)=wblpdf(K,3.5,3); % given that K is ofter 3~4
% kappa_1=param(2); % Unit resource
% % Gamma prior for unit resource
% p0(2)=gampdf(kappa_1,3,5);
% Nparam=2;
% 
% if ~isfield(Input,'Variants') % No Variants
%     Input.Variants={};
% end
% if any(strcmp(Input.Variants,'ResponseNoise'))
%     kappa_r=param(3); % Response variability
%     % Gamma prior for response noise
%     p0(3)=gampdf(kappa_r,3,5);
% end
% if any(strcmp(Input.Variants,'Bias'))
%     Nparam=Nparam+1;
%     bias=param(Nparam); % Mean bias
%     % Gaussian prior for bias
%     p0(Nparam)=normpdf(bias, 0, 1);
% end
% if any(strcmp(Input.Variants,'BiasF'))
%     Nparam=Nparam+1;
%     biasF=param(Nparam); % Fluctuation of bias
%     % Gaussian prior for the fluctuation of bias
%     p0(Nparam)=normpdf(biasF, 0, 5);
% end
% if any(strcmp(Input.Variants,'PrecF'))
%     Nparam=Nparam+1;
%     precF=param(Nparam); % Fluctuation of precision
%     % Gaussian prior for the fluctuation of precision
%     p0(Nparam)=normpdf(precF, 0, 1);
% end
% if any(strcmp(Input.Variants,'Swap'))
%     Nparam=Nparam+1;
%     s=param(Nparam); % Swap rate
%     % Gaussian prior for the swap rate
%     p0(Nparam)=normpdf(s, 0.5, 1);
% end
% 
% % Construct joint distribution
% % Consider independent parameters here
% % We think it's generally acceptable for prior definition,
% % tho it's usually not the actual case
% p=1;
% for i=1:Nparam
%     p=p*p0(i);
% end
p=1;

end
