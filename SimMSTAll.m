%% Simulate MST


%%
Ntrial = 192;
Nsubj = 200;
ModelList = {'SDT1', 'SDT2', 'SDT2', 'HT', 'HT', 'HT', 'HT', 'HT'};
Nmodel = length(ModelList);
param{1} = [1.98 0.99 0.91 1.17]; % SDT-1; d1 d2 s1 s2
param{2} = [2.41 0.76 0.96 0.79 0]; % SDT-2; d1 d2 s1 s2 0
param{3} = [2.41 1.30 0.96 1.33 2.87]; % SDT-2-N; d1 d2 s1 s2 d3
% param{4} = [0.93 0.31 0.38 0.47]; % SDT-1-U; d1 d2 s1 s2
% param{5} = [1.39 0.53 0.59 0.53 0]; % SDT-2-U; d1 d2 s1 s2 0
% param{6} = [1.38 0.71 0.58 0.73 1.24]; % SDT-2-U-n; d1 d2 s1 s2 d3
param{4} = [0.7987    0.0491    0.2253  0.4778 0 0.7987]; % HT-2; pmem, gold, gsim, sigmal, pnew, psim
param{5} = [0.7380    0.2747    0.2077  1 0  0.2979]; % HT-1;
param{6} = [0.7202    0.1563    0.4701    0.4654    0.6272 0.7202]; % HT-2-N;
param{7} = [0.6633    0.4186    0.4399    1 0.7584    0.0710]; % HT-1-N;
param{8} = [0.6404    0.4682    0.2732    1 0.4720    0.4720]; % HT-1-LDI;

% for i = 1:Nmodel
%     param{i} = median(FitResults_MST{i}.Param);
% end

%% Simulation
Data_MST_Sim = cell(Nsubj, Nmodel);
for subj = 1:Nsubj
    for m = 1:Nmodel
        [haha.Nresp,~] = SimMST(ModelList{m}, param{m}, round(Ntrial/3)*ones(3,1));
        haha.Ntrial = Ntrial;
        Data_MST_Sim{subj, m} = haha;
    end
end

%% Fit
% Define models
FitModelSpace={'EVSD3','MST_SD_S2','MST_SD_S2','MST_HT','MST_HT','MST_HT','MST_HT','MST_HT'}; % define model space
Nfit_models=length(FitModelSpace);
FitVariantSpace={{''},{''},{'Newness Separation'},{'Discrimination Threshold'},{'Similar Threshold'},...
    {'Discrimination Threshold','Newness Threshold'},{'Newness Threshold','Similar Threshold'},{'LDI'}}; % define variant space

% Fit models
FitResults_MST_Sim = cell(Nfit_models, Nfit_models);
for model_real=1:Nfit_models
    Config_MA.Data = Data_MST_Sim(:, model_real);
    for model_fit=1:Nfit_models
        if isempty(cell2mat(FitVariantSpace{model_fit}))
            Variants_Display='None';
        else
            Variants_Now=FitVariantSpace{model_fit};
            Variants_Display=cell(1,length(Variants_Now));
            for j=1:length(Variants_Now)
                Variants_Display{j}=[Variants_Now{j},' '];
            end
            Variants_Display=cell2mat(Variants_Display);
        end
        fprintf('\nModel: %s, Variants: %s\n',FitModelSpace{model_fit},Variants_Display) % Progress
        % Specify model
        Config_MA.Model.Model=FitModelSpace{model_fit};
        % Specify variants
        Config_MA.Model.Variants=FitVariantSpace{model_fit};
        % Specify the criteria of interest
        Config_MA.Criteria={'LLH','BIC','AIC','AICc'}; % Calculate criteria of interest
        Config_MA.FitOptions.Algorithm='fmincon: interior-point'; % Change the optimization algorithm (Default: 'DE-MCMC')
        % Configuration
        MA=Configuration_BMW(Config_MA);
        % Estimation
        MA=ModelFit_BMW(MA);
        FitResults_MST_Sim{model_real, model_fit} = MA;
    end
end

% Parameter Recovery
for i = 1:Nfit_models
    fit_now =FitResults_MST_Sim{i,i};
    if i == 2
        fit_now.Param = [fit_now.Param, nan(Nsubj,1)];
    elseif i == 4
        fit_now.Param = [fit_now.Param, nan(Nsubj,2)];
    elseif i == 5
        hahaParam(:,1:3) = fit_now.Param(:,1:3);
        hahaParam(:,4:6) = [nan(Nsubj,2), fit_now.Param(:,4)];
        fit_now.Param = hahaParam;
    elseif i == 6
        fit_now.Param = [fit_now.Param, nan(Nsubj,1)];
    elseif i == 7
        hahaParam(:,1:3) = fit_now.Param(:,1:3);
        hahaParam(:,4:6) = [nan(Nsubj,1), fit_now.Param(:,4:5)];
        fit_now.Param = hahaParam;
    elseif i == 8
        hahaParam(:,1:3) = fit_now.Param(:,1:3);
        hahaParam(:,4:6) = [nan(Nsubj,2), fit_now.Param(:,4)];
        fit_now.Param = hahaParam;
    end
    for j = 1:length(param{i})
        if ~any(isnan(fit_now.Param(:,j)))
        subplot(Nfit_models,6, j + 6*(i-1))
        param_fit = fit_now.Param(:,j);
        histogram(param_fit)
        yticklabels({''})
        hold on
        plot([param{i}(j), param{i}(j)], [0 60],'r--')
        end
    end
end

% Model Recovery
BICM1 = nan(Nfit_models, Nfit_models);
BICM2 = nan(Nfit_models, Nfit_models);
AICM1 = nan(Nfit_models, Nfit_models);
AICcM1 = nan(Nfit_models, Nfit_models);
for i = 1:Nfit_models
    MC = ModelComparison_BMW(FitResults_MST_Sim(i, :));
    BICM1(i,:) = MC{1}.BIC.pBestModel;
    BICM2(i,:) = MC{1}.BIC.ModelFreq;
    AICM1(i,:) = MC{1}.AIC.pBestModel;
    AICcM1(i,:) = MC{1}.AICc.pBestModel;
end

hex_end = {'#003f5c',...
'#ff8d26'};
rgb_end = hex2rgb(hex_end); 
Ncolor = 10;
% linear interpolation
R = linspace(rgb_end(1,1), rgb_end(2,1), Ncolor); 
G = linspace(rgb_end(1,2), rgb_end(2,2), Ncolor); 
B = linspace(rgb_end(1,3), rgb_end(2,3), Ncolor); 
rgb_all = [R' G' B'];
% colormap(rgb_all);  % create colormap

MName = {'SDT-1', 'SDT-2', 'SDT-2-n', 'HT-2', 'HT-1', 'HT-2-N', 'HT-1-N', 'HT-1-LDI'};
figure(1)
colormap(rgb_all);  % create colormap
imagesc(BICM1)
xticklabels(MName)
yticklabels(MName)
colorbar
caxis([0 1]);
axis square
figure(2)
colormap(rgb_all);  % create colormap
imagesc(BICM2)
xticklabels(MName)
yticklabels(MName)
colorbar
caxis([0 1]);
axis square
figure(3)
colormap(rgb_all);  % create colormap
imagesc(AICM1)
xticklabels(MName)
yticklabels(MName)
colorbar
caxis([0 1]);
axis square
figure(4)
colormap(rgb_all);  % create colormap
imagesc(AICcM1)
xticklabels(MName)
yticklabels(MName)
