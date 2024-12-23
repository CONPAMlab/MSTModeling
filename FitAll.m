
dataname1 = 'Data_WM_LTM_BPSO_All3.mat';
load(dataname1)

dataname2 = 'Data_LeeStark.mat';
load(dataname2)

% merge data
Data_All = [Data_BPSO; Data_LeeStark];

%% BPSO & Else
% PSI_all = zeros(Nsubj, 1);
% CSPSI_all = zeros(Nsubj, 1);
% LDI_all = zeros(Nsubj, 1);
% CRS_all = zeros(Nsubj, 1);
% PCS_all = zeros(Nsubj, 1);
% cSD_WM = zeros(Nsubj, 1);
% cSD_LTM = zeros(Nsubj, 1);
% for subj_id = 1:Nsubj
%     pLO(subj_id) = Data_BPSO{subj_id}.LO/sum(Data_BPSO{subj_id}.Nresp(2,:));
%     pLN(subj_id) = Data_BPSO{subj_id}.LN/sum(Data_BPSO{subj_id}.Nresp(2,:));
%     pLS(subj_id) = Data_BPSO{subj_id}.LS/sum(Data_BPSO{subj_id}.Nresp(2,:));
%     pTN(subj_id) = Data_BPSO{subj_id}.TN/sum(Data_BPSO{subj_id}.Nresp(3,:));
%     PSI_all(subj_id) = Data_BPSO{subj_id}.PSI;
%     CRS_all(subj_id) = Data_BPSO{subj_id}.CRS;
%     PCS_all(subj_id) = Data_BPSO{subj_id}.PCS;
%     CSPSI_all(subj_id) = Data_BPSO{subj_id}.CSPSI;
%     LDI_all(subj_id) = Data_BPSO{subj_id}.LDI;
%     error_WM = CircDist_BMW('Diff', Data_WM{subj_id}.response, Data_WM{subj_id}.sample, 180);
%     cSD_WM(subj_id) = circ_std(error_WM/180*pi)/pi*180;
% %     if valid_LTM_id(subj_id) == 1
% %         error_LTM = CircDist_BMW('Diff', Data_WM{subj_id}.response, Data_WM{subj_id}.sample, 360);
% %         cSD_LTM(subj_id) = circ_std(error_LTM/180*pi)/pi*180;
% %     end
% end

for subj_id = 1:Nsubj
    LDI_data(subj_id) = Data_All{subj_id}.Nresp(2,2)./(sum(Data_All{subj_id}.Nresp(2,:))) - ...
        Data_All{subj_id}.Nresp(1,2)./(sum(Data_All{subj_id}.Nresp(1,:)));
end

Nsubj = length(Data_All);


% Define models
% FitModelSpace={'EVSD3','MST_SD_S2','MST_SD_S2','EVSD3_LDI','MST_SD_S2_LDI','MST_SD_S2_LDI','MST_HT','MST_HT','MST_HT','MST_HT'}; % define model space
% Nfit_models=length(FitModelSpace);
% FitVariantSpace={{''},{''},{'Newness Separation'},{''},{''},{'Newness Separation'},{'Discrimination Threshold'},{'Similar Threshold'},...
%     {'Discrimination Threshold','Newness Threshold'},{'Newness Threshold','Similar Threshold'}}; % define variant space

FitModelSpace={'EVSD3','MST_SD_S2','MST_SD_S2','MST_HT','MST_HT','MST_HT','MST_HT','MST_HT'}; % define model space
Nfit_models=length(FitModelSpace);
FitVariantSpace={{''},{''},{'Newness Separation'},{'Discrimination Threshold'},{'Similar Threshold'},...
    {'Discrimination Threshold','Newness Threshold'},{'Newness Threshold','Similar Threshold'},{'LDI'}}; % define variant space


% Fit models
Config_MA.Data=Data_All;
FitResults_MST = cell(Nfit_models, 1);
for model_id=1:Nfit_models
    if isempty(cell2mat(FitVariantSpace{model_id}))
        Variants_Display='None';
    else
        Variants_Now=FitVariantSpace{model_id};
        Variants_Display=cell(1,length(Variants_Now));
        for j=1:length(Variants_Now)
            Variants_Display{j}=[Variants_Now{j},' '];
        end
        Variants_Display=cell2mat(Variants_Display);
    end
    fprintf('\nModel: %s, Variants: %s\n',FitModelSpace{model_id},Variants_Display) % Progress
    % Specify model
    Config_MA.Model.Model=FitModelSpace{model_id};
    % Specify variants
    Config_MA.Model.Variants=FitVariantSpace{model_id};
    % Specify the criteria of interest
    Config_MA.Criteria={'LLH','BIC','AIC','AICc'}; % Calculate criteria of interest
    Config_MA.FitOptions.Algorithm='fmincon: interior-point'; % Change the optimization algorithm (Default: 'DE-MCMC')
    % Configuration
    MA=Configuration_BMW(Config_MA);
    % Estimation
    MA=ModelFit_BMW(MA);
    FitResults_MST{model_id} = MA;
end

dataResp = nan(3,3,153);
for i = 1:153
    dataResp( :,:,i) = Data_All{i}.Nresp ./ repmat(sum(Data_All{i}.Nresp,2), [1,3]);
end

PPdata = [reshape(mean(dataResp(:,:,1)),[1,3]); reshape(mean(dataResp(:,:,2)),[1,3]); reshape(mean(dataResp(:,:,3)),[1,3])];
PPerr = [reshape(std(dataResp(:,:,1)),[1,3])/sqrt(153); reshape(std(dataResp(:,:,2)),[1,3])/sqrt(153);...
    reshape(std(dataResp(:,:,3)),[1,3])/sqrt(153)];
% PPx=categorical({'\it CL','\it DM'});
figure
hBar = bar(PPdata);
hBar(1).XData
hold on
for k1 = 1:size(PPdata,2)
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');   % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
end
hold on
errorbar(ctr, ydt, PPerr.', '.k')
PPx=categorical({'Foil','Lure', 'Target'});
set(gca, 'XTickLabel', PPx)
box off
set(gca, 'LineWidth',2)

% Model Comparison
% AllIC = zeros(Nsubj, Nfit_models);
% for m = 1:Nfit_models
%     AllIC(:, m) = FitResults_MST{m}.BIC;
% end
% imagesc(AllIC)
MC = ModelComparison_BMW(FitResults_MST(1:end)');


config.Criteria = 'BIC';
Visualization_BMW(MC{1}, 'ModelSelectionColormap',config)

config.Criteria = 'AIC';
Visualization_BMW(MC{1}, 'ModelSelectionColormap',config)

MName = {'SDT-1', 'SDT-2', 'SDT-2-n', 'HT-2', 'HT-1', 'HT-2-N', 'HT-1-N', 'HT-1-LDI'};
% subplot(1,2,1)
% plot(1:7, MC{1}.BIC.ModelFreq,'.-', 'LineWidth',2)
% xlim([0.5 7.5])
% xticks(1:7)
% xticklabels(MName)
% xtickangle(45)
% box off
% subplot(1,2,2)
% plot(1:7, MC{1}.BIC.EP,'.-','LineWidth',2)
% xlim([0.5 7.5])
% xticks(1:7)
% xticklabels(MName)
% xtickangle(45)
% box off
pBestM_BIC = zeros(Nfit_models,1);
pBestM_AIC = zeros(Nfit_models,1);
for i = 1:Nfit_models
    pBestM_BIC(i) = sum(MC{1}.BIC.BestModel == i)/153;
    pBestM_AIC(i) = sum(MC{1}.AIC.BestModel == i)/153;
end

plot(1:Nfit_models, pBestM_AIC,'*-', 'LineWidth',2)
hold on
plot(1:Nfit_models, pBestM_BIC,'*--', 'LineWidth',2)
plot(MC{1}.BIC.ModelFreq,'*:', 'LineWidth',2)
plot(MC{1}.BIC.EP,'*-', 'LineWidth',2)
xlim([0.5 Nfit_models + 0.5])
xticks(1:Nfit_models)
xticklabels(MName)
xtickangle(45)
box off
axis square
set(gca, 'LineWidth',2)

subplot(1,3,1)
fam = [1 1 1 2 2 2 2 2]; % sdt vs. threshold
% fam = [1 2 2 2 1 2 1 0]; % single vs. dual-process
% fam = [0 0 0 1 1 2 2 0]; % Novelty detection
MC = ModelComparison_BMW(FitResults_MST(1:Nfit_models)',fam);
pBestFam_BIC = zeros(max(fam),1);
pBestFam_AIC = zeros(max(fam),1);
for i = 1:max(fam)
    pBestFam_BIC(i) = sum(MC{1}.BIC.BestFamily == i)/153;
    pBestFam_AIC(i) = sum(MC{1}.AIC.BestFamily == i)/153;
end

plot(pBestFam_AIC,'*-', 'LineWidth',2)
hold on
plot(pBestFam_BIC,'*--', 'LineWidth',2)
plot(MC{1}.BIC.FamilyFreq,'*:', 'LineWidth',2)
plot(MC{1}.BIC.EP_Family,'*-', 'LineWidth',2)
xlim([0.7 max(fam)+.3])
xticks(1:max(fam))
xticklabels({'SDT', 'Threshold'})
xtickangle(45)
box off
ylim([0 1])
set(gca, 'LineWidth',2)

subplot(1,3,2)
% fam = [1 1 1 2 2 2 2 0]; % sdt vs. threshold
fam = [1 2 2 2 1 2 1 1]; % single vs. dual-process
% fam = [0 0 0 1 1 2 2 0]; % Novelty detection
MC = ModelComparison_BMW(FitResults_MST(1:Nfit_models)',fam);
pBestFam_BIC = zeros(max(fam),1);
pBestFam_AIC = zeros(max(fam),1);
for i = 1:max(fam)
    pBestFam_BIC(i) = sum(MC{1}.BIC.BestFamily == i)/153;
    pBestFam_AIC(i) = sum(MC{1}.AIC.BestFamily == i)/153;
end

plot(pBestFam_AIC,'*-', 'LineWidth',2)
hold on
plot(pBestFam_BIC,'*--', 'LineWidth',2)
plot(MC{1}.BIC.FamilyFreq,'*:', 'LineWidth',2)
plot(MC{1}.BIC.EP_Family,'*-', 'LineWidth',2)
xlim([0.7 2.3])
xticks(1:2)
xticklabels({'1-D', '2-D'})
xlabel('Mnemonic Evidence')
xtickangle(45)
box off
ylim([0 1])
set(gca, 'LineWidth',2)

subplot(1,3,3)
% fam = [1 1 1 2 2 2 2 0]; % sdt vs. threshold
% fam = [1 2 2 2 1 2 1 0]; % single vs. dual-process
fam = [0 0 0 1 1 2 2 2]; % Novelty detection
MC = ModelComparison_BMW(FitResults_MST(1:Nfit_models)',fam);
pBestFam_BIC = zeros(max(fam),1);
pBestFam_AIC = zeros(max(fam),1);
for i = 1:max(fam)
    pBestFam_BIC(i) = sum(MC{1}.BIC.BestFamily == i)/153;
    pBestFam_AIC(i) = sum(MC{1}.AIC.BestFamily == i)/153;
end

plot(pBestFam_AIC,'*-', 'LineWidth',2)
hold on
plot(pBestFam_BIC,'*--', 'LineWidth',2)
plot(MC{1}.BIC.FamilyFreq,'*:', 'LineWidth',2)
plot(MC{1}.BIC.EP_Family,'*-', 'LineWidth',2)
xlim([0.7 2.3])
xticks(1:2)
xticklabels({'No', 'Yes'})
xlabel('Novelty Detection')
xtickangle(45)
box off
ylim([0 1])
set(gca, 'LineWidth',2)

%% Correlation
plot(FitResults_MST{1}.Param(:,2)./(FitResults_MST{1}.Param(:,1) + FitResults_MST{1}.Param(:,2)), LDI_data','.');
xlim([0 .8])

plot(FitResults_MST{5}.Param(:,4), LDI_data','.');
xlim([0 .8])

plot(FitResults_MST{8}.Param(:,4), LDI_data','.');
xlim([0 .8])


%% Check prediction
mean_dataM = mean(dataResp, 3);
sd_dataM = std(dataResp, [], 3);
Ntrial = 192;
% Simulate SDT-1
p_now1 = zeros(3,3,Nsubj);
for i = 1:Nsubj
    [~, p_now1(:,:,i)] = SimMST('SDT1', FitResults_MST{1}.Param(i,:), round(Ntrial/3)*ones(3,1));
end
pred_SDT1 = mean(p_now1, 3);
std_SDT1 = std(p_now1,[],3);
% Simulate HT-1-LDI
p_now2 = zeros(3,3,Nsubj);
for i = 1:Nsubj
    [~, p_now2(:,:,i)] = SimMST('HT', [FitResults_MST{8}.Param(i,1:3), 1, FitResults_MST{8}.Param(i,4), FitResults_MST{8}.Param(i,4)], round(Ntrial/3)*ones(3,1));
end
pred_HT = mean(p_now2, 3);
std_HT = std(p_now2,[],3);
% % Simulate HT-2-N
% p_now3 = zeros(3,3,Nsubj);
% for i = 1:Nsubj
%     [~, p_now3(:,:,i)] = SimMST('HT', [FitResults_MST{6}.Param(i,:), FitResults_MST{6}.Param(i,1)], round(Ntrial/3)*ones(3,1));
% end
% pred_HT2N = mean(p_now3, 3);

mean_y = [mean_dataM(1:6); pred_SDT1(1:6); pred_HT(1:6)];
var_y = [sd_dataM(1:6); std_SDT1(1:6); std_HT(1:6)]./sqrt(Nsubj);


PPdata = mean_y';
PPerr = var_y';
% PPx=categorical({'\it CL','\it DM'});
figure
hBar = bar(PPdata);
hBar(1).XData
hold on
for k1 = 1:size(PPdata,2)
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');   % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
end
hold on
errorbar(ctr, ydt, PPerr.', '.k')
PPx=categorical({'p(Old | Target)','p(Similar | Target)', 'p(New | Target)',...
    'p(Old | Lure)','p(Similar | Lure)', 'p(New | Lure)'});
set(gca, 'XTickLabel', PPx)
box off
set(gca, 'LineWidth',2)

