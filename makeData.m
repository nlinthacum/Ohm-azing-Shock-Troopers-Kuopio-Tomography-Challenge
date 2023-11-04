%script to create test data with missing data values
categoryNbr = 6;
sample_idx = 1;

    
load(['/Users/nick/Documents/Dbar_otsu/InputData1/data' num2str(sample_idx) '.mat'])

vincl = true(31,76); %which measurements are used in the inversion - 31 voltage values for each of the 76 current injections
rmind = 1:2*(categoryNbr - 1);
load('/Users/nick/Documents/Dbar_otsu/InputData1/ref.mat')
for ii=1:76
    if any(Injref(rmind,ii) ~= 0)
        vincl(:,ii) = false; %remove all voltage measurements for this current injection
    end
    vincl(rmind,ii) = false; %remove all voltage measurements collected using these electrodes
end
vincl = vincl(:);


nonzeroColumns = any(Inj(1:2*(categoryNbr - 1), :) ~= 0);

% Delete columns with nonzero values in row 1
Inj(:, nonzeroColumns) = nan;
Inj(1:2*(categoryNbr - 1), :) = nan;

Uel(vincl == 0) = nan; % = Uel(vincl == 1);
save(['InputData2/data' num2str(sample_idx)], 'Inj', 'Uel')


% load('/Users/nick/Documents/Dbar_otsu/InputData/ref.mat')
% 
% nonzeroColumns = any(Injref(1:2*(categoryNbr - 1), :) ~= 0);
% 
% % Delete columns with nonzero values in row 1
% Injref(:, nonzeroColumns) = nan;
% Injref(1:2*(categoryNbr - 1), :) = nan;
% 
% Uelref(vincl == 0) = nan;
% save('InputData3/ref', 'Injref', 'Uelref')
% 
