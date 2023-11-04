function result = DbarRecon(inputFolder, outputFolder, categoryNbr)
%==========================================================================
% This script runs the Dbar algorithm, calling the necessary functions to
% compute the approx scattering transform texp and solve the Dbar equation.
%
% Authors:      Melody Alsaker, Nick Linthacum, and Siiri Rautio
% Date Modified:            September 2014, October 2023
%  
%==========================================================================
%clear all;
input_dir = inputFolder;%'InputData/';

%find number of samples to reconstruct
S = dir(fullfile(inputFolder,'**','*.mat'));






for sample = 1:numel(S)-1 %loop over all the data samples -1 since ref.mat is included
%The sample number refers to the number in its name

   
    total_gamma = []; 
    
    %=========================== Load external data ===========================
% Vref is the data set with nothing in the tank (reference data) 
% Vmulti is the data set with something in the tank (measured data)
% There is only 1 frame for this data
load([input_dir '/' 'ref.mat'])

sample_name = [input_dir '/' 'data' num2str(sample)];

load(sample_name)
    
%=========================== Modify Data based on difficulty level ===========================
vincl = true(31,76); %which measurements are used in the inversion - 31 voltage values for each of the 76 current injections
rmind = 1:2*(categoryNbr - 1);
for ii=1:76
    if any(Injref(rmind,ii) ~= 0)
        vincl(:,ii) = false; %remove all voltage measurements for this current injection
    end
    vincl(rmind,ii) = false; %remove all voltage measurements collected using these electrodes
end
vincl = vincl(:);

nonzeroColumns = any(Injref(1:2*(categoryNbr - 1), :) ~= 0);

% Delete columns with nonzero values in row 1
Inj(:, nonzeroColumns) = [];
Inj(1:2*(categoryNbr - 1), :) = [];

Uel = Uel(vincl == 1);

% Delete columns with nonzero values in row 1
Injref(:, nonzeroColumns) = [];
Injref(1:2*(categoryNbr - 1), :) = [];

Uelref = Uelref(vincl == 1);


num_itr = floor(size(Inj, 2) / rank(Inj));

starting_idx = [1:rank(Inj):size(Inj,2)];%%starting current injection, the list 
starting_idx = starting_idx(1:num_itr);
%is used to do recons with more data
for starting_i = 1:length(starting_idx) %loop through the starting current injection to use all of the given data


[m,n] = size(Inj);
L = m - 1;%num electrodes 
num_cur_injections = rank(Inj);
cur_start = starting_idx(starting_i);
data_start = starting_idx(starting_i);
data_range = (((data_start-1)*L)+1):(((data_start-1)*L+1) + (L*num_cur_injections)-1);


Vref = Uelref*1000; 

Vref = Vref(data_range); 

Vmulti = Uel(data_range)*1000;
Vmulti = Vmulti';

[total_num_frames,~] = size(Vmulti); %here, we are only usng 1 frame

%========================Set up numerical parameters=======================
imidx = 1; % Defines which frames we want to reconstruct
num_frames = numel(imidx);      % Number of frames we want to reconstruct  

% Grid and computational parameters
s = 4.3; % Defines size of k-grid square. Must be > 2
s_str = ['recon_z_grid_4p3_sample' num2str(sample)];
M = 16;              % Size of k-grid is M x M  
h = 2*s/(M-1);       % k-grid step size 
trunc = s;           % Truncation radius for k-grid
hh = 0.01;% was 0.03 before       % z-grid step size

% Physical parameters
                               % Number of electrodes
dtheta = 2*pi/L;                       % Ang. dist. between electrode ctrs
etheta = pi/2:dtheta:5*pi/2-dtheta;    % Angs. corres. to electrode ctrs
eheight = 0.013; ewidth = 0.0254;      % Electrode height, width (meters)
%eheight = 0.07; ewidth = 0.025;    
Radius = 0.23; %0.15;                         % Radius of circular domain (meters) 
efrac = (ewidth*L)/(2*pi*Radius);      % Frac. of bndry covered by electrodes
eArea2 = efrac * dtheta * ...          % Area of electrode (m^2)
           Radius * eheight;        
skip =   1; % Which skip pattern
numCP = L - gcd(skip+1,L);                  % Num. of lin. indep. C.P.s

% Get Current amplitude CurrAmp (mA)
CurrAmp = .5*(mean(mean(Inj)) + mean(mean(Inj)));%.5*(mean(mean(curr_in_data)) + mean(mean(curr_out_data)));
CurrAmp = CurrAmp / 2;
%======================Set up computational grids==========================

%..........................................................................
% Construct mesh of z-values representing physical domain. We can throw out
% z-values outside the domain; these will not be needed in the computation.
%..........................................................................
xx = -1:hh:1;
N = numel(xx);
[z1,z2] = meshgrid(xx,xx);
z = z1 + 1i*z2;
z = reshape(z,N*N,1);
zidx = find(abs(z) <= 1);           % Indices of z-values in domain
z = z(zidx);                        % Get rid of values outside domain.
conjz = conj(z);
numz = numel(z);

%..........................................................................
% Construct computational grid with M x M elements & complex variable k
% We will need to eliminate the k-values outside the truncation radius and
% within a small radius around k=0, but we will need to keep the original
% k-grid for some computations.
%..........................................................................
x = -s:h:s;        
[K1, K2] = meshgrid(x,x);  
k = K1 + 1i*K2;
numk = M*M;
kidx = find(abs(k)<trunc & abs(k)>0.1); % Indices of k-vals in trunc area
ktrunc = k(kidx);                       % k-vals inside trunc area (vector)
numktrunc = numel(ktrunc);              % Number of k-vals in trunc area
conjktrunc = conj(ktrunc);              % conj of trunc k-vals (vector)
conjk= conj(k);

% The grid for the Green's function beta needs to be larger to accomodate 
% the convolution.
xBig            = [-(s+((M/2):-1:1)*h),x,s+(1:(M/2))*h];
[K1Big, K2Big]  = meshgrid(xBig,xBig);
k_Big           = K1Big + 1i*K2Big;


%======================Define Green's function beta========================

beta = 1 ./(pi * k_Big);
beta(M+1,M+1)=0;   % Set beta(0,0) = 0 to avoid singularity.
beta = h^2*beta;   % Multiply by h^2 for when we compute the convolution

% Take the fast fourier transform of the Green's function. 
% This is an averaging (Andreas's form)
p1 = shiftr(1:(2*M),0,M+1,1);
p  = shiftr(1:(2*M),0,M,1);
fft_beta  = fft2(beta(p,p)+beta(p1,p1)+beta(p1,p)+beta(p,p1))/4;

%fft_beta = fft2(beta);  % This doesn't work!!!

%======================= Construct Current Matrix J =======================

skip_vector = zeros(L,1);
skip_vector(1) = CurrAmp;
skip_vector(2+skip) = -CurrAmp;
% J = gallery('circul', skip_vector)';

J = Inj;
J = J(:, cur_start:(cur_start +num_cur_injections)-1);
numCP = num_cur_injections;

J = J(1:L,1:numCP);
Q = sparse(GS_orthonorm(J));    % Orthonormalize
CoB = Q.'*J;            % Coeffs for change of basis for voltage matrices
J = Q;                
J = J/eArea2;

%================ Construct DN matrix for reference data ==================

% Voltage matrix. Columns are CP and each row is an electrode
Vref = reshape(Vref,L, num_cur_injections);  
 

% Normalize the entries so that the voltages sum to zero in each col.
adj = sum(Vref)/L;
Vref = Vref - adj(ones(L,1),:);%Nick change Vref - adj(ones(L,1),:);

% Scale columns to match L2-normalization of the C.P. matrix J
Vref = Vref / CoB;

refLambda = inv(Vref' * J);   % refLambda is the DN map, size numCP x numCP 

%============= Precompute some values necessary for Dbar eqn ==============

%..........................................................................
% DC encodes information for the coeffs c, d used in scattering transform.
% The coeffs. c, d are given by:
% c = J' * exp(1i*ktrunc(ii)*zee).' and 
% d = J' * exp(1i*conjktrunc(ii)*conjzee).'
% These would be vectors of length numCP whose entries are the coeffs from
% the eigenfxn expansion theorem. Encoding the information into the matrix
% DC avoids loops later when the scattering transform is computed.
%..........................................................................
zee = (cos(etheta) + 1i * sin(etheta)); % Centers of electrodes
conjzee = conj(zee);
DC = zeros(numCP, numCP , numktrunc);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
for ii = 1:numktrunc
    Jexp1 = (exp(1i*ktrunc(ii)*zee) * J);
    Jexp2 = (exp(1i*conjktrunc(ii)*conjzee) * J).';
    DC(:,:,ii) = Jexp2(:,ones(1,numCP)) .* Jexp1(ones(numCP,1),:);
end

%..........................................................................
% EXP encodes the factor exp(-i(kz+conj(kz))) / (4pi*conj(k)) used in Dbar 
% eqn. This will be multiplied by the scattering transform later to form
% the pointwise multiplication operator TR.
%..........................................................................
EXP = zeros(M,M, numz);
for ii = 1:numz
   EXP(:,:,ii) = exp(-1i*(k*z(ii) + conjk*conjz(ii)))./ conjk; 
end
EXP = EXP / (4*pi);

% Construct rhs of eqn DBO*m = 1. We also use rhs for the init guess.
rhs = [ones(numk,1); zeros(numk,1)];


%======Loop through all data sets, reconstruct conductivity for each ======

% gamma is the conductivity we will reconstruct. Outside domain will be NaN
gamma = ones(num_frames,N*N) * NaN;
gamma_best = 1;
jj = 1; %only 1 frame 
    %================= Construct DN matrix for measured data ==============
    V = Vmulti(imidx(jj),:);
    V  = reshape(V, L,num_cur_injections);
    V  = V(:,1:numCP);
    Vorig = V;
    % Normalize the entries so that the voltages sum to zero in each col.
    adj = sum(V)/L;
    V = V - adj(ones(L,1),:);
    
    % Scale columns to match L2-normalization of the C.P. matrix J
    V = V / CoB;
  
    
    Lambda = inv(V' * J );
    
    dLambda = Lambda - refLambda;
    
   
    %==================Compute approx. scattering transform================
    
    texp = zeros(numk,1);
    %..................................................................
    % Below is a fast computation for the sum:
%         c = zeros(numktrunc, numCP);
%         d = zeros(numktrunc, numCP);
%         for ii = 1:numktrunc
%         c(ii,:) = (exp(1i*ktrunc(ii)*zee) * J ).' ;
%         d(ii,:) = (exp(1i*conjktrunc(ii)*conjzee) * J ).';
%             for EM = 1:numCP
%                 for em = 1:numCP
%                     texp(kidx(ii)) = texp(kidx(ii)) + d(ii,EM) * c(ii,em) * dLambda(EM,em);
%                 end
%             end
%         end

    % The coeffs. c, d are vectors of length numCP whose entries are the 
    % coeffs from the eigenfxn expansion thm.
    % These have been encoded into the matrix DC computed earlier.
    %..................................................................
    
    texp(kidx) = squeeze(sum(sum(DC .* repmat(dLambda,[1,1,numktrunc]),1),2));

    texp = texp * Radius * dtheta /(gamma_best) ; % Scaling factors
    
    % This is the pointwise multiplication operator used in the Dbar eqn.
    TR = repmat(reshape(texp,M,M),[1,1,numz]) .* EXP;
    
    %==========================Solve Dbar Equation=========================   
    
    % Loop through all z values in domain
    for ii = 1:numz
        T = TR(:,:,ii);
        [m,~] = bicgstab(@DBO, rhs, 1e-5, 10, [], [], rhs, M, T, fft_beta);
        
        %------------------------------------------------------------------
        % Use this version if you are using the nested function DBop
        % defined within this file. This will improve runtime. 
        % If this option is used, DBop must be un-commmented below,
        % and this entire file must be run as a function rather than
        % as a script.
        %------------------------------------------------------------------
        %[m,~] = bicgstab(@DBop, rhs, 10, 1e-5, 10, [],[],rhs);
        %------------------------------------------------------------------
        
        sqrtgamma = m((numk+M)/2 +1) + 1i * m( (3*numk + M)/2 + 1);
        gamma(jj,zidx(ii)) = sqrtgamma * sqrtgamma;
    end
    

%==========================================================================

% Make each conductivity distribution into a matrix, one for each image.
tmp_gamma = reshape(gamma,num_frames,N,N)*gamma_best;  

 total_gamma = cat(1, total_gamma,real(tmp_gamma));  
end %end looping through all starting points
gamma = nanmean(total_gamma, 1);





%==========================================================================
% Nested function for the operator used in Dbar Eqn. This replaces the
% function call to the m-file DBO.m, thus reducing overhead caused by
% calls to external files. Comment this out if you're using the version
% that calls DBO.m in the solution to the Dbar equation (desirable if you
% are running this file as a script rather than as a function).
%==========================================================================
% function result = DBop(f)
%     f = f(1:numk) + 1i * f(numk+1:2*numk);
% 
%     % Construct conj(matf) .* T with zero padding to accommodate convolution
%     temp_fT = conj(reshape(f,M,M)) .* T; 
%     temp_fT_Big  = zeros(2*M);
%     temp_fT_Big((M/2+1):(3*M/2), (M/2+1):(3*M/2)) = temp_fT;  % zero padding
% 
%     % Compute action of operator on f. Note, the h^2 is already included in
%     % beta, but it could go here instead
%     temp_fT_Big = ifftn(fft_beta.*fftn(temp_fT_Big));
%     tmp = temp_fT_Big((M/2+1):(3*M/2), (M/2+1):(3*M/2)); %Remove zero padding
%     f = f - tmp(:);
% 
%     %Stack real and imaginary parts
%     result = [real(f); imag(f)];   
% end
%==========================================================================


num_frames_to_plot = num_frames;  % Use num_frames to plot all
gam_real = real(gamma);
datamin = min(min(min(gam_real)));
datamax = max(max(max(gam_real)));
datarange = datamax-datamin;
colorbartrunc = 0.0;
datamin = datamin + colorbartrunc * datarange;
datamax = datamax - colorbartrunc * datarange;

%% ========= Use this to display or and/or save individual images ===========
for jj = 1:num_frames_to_plot
    %h = figure('visible', 'off');  % Suppress display to screen
%     save(['Dbar', num2str(trunc), '_M', num2str(M),'_zstep', num2str(sample), '.mat'],'gam_real' );
    h=figure('visible','off');
    % imagesc(xx,xx,squeeze(gam_real(jj,:,:)),[datamin, datamax]);
    imagesc(xx,xx,squeeze(gam_real(jj,:,:)),[datamin, datamax]);
    set(gca, 'Ydir', 'normal');
    title(['Frame number ', num2str(imidx(jj))]);
    colorbar;
    axis([-1 1 -1 1 ]);
    axis square;
    colormap parula;
   
    %savefig(s_str);
    %close all;
    %print(h,'-djpeg', ['frame_', num2str(imidx(jj)), '.jpg']); % Save image
end
% ==========================================================================

%========== Use this to save workspace variables ==========================
%     save(['Dbar', num2str(trunc), '_M', num2str(M),'_zstep', num2str(sample), '.mat'],'gam_real' );
%end

%% ============== Thresholding with Otsu =================================
% Pre-process
result = flip(squeeze(gam_real(jj,:,:)));
result_256 = imresize(result,[256,256]);

% Treshold the image histogram using Otsu's method
[level,x] = Otsu2(result_256,256);

reconstruction = zeros(size(result_256));

ind1 = find(result_256 < x(level(1)));
ind2 = find(result_256 >= x(level(1)) & result_256 <= x(level(2)));
ind3 = find(result_256 > x(level(2)));

% Check which class is the background (assumed to be the one with the
% most pixels in it).
[bgnum,bgclass] = max([length(ind1) length(ind2) length(ind3)]);
switch bgclass
    case 1 %background is the class with lowest values - assign other two classes as conductive inclusions
        reconstruction([ind2]) = 2;
        reconstruction([ind3]) = 2;
    case 2 %background is the middle class - assign the lower class as resistive inclusions and the higher class as conductive
        reconstruction([ind1]) = 1;
        reconstruction([ind3]) = 2;
    case 3 %background is the class with highest values - assign the other two classes as resistive inclusions
        reconstruction([ind1]) = 1;
        reconstruction([ind2]) = 1;
end

% Show Otsu thresholded result
figure('visible','off')
imagesc(reconstruction)
axis equal
axis off
%end
% newVarName = ['reconstruction' num2str(sample)];
% eval([newVarName ' = reconstruction;']);

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
    
save([outputFolder '/' num2str(sample) '.mat'], 'reconstruction');

% texpmat = reshape(texp,M,M);
% 
% figure
% imagesc(real(texpmat))
% title('real part of texp');
% colorbar
% 
% figure
% imagesc(imag(texpmat))
% title('imag part of texp');
% colorbar
% runtime = toc


end %end loop
disp("Execution Done")
end %end function



