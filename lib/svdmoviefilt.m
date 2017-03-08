function [mDenoised, tModes] = svdmoviefilt(mMovieRaw, nMcomp)
% SVDMOVIEFILT Filter movie by singular value decomposition
%   [OUT, Md] = moviesvd(MOV, M), where MOV is a X-by-Y-by-FRAMES matrix.
%   M is maximum number of SVD modes to compute.
%
%   OUT is a matrix of same size as MOV, containing only Md modes that
%   account for more variability in the data than estimated from shuffled
%   data (null hypothesis).
%
% Notes:
%   Lines in the pixels vs time plot may be due to filtering of frames,
%   i.e. is not a bug in this code.
%

%   Copyright, Per M Knutsen <pmknutsen@gmail.com>, 2015
%


% TODO
%   Add as option to return 99% var data

nK = min([50 min(size(mMovieRaw))]); % largest singular value computed

% Parse inputs
bDebug = 1;
nNumFrames = size(mMovieRaw, 3);
vImgSize = [size(mMovieRaw, 1) size(mMovieRaw, 2)];

% Reshape data from 3D to 2D by concatenating images as row-vectors across columns and time across rows
vPixels = mMovieRaw(:);
mReshaped = reshape(vPixels, [prod(vImgSize) nNumFrames])';
mShuffled = reshape(vPixels(randperm(length(vPixels))), [prod(vImgSize) nNumFrames])';

% Compute SVD of data
disp('Computing SVDs of data (null hypothesis)...')
[u,s,v] = svds(mReshaped, nK);

% Compute SVD of shuffled data (for null hypothesis)
disp('Computing SVDs of shuffled data (null hypothesis)...')
[~,ss,~] = svds(mShuffled, min([nK size(s, 1)]));

% Get the N most significant modes (compared to null hypothesis)
nMdisp = find(diag(s.^2) > diag(ss.^2), 1, 'last');
vEigenvalues = diag(s);
vEigenProp = vEigenvalues ./ sum(vEigenvalues) * 100;

% Get data on modes that is returned
tModes.number_of_modes = nMdisp;
tModes.eigenprop_sum = sum(vEigenProp(1:nMdisp)); % total percent variability accounted for

% Denoise by keeping only first nMdisp modes
mReshapedDenoised = u(:,1:nMdisp) * s(1:nMdisp, 1:nMdisp) * v(:, 1:nMdisp)';

% Reshape denoised data back to original dimensions
mDenoised = reshape(mReshapedDenoised', [vImgSize nNumFrames]);

%% Display singular values
if bDebug
    vClim = [min(mMovieRaw(:)) max(mMovieRaw(:))];
    
    hFig = figure('Name', 'Denoised Time Series', 'numbertitle', 'off');
    hAx(1) = subplot(1, 3, 1);
    hold(hAx(1), 'on')
    plot(diag(s.^2), 'ko-.', 'markerSize', 12);
    plot(diag(ss.^2), 'ro-.', 'markerSize', 12);
    axis(hAx(1), 'tight')
    set(hAx(1), 'yscale', 'log')
    ylabel('\lambda_i^2') % eigenvalue, lambda square
    xlabel('Ranked index (i)') % index of eigenvalue (mode number)

    hAx = [];
    hAx(end+1) = subplot(1, 3, 2);
    hImg(1) = imagesc(mReshaped);
    title('Original')
    
    hAx(end+1) = subplot(1, 3, 3);
    hImg(end+1) = imagesc(mReshapedDenoised);
    title(sprintf('Reconstructed with %d modes', nMdisp));
    
    set(hAx, 'clim', vClim)
    axis(hAx, 'off')
    colormap(hAx(1), copper)

    % Display first Md modes
    hFig = figure('name', 'Rank-ordered modes', 'numbertitle', 'off');
    nspcols = ceil(sqrt(nMcomp));
    nsprows = ceil(nMcomp / nspcols);
    
    % TODO DISPLAY ASO TIMES MODES!!
        
    ny = vImgSize(1);
    nx = vImgSize(2);
    nt = nNumFrames;
    
    hAx = [];
    for iSP = 1:nMcomp
        hAx(end+1) = subplot(nsprows, nspcols, iSP);
        imagesc(reshape(v(:, iSP), ny, nx));
        axis image off
        title(sprintf('%3.2f %%', vEigenProp(iSP)))
    end
    set(hAx, 'clim', [-.1 .1])
    colormap(hAx(1), copper)
end
%%

return

