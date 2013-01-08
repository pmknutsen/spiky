%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    % Calculate an estimate of the crosscorrelation    %
%    % function for either local field or spikes        %
%    %                                                  %
%    % Note on methods used for spikes...               %
%    % Uses standard histogram method.  Use method 2 to %
%    % remove the shuffle correction i.e. the auto-corr %
%    % recalculated with the trials shuffled.  This is  %
%    % similar to, but not the same as removing the     %
%    % value expected under a Poisson null (as the rate %
%    % may vary).  Methods 3 and 4 both subtract the    %
%    % expected values for a homogemeous Poisson null   %
%    % Method 3 divides by the standard deviation       %
%    % expected under a Poisson null                    %
%    %                                                  %
%    % In keeping with convension in neuroscience we use%
%    % non area preserving auto-cross correlations i.e. %
%    % they are based on counts of intervals rather and %
%    % are not normalized by the data length            %
%    %                                                  %
%    % Note that the cross correlation between spikes   %
%    % and the local field is the spike triggered ave   %
%    % and so this is treated separately (by sta.m)     %
%    %                                                  %
%    % A common clock smp is assumed for lfp data       %
%    %                                                  %
%    % INPUT                                            %
%    %                                                  %
%    % data1   - spike or local field (required)        %
%    % data2   - spike or local field (required)        %
%    % smp    - lfp sample times |range of tau to plot  %
%    % plt    - plot ('n' | color)    default 'n'       %
%    % T      - time interval (default all)             %
%    % method  /1 = biased (normal autocorrelation)     %
%    % (if lfp)|2 = unbiased (divide by triangle fn)    %
%    %         \3 = correlation coeffient (1 at 0)      %
%    %         /1 = basic (normal autocorrelation)      %
%    % (if spk)|2 = shuffle correction subtracted       %
%    %         |3 = homogeneous poisson subtracted      %
%    %         \4 = hom poisson subtracted and norm -   %
%    %              by the standard deviation under     %
%    %              Poisson null                        %
%    % if -ve then do an adaptively smoothed kernel     %
%    % estimate rather than a histogram                 %
%    %                                                  %
%    % err    - 1   errorbars                           %
%    %        - 0   no errorbars                        %
%    %                                                  %
%    % Nbins  number of bins (default -1 = auto)        %
%    %        or mean std dev of kernel estimate if     %
%    %        method if -ve                             %
%    %                                                  %
%    % OUTPUT                                           %
%    %                                                  %
%    % mu       function                                %
%    % tau      lag times                               %
%    % E        error bar (standard deviation)          %
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[mu,tau,E] = ccorr(data1,data2,smp,plt,T,method,err,Nbins)

if isempty(data1); error('I need data'); end
if isempty(data2); error('I need data'); end

flag1 = 0;if length(find(diff(data1(1,:))<0)) > 3;flag1 = 1;end
flag2 = 0;if length(find(diff(data2(1,:))<0)) > 3;flag2 = 1;end

if flag1 ~= flag2; error(['data1 and data2 should be the ' ...
        'same data type : use sta instead']);
else flag = flag1;end

if nargin < 4; plt = 'n'; end
if nargin < 5
    if flag; T = [min(smp) max(smp)];
    else; T = [min(data1(:,1)) max(max(data1))]; end
end
if nargin < 6; method = 1;end
if nargin < 7; err = 0; end
if nargin < 8; Nbins = -1;end

if isempty(plt); plt = 'n'; end
if isempty(T);if flag; T = [min(smp) max(smp)];else; T = [min(data1(:,1)) max(max(data1))]; end;end
if isempty(method); method = 1;end
if isempty(err); err = 0; end
if isempty(Nbins); Nbins = -1;end

if flag == 0 & isempty(smp); smp = [-(T(2)-T(1)) (T(2)-T(1))];end
if nargin < 3 & flag == 1; error('sample times required for lfp');end

NT = length(data1(:,1));
if NT < 4; err = 0; end

if flag                           % continuous process data...
    indx = find(smp>T(1) & smp < T(2));
    smp = smp(2)-smp(1);
    data1 = data1(:,indx);
    data2 = data2(:,indx);
    N = length(data1(1,:));
    ccorr = zeros(NT,2*N-1);
    for n=1:NT
        data1(n,:) = detrend(data1(n,:));
        data2(n,:) = detrend(data2(n,:));
        tau = smp*(-(N-1):(N-1));
        if method == 1;              %  biased
            ccorr(n,:) = xcorr(data1(n,:),data2(n,:),'biased');
        elseif method == 2           % unbiased
            ccorr(n,:) = xcorr(data1(n,:),data2(n,:),'unbiased');
        elseif method == 3           % coeff
            ccorr(n,:) = xcorr(data1(n,:),data2(n,:),'coeff');
        end
    end
else                             % spike data...

    %  trim data...

    datatmp = zeros(size(data1));
    Np1 = zeros(1,NT);
    for n=1:NT
        indx = find(T(1)<data1(n,:) & T(2)>data1(n,:) & data1(n,:) ~= 0);
        Np1(n) = length(indx);
        datatmp(n,1:Np1(n)) = data1(n,indx) - T(1);
    end
    data1 = datatmp(:,1:max(Np1));

    datatmp = zeros(size(data2));
    Np2 = zeros(1,NT);
    for n=1:NT
        indx = find(T(1)<data2(n,:) & T(2)>data2(n,:) & data2(n,:) ~= 0);
        Np2(n) = length(indx);
        datatmp(n,1:Np2(n)) = data2(n,indx) - T(1);
    end
    data2 = datatmp(:,1:max(Np2));

    % figure out the bin density and grid...

    if Nbins == -1
        Nbins = max(floor((sum(Np1.*Np2) - NT)/20),2);
    end

    DT = T(2)-T(1);
    spcfac = 5;            % density of bins for kernel smoother (bins per sig)
    if method > 0; Ngrid = Nbins;end
    if method < 0; Ngrid = spcfac*Nbins;end

    tau = linspace(smp(1),smp(2),Ngrid);
    ccorr = zeros(NT,Ngrid);
    mn_N = zeros(NT,Ngrid);
    std_N = zeros(NT,Ngrid);
    spc = tau(2)-tau(1);

    for n=1:NT
        intervals = getintervals(data1(n,1:Np1(n)),data2(n,1:Np2(n)));
        indx = find(intervals > smp(1) & intervals < smp(2));
        intervals = intervals(indx);
        if method > 0
            ccorr(n,:) = hist(intervals,tau);
            %      plot(tau,ccorr(n,:),'r-');hold on
        else
            ccorr(n,:)= psth(intervals,-((smp(2)-smp(1))/Nbins),'n',[smp(1) smp(2)],0,tau);
            fac = spcfac*length(intervals)/sum(ccorr(n,:)); % since grid is spcfac x denser
            ccorr(n,:) = fac*ccorr(n,:);                    % for the kernel smoothed
            %      plot(tau,ccorr(n,:),'r-');hold on
        end

        %   different normalizations...

        if method < 0; spc0 = spc*spcfac; else; spc0 = spc;end
        mn_N(n,:) = length(intervals)*(DT-abs(tau))*spc0/(DT^2);

        std_N(n,:) = sqrt(mn_N(n,:));
        if abs(method) == 3                     % subtract the mean under
            ccorr(n,:) = ccorr(n,:) - mn_N(n,:);  % Poisson null
        elseif abs(method) == 4
            ccorr(n,:) = (ccorr(n,:) - mn_N(n,:))./std_N(n,:);
        end
    end

    if abs(method) == 2  % calc cross-corr with other trials...
        shuffle_corr = 0;
        for n = 1:NT
            for nn=n+1:NT
                intervals = ...
                    getintervals(data1(n,1:Np1(n)),data2(nn,1:Np2(nn)));
                indx = find(intervals > smp(1) & intervals < smp(2));
                intervals = intervals(indx);
                if method > 0
                    shcorr = hist(intervals,tau);
                else
                    shcorr= psth(intervals,-((smp(2)-smp(1))/Nbins),'n',[smp(1) smp(2)],0,tau);
                    fac = spcfac*length(intervals)/sum(shcorr);  % since the grid is 5x denser
                    shcorr = fac*shcorr;                         % for the kernel smoothed
                end
                shuffle_corr = shuffle_corr + shcorr;
            end
        end
        shuffle_corr = shuffle_corr/(NT*(NT-1));
        for n = 1:NT
            ccorr(n,:) = ccorr(n,:) - shuffle_corr;
        end
    end
end
if NT > 1;mu = mean(ccorr); else; mu = ccorr;end

% do the errorbars...
if err
    Nboot = 20;
    m = 0;
    s = 0;
    for b=1:Nboot
        indx = floor(NT*rand(1,NT)) + 1;
        tmp = mean(ccorr(indx,:));
        m = m + tmp;
        s = s + tmp.^2;
    end
    E = sqrt((s/Nboot - m.^2/Nboot^2));
else
    E = [];
end

if ~strcmp(plt,'n')
    if flag
        plot(tau,mu,plt)
        hold on
        if err == 1
            plot(tau,mu + 2*E,'g')
            plot(tau,mu - 2*E,'g')
        end
        line(get(gca,'xlim'),[0 0],'color','k')
        line([0 0],get(gca,'ylim'),'color','k')
        hold off
    else
        if method > 0
            bar(tau,mu);
        else
            plot(tau,mu,'r-');
        end
        line([0 0],get(gca,'ylim'),'color','k')
        hold on
        if err == 1
            plot(tau,2*E,'g')
        end
        hold off
    end
end

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Routine used by the auto and cross correlation routines%
%     % Returns all possible intervals between X and Y         %
%     % Exact zeros are excluded                               %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[int] = getintervals(A,B)

NA = length(A);
NB = length(B);

NaNb1 = repmat(A,NB,1);
NaNb2 = repmat(B',1,NA);
int = reshape(NaNb1 - NaNb2,1,NA*NB);
int = sort(int(find(int ~=0)));
