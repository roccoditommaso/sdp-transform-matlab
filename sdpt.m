function varargout = sdpt(x, fmin, fmax, fs, varargin)
%SDPT  Stockwell-Ditommaso-Ponzo (SDP) time-frequency representation.
%
%   ---------------------------------------------------------------------------
%   QUICK GUIDE
%   ---------------------------------------------------------------------------
%   S = SDPT(X,FMIN,FMAX,FS)          magnitude between FMIN and FMAX [Hz]
%   [S,F,T] = SDPT(...)               nodes F [Hz] and time T [s]
%   [S,F,T,INFO] = SDPT(...)          + settings, rational labels, ridge
%   Any call plots the record and the map unless 'Plot' is set to false.
%   SDPT                              no input: runs a self-contained demo
%
%   Every option is a name-value pair and every default is collected in the
%   block DEFAULT SETTINGS at the top of the code, immediately below this help.
%   Change a default there once, or override it in the call; you never need to
%   read or edit the body of the function.
%
%     NAME            DEFAULT       MEANING
%     'Conductor'     [] (auto)     D, largest admissible denominator. It fixes
%                                   the support together with (FS,[FMIN FMAX])
%                                   and is the only knob for resolution. Auto
%                                   means: as many nodes as the conventional
%                                   Stockwell transform on the same band, and
%                                   never fewer than needed to reach FMIN.
%                                   Declare it explicitly when several records
%                                   must share the same support.
%     'Width'         1             k > 0, Gaussian window width in units of the
%                                   original Stockwell window (sigma_t = k/f).
%                                   k < 1 sharpens time, k > 1 sharpens frequency.
%     'Ridge'         true          extract the instantaneous-frequency ridge.
%                                   Set false to skip it (slightly faster).
%     'RidgePenalty'  0.5           smoothness of the ridge, >= 0, in units of
%                                   (nodes per stored time sample)^2. 0 gives the
%                                   plain maximum of each column; larger values
%                                   forbid abrupt jumps. The default costs
%                                   nothing on a clean signal and roughly halves
%                                   the ridge error at 0 dB; raise it for noisy
%                                   records with a slowly moving component,
%                                   lower it when the frequency changes fast.
%     'RidgeBand'     []            [flo fhi] restricts the ridge search to a
%                                   sub-band, without changing the map.
%     'Output'        'magnitude'   'magnitude' or 'complex' (keeps the absolute
%                                   phase reference of the Stockwell transform).
%     'TimeStep'      1             store one time sample every TimeStep samples.
%     'MaxNodes'      20000         safety limit on the number of nodes.
%     'Block'         128           nodes evaluated per batch: speed and memory
%                                   only, never the result.
%     'Plot'          true          show the figure: the record on top and the
%                                   SDP map with its ridge below, on time axes of
%                                   identical geometry and linked to each other.
%                                   Set false to compute without plotting.
%
%   ---------------------------------------------------------------------------
%   WHAT THIS TRANSFORM IS FOR
%   ---------------------------------------------------------------------------
%   The SDP representation is a Stockwell transform evaluated on a frequency
%   support declared before the analysis and independent of the length of the
%   record: the primitive fractions a/d with d <= D (a Farey sequence) contained
%   in the analysed band. Two records of the same process, of different
%   duration, therefore produce coefficients that occupy the same positions and
%   can be compared, stored and differenced directly, with no interpolation onto
%   a common grid. The support is fixed by (FS, D, [FMIN FMAX]) alone.
%
%   ---------------------------------------------------------------------------
%   INPUT / OUTPUT
%   ---------------------------------------------------------------------------
%     X      signal: a vector, or a matrix with samples along the rows and
%            channels along the columns (a row vector is transposed).
%     FMIN   lowest frequency to be represented [Hz], > 0.
%     FMAX   highest frequency to be represented [Hz], <= FS/2.
%     FS     sampling frequency [Hz].
%
%     S      M-by-NT, or M-by-NT-by-C for C channels. Singleton dimensions are
%            dropped, as MATLAB does with any array.
%     F      M-by-1 frequency nodes [Hz], ascending and not uniformly spaced.
%     T      1-by-NT time vector [s].
%     INFO   struct with fields:
%              ridge      NT-by-C ridge frequency [Hz] ([] if 'Ridge' is false)
%              ridgeIdx   NT-by-C index of the ridge node
%              a, d, q    rational label of each node (F = a/d*FS)
%              edges      M+1 Voronoi cell edges [Hz]
%              weights    M-by-1 cell widths df/f, for band integration
%              D, k, band, fs, N, channels, nodes, timeStep, output
%
%   ---------------------------------------------------------------------------
%   EXAMPLES
%   ---------------------------------------------------------------------------
%     %% a chirp, everything by default (map + ridge)
%     fs = 200; t = (0:1/fs:10)'';
%     x  = sin(2*pi*(5*t + 0.5*t.^2));
%     sdpt(x,2,30,fs);
%
%     %% read the ridge back
%     [~,~,t,info] = sdpt(x,2,30,fs,'Plot',false);
%     plot(t,info.ridge); xlabel('Time [s]'); ylabel('Ridge [Hz]');
%
%     %% two channels, declared support, decimated storage, no ridge
%     y = [x, sin(2*pi*12*t)];
%     [S,f,t] = sdpt(y,2,30,fs,'Conductor',150,'TimeStep',4,'Ridge',false);
%
%     %% a stiffer ridge, searched only below 20 Hz
%     sdpt(x,2,30,fs,'RidgePenalty',2,'RidgeBand',[2 20]);
%
%   ---------------------------------------------------------------------------
%   REFERENCE
%   ---------------------------------------------------------------------------
%     Ditommaso R., Ponzo F.C. (2026). Time-frequency representation of
%     nonlinear structural response through rationally sampled Stockwell
%     analysis. Mechanical Systems and Signal Processing.
%     DOI: https://doi.org/10.1016/j.ymssp.2026.114920
%
%   A short illustrated guide (SDPT_guide.pdf) is distributed with the function.
%
%   Please cite the reference above if you use this function in your research.
%
%   FUNDING
%     This study was funded by the Italian Civil Protection Department within the
%     project DPC-ReLUIS 2024-26 - WP6 "Monitoring and satellite data".
%
%   Authors: Rocco Ditommaso, Felice Carlo Ponzo
%            Department of Engineering, University of Basilicata, Potenza, Italy
%
%   LIMITATIONS
%     The convolution is performed with the FFT and is therefore circular, as in
%     the usual implementations of the Stockwell transform. The Gaussian window
%     spans about 8k/f seconds, so near FMIN a transient close to one end of the
%     record leaks into the other end. Keep FMIN well above 8k/T, where T is the
%     record duration, or disregard the first and last few window lengths; a
%     warning is issued when the window becomes comparable to the record. The
%     same applies to the ridge, which is extracted from the map as it is.
%
%   REQUIREMENTS
%     MATLAB is required. No toolbox is required: the transform uses only
%     core functions. The colormap falls back to PARULA where TURBO is absent.
%
%   See also FFT, SPECTROGRAM, CWT.

% Copyright (c) 2026, Rocco Ditommaso and Felice Carlo Ponzo
% Licensed under the BSD 3-Clause License; see LICENSE.txt.
% See README.md for the associated research article.

%% ===========================================================================
%  DEFAULT SETTINGS  --  the only block a user may want to edit
%  Each entry is the default of the name-value pair with the same name, and can
%  always be overridden in the call without editing anything here.
%  ============================================================================
DEF.Conductor    = [];            % D. [] = automatic (see QUICK GUIDE)
DEF.Width        = 1;             % k, window width; 1 = original Stockwell window
DEF.Ridge        = true;          % extract the ridge
DEF.RidgePenalty = 0.5;           % ridge smoothness, >= 0 (0 = column maximum)
DEF.RidgeBand    = [];            % [flo fhi] to restrict the ridge search
DEF.Output       = 'magnitude';   % 'magnitude' or 'complex'
DEF.TimeStep     = 1;             % store one sample every TimeStep samples
DEF.MaxNodes     = 20000;         % safety limit on the number of nodes
DEF.Block        = 128;           % nodes per batch (speed/memory only)
DEF.Plot         = true;          % show the figure (record + map); false to skip
%  ============================================================================
%  END OF SETTINGS  --  nothing below needs to be edited to use the function
%  ============================================================================

if nargin == 0, localDemo(); return, end
narginchk(4,inf);

%% ---------------------------------------------------------------- validation
validateattributes(x   ,{'numeric'},{'nonempty','real','finite','2d'},mfilename,'X',1);
validateattributes(fmin,{'numeric'},{'scalar','real','positive','finite'},mfilename,'FMIN',2);
validateattributes(fmax,{'numeric'},{'scalar','real','positive','finite'},mfilename,'FMAX',3);
validateattributes(fs  ,{'numeric'},{'scalar','real','positive','finite'},mfilename,'FS',4);

if isrow(x), x = x(:); end
x = double(x); [N,C] = size(x);

if fmax <= fmin
    error('sdpt:band','FMAX (%g Hz) must be larger than FMIN (%g Hz).',fmax,fmin);
end
if fmax > fs/2 + eps(fs)
    error('sdpt:nyquist','FMAX (%g Hz) exceeds the Nyquist frequency (%g Hz).',fmax,fs/2);
end
if fmin < fs/N
    error('sdpt:tooLow',['FMIN (%g Hz) is below the frequency resolution of the ' ...
          'record (fs/N = %g Hz). Use a longer record or raise FMIN.'],fmin,fs/N);
end

p = inputParser; p.FunctionName = mfilename;
addParameter(p,'Conductor'   ,DEF.Conductor   ,@(v)isempty(v)||(isscalar(v)&&isreal(v)&&v>=2&&isfinite(v)));
addParameter(p,'Width'       ,DEF.Width       ,@(v)isscalar(v)&&isreal(v)&&v>0&&isfinite(v));
addParameter(p,'Ridge'       ,DEF.Ridge       ,@(v)isscalar(v)&&(islogical(v)||isnumeric(v)));
addParameter(p,'RidgePenalty',DEF.RidgePenalty,@(v)isscalar(v)&&isreal(v)&&v>=0&&isfinite(v));
addParameter(p,'RidgeBand'   ,DEF.RidgeBand   ,@(v)isempty(v)||(numel(v)==2&&v(2)>v(1)));
addParameter(p,'Output'      ,DEF.Output      ,@(v)any(strcmpi(v,{'magnitude','complex'})));
addParameter(p,'TimeStep'    ,DEF.TimeStep    ,@(v)isscalar(v)&&v>=1&&v==fix(v));
addParameter(p,'MaxNodes'    ,DEF.MaxNodes    ,@(v)isscalar(v)&&v>=1);
addParameter(p,'Block'       ,DEF.Block       ,@(v)isscalar(v)&&v>=1&&v==fix(v));
addParameter(p,'Plot'        ,DEF.Plot        ,@(v)isempty(v)||islogical(v)||isnumeric(v));
parse(p,varargin{:}); opt = p.Results;

doPlot = opt.Plot; if isempty(doPlot), doPlot = true; end
doPlot     = logical(doPlot);
doRidge    = logical(opt.Ridge);
complexOut = strcmpi(opt.Output,'complex');

%% ------------------------------------------------- rational support (Farey)
b1 = fmin/fs;  b2 = fmax/fs;
Dmin = ceil(1/b1);                    % below this, FMIN has no node at all
if isempty(opt.Conductor)
    D = max(Dmin, round(pi*sqrt(N/3)));   % ~ as many nodes as the FFT grid
else
    D = round(opt.Conductor);
    if D < Dmin
        warning('sdpt:conductorTooSmall', ...
            ['Conductor D = %d cannot represent FMIN = %g Hz: the lowest node of ' ...
             'the support is fs/D = %g Hz. Raising D to %d.'],D,fmin,fs/D,Dmin);
        D = Dmin;
    end
end

[q,a,d] = localFarey(D,b1,b2);
M = numel(q);
if M == 0
    error('sdpt:emptySupport',['The support is empty. Increase Conductor or ' ...
          'widen the band [%g %g] Hz.'],fmin,fmax);
end
if M > opt.MaxNodes
    error('sdpt:tooManyNodes',['The support has %d nodes, above MaxNodes = %d. ' ...
          'Lower Conductor (now %d), narrow the band, or raise MaxNodes.'],M,round(opt.MaxNodes),D);
end
f = q(:)*fs;

%% ------------------------------------------------ Voronoi cells and weights
mid = (f(1:end-1) + f(2:end))/2;
if M == 1
    edges = [fmin; fmax];
else
    edges = [max(fmin,f(1)-(mid(1)-f(1))); mid; min(fmax,f(end)+(f(end)-mid(end)))];
end
weights = (edges(2:end) - edges(1:end-1))./f;      % df/f, the Stockwell measure

if 8/fmin > N/fs
    warning('sdpt:windowLongerThanRecord', ...
        ['At FMIN = %g Hz the Gaussian window spans %.3g s, comparable to the ' ...
         'record (%.3g s). Because the convolution is circular, the lowest ' ...
         'frequencies will be contaminated by wrap-around; raise FMIN or use a ' ...
         'longer record.'],fmin,8*opt.Width/fmin,N/fs);
end

%% --------------------------------------------------------------- transform
tIdx = 1:opt.TimeStep:N;  NT = numel(tIdx);
bytes = M*NT*C*8*(1+complexOut);
if bytes > 2e9
    warning('sdpt:largeOutput',['The output needs about %.1f GB. Consider a ' ...
            'smaller Conductor or a larger TimeStep.'],bytes/1e9);
end

if complexOut, S = complex(zeros(M,NT,C)); else, S = zeros(M,NT,C); end
n  = (0:N-1)';
kk = [0:ceil(N/2)-1, -floor(N/2):-1]';
nu = kk/N;  k2 = opt.Width^2;
blk = min([opt.Block, M, max(1,floor(3e7/N))]);   % cap the N-by-blk intermediate

for c = 1:C
    xc = x(:,c);
    for i0 = 1:blk:M
        j  = i0:min(i0+blk-1,M);
        qb = reshape(q(j),1,[]);
        % exact spectral shift: the modulation is exact for any rational node,
        % so the coefficient never depends on the FFT grid of the record
        Y = fft(xc.*exp(-2i*pi*(n*qb)));
        G = exp(-2*pi^2*k2*(nu.^2)*(1./qb.^2));   % Stockwell window, G(0) = 1
        V = ifft(Y.*G);
        V = V(tIdx,:).';
        if complexOut, S(j,:,c) = V; else, S(j,:,c) = abs(V); end
    end
end

%% ------------------------------------------------------------------- ridge
ridge = []; ridgeIdx = [];
if doRidge
    sel = true(M,1);
    if ~isempty(opt.RidgeBand)
        sel = f >= opt.RidgeBand(1) & f <= opt.RidgeBand(2);
        if ~any(sel)
            warning('sdpt:emptyRidgeBand','RidgeBand contains no node: using the full band.');
            sel = true(M,1);
        end
    end
    fsel = f(sel);
    ridgeIdx = zeros(NT,C); ridge = zeros(NT,C);
    for c = 1:C
        ii = localRidge(abs(S(sel,:,c)), opt.RidgePenalty);
        ridgeIdx(:,c) = ii(:); ridge(:,c) = fsel(ii);
    end
    idxAll = find(sel); ridgeIdx = idxAll(ridgeIdx);   % indices into the full support
    ridgeIdx = reshape(ridgeIdx,NT,C);
end

S = squeeze(S);
t = (tIdx-1)/fs;

info = struct('ridge',ridge,'ridgeIdx',ridgeIdx,'a',a(:),'d',d(:),'q',q(:), ...
              'edges',edges,'weights',weights,'D',D,'k',opt.Width, ...
              'band',[fmin fmax],'fs',fs,'N',N,'channels',C,'nodes',M, ...
              'timeStep',opt.TimeStep,'output',lower(opt.Output), ...
              'ridgePenalty',opt.RidgePenalty, ...
              'reference','Ditommaso & Ponzo (2026), MSSP, doi:10.1016/j.ymssp.2026.114920');

if doPlot, localPlot(S,f,t,info,x(tIdx,1)); end

varargout = {S,f,t,info};
varargout = varargout(1:max(nargout,1));
if nargout == 0, varargout = {}; end
end

%% ========================================================= local functions
function [q,a,d] = localFarey(D,b1,b2)
%LOCALFAREY  Primitive fractions a/d, d <= D, contained in [b1,b2].
cap = 0;
for dd = 1:D, cap = cap + max(0,floor(b2*dd) - ceil(b1*dd) + 1); end
a = zeros(cap,1); d = zeros(cap,1); m = 0;
for dd = 1:D
    lo = max(1,ceil(b1*dd)); hi = floor(b2*dd);
    for aa = lo:hi
        if gcd(aa,dd) == 1, m = m+1; a(m) = aa; d(m) = dd; end
    end
end
a = a(1:m); d = d(1:m);
[q,i] = unique(a./d);
a = a(i); d = d(i);
end

function idx = localRidge(A,lambda)
%LOCALRIDGE  Ridge of a time-frequency map by dynamic programming.
%   Maximises the accumulated amplitude while penalising jumps between
%   consecutive columns by lambda*(delta node index)^2. lambda = 0 returns the
%   plain column-wise maximum. The quadratic penalty is minimised with the
%   lower-envelope algorithm, so the cost is linear in the number of nodes.
[M,NT] = size(A);
if M == 1, idx = ones(NT,1); return, end
A = A./max(A(:) + realmin);
cost = -log(A + 1e-12);                 % small amplitude, large cost
if lambda < 1e-12                       % also catches subnormal values, for which
    [~,idx] = min(cost,[],1);           % the envelope recursion would divide by ~0
    idx = idx(:); return
end
Cacc = cost(:,1);
back = zeros(M,NT,'uint32');
for t = 2:NT
    [Cmin,arg] = localLowerEnvelope(Cacc,lambda);
    back(:,t)  = uint32(arg);
    Cacc       = cost(:,t) + Cmin;
end
idx = zeros(NT,1);
[~,idx(NT)] = min(Cacc);
for t = NT:-1:2, idx(t-1) = double(back(idx(t),t)); end
end

function [Cmin,arg] = localLowerEnvelope(C,lambda)
%LOCALLOWERENVELOPE  min_j C(j) + lambda*(i-j)^2 for every i, in O(M).
%   Felzenszwalb & Huttenlocher generalised distance transform.
M = numel(C);
v = ones(M,1); z = -inf(M+1,1); z(1) = -inf; z(2) = inf; k = 1;
for i = 2:M
    s = ((C(i) + lambda*i^2) - (C(v(k)) + lambda*v(k)^2)) / (2*lambda*(i - v(k)));
    while k > 1 && s <= z(k)
        k = k - 1;
        s = ((C(i) + lambda*i^2) - (C(v(k)) + lambda*v(k)^2)) / (2*lambda*(i - v(k)));
    end
    k = k + 1; v(k) = i; z(k) = s; z(k+1) = inf;
end
Cmin = zeros(M,1); arg = zeros(M,1); k = 1;
for i = 1:M
    while z(k+1) < i, k = k + 1; end
    Cmin(i) = lambda*(i - v(k))^2 + C(v(k));
    arg(i)  = v(k);
end
end

function localPlot(S,f,t,info,xs)
%LOCALPLOT  Record on top, SDP map with its ridge below.
%   The two axes are given explicit positions with identical left edge and
%   width, and the colorbar is placed by hand outside them, so that the two
%   time axes have exactly the same geometry and remain directly comparable.
%   They are also linked, so zooming or panning one moves the other.
A = abs(S);
if ~ismatrix(A)
    warning('sdpt:plotFirstChannel','Only the first of %d channels is displayed.',size(S,3));
    A = A(:,:,1);
end

LEFT = 0.11; WIDTH = 0.74; CBW = 0.025;      % shared geometry
BOT2 = 0.100; H2 = 0.585;                    % map
BOT1 = 0.730; H1 = 0.185;                    % record

fig = figure('Color','w','Name','SDP representation');

ax2 = axes('Parent',fig,'Position',[LEFT BOT2 WIDTH H2]);
dt  = info.timeStep/info.fs;
h = pcolor(ax2,[t(:).', t(end)+dt], info.edges, [A, A(:,end); A(end,:), A(end,end)]);
set(h,'EdgeColor','none'); axis(ax2,'xy'); hold(ax2,'on');
if ~isempty(info.ridge)
    plot(ax2,t,info.ridge(:,1),'w','LineWidth',1.4);
end
try colormap(ax2,turbo); catch, colormap(ax2,parula); end
xlabel(ax2,'Time [s]'); ylabel(ax2,'Frequency [Hz]');
ylim(ax2,info.band);
if numel(t) > 1, XL = [t(1) t(end)]; else, XL = [t(1) t(1)+dt]; end   % a single
xlim(ax2,XL);                                    % stored sample still needs a range
cb = colorbar(ax2); set(get(cb,'Label'),'String','|SDP|');
set(cb,'Position',[LEFT+WIDTH+0.025 BOT2 CBW H2]);
set(ax2,'Position',[LEFT BOT2 WIDTH H2]);   % colorbar must not resize the axes

ax1 = axes('Parent',fig,'Position',[LEFT BOT1 WIDTH H1]);
plot(ax1,t,xs,'Color',[0 0.2 0.55],'LineWidth',0.6);
ylabel(ax1,'x(t)'); grid(ax1,'on');
xlim(ax1,XL); set(ax1,'XTickLabel',[]);
title(ax1,sprintf('SDP  --  D = %d,  %d nodes in [%.3g, %.3g] Hz,  k = %g', ...
      info.D,info.nodes,f(1),f(end),info.k));

linkaxes([ax1 ax2],'x');
set(fig,'UserData',struct('axRecord',ax1,'axMap',ax2,'colorbar',cb));
end

function localDemo()
%LOCALDEMO  Self-contained example, run when sdpt is called with no input.
fs = 200; t = (0:1/fs:10)';
fi = 4 + 6./(1+exp(-(t-5)));                 % 4 Hz -> 10 Hz transition
x  = sin(2*pi*cumtrapz(t,fi)) + 0.05*randn(size(t));
fprintf('sdpt demo: %g s at %g Hz, instantaneous frequency from 4 to 10 Hz.\n',t(end),fs);
[~,f,tt,info] = sdpt(x,2,20,fs);             % the standard figure appears here
fprintf('support: D = %d, %d nodes between %g and %g Hz.\n',info.D,info.nodes,f(1),f(end));
err = 1000*sqrt(mean((info.ridge - interp1(t,fi,tt(:))).^2));
fprintf('ridge RMSE against the true instantaneous frequency: %.1f mHz\n',err);
figure('Color','w','Name','sdpt demo -- ridge check');
plot(tt,info.ridge,'LineWidth',1.3); hold on; plot(t,fi,'k--','LineWidth',1.1);
xlabel('Time [s]'); ylabel('Frequency [Hz]'); grid on; axis tight
legend('SDP ridge','true instantaneous frequency','Location','best');
title('demo -- extracted ridge against the truth');
end
