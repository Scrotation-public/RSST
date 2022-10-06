% RSST Final Analysis User code
% last edit 10/4

function [sst,f]  = rsst(x,waveParameter,varargin)
narginchk(1,8);
nbSamp = numel(x);
x = x(:)';
validateattributes(x,{'double'},{'row','finite','real'},'wsst','X');
if numel(x)<4
    error(message('Wavelet:synchrosqueezed:NumInputSamples'));
end
params = parseinputs(nbSamp,waveParameter,varargin{:});
nv = params.nv;
noct = params.noct;
% Create scale vector
na = noct*params.nv;


% If sampling frequency is specified, dt = 1/fs
if (isempty(params.fs) && isempty(params.Ts))
    % The default is 1 for normalized frequency
    dt = params.dt;
    Units = '';
elseif (~isempty(params.fs) && isempty(params.Ts))
    % Accept the sampling frequency in hertz
    fs = params.fs;
    dt = 1/fs;
    Units = '';
elseif (isempty(params.fs) && ~isempty(params.Ts))
    % Get the dt and Units from the duration object
    [dt,Units] = wavelet.internal.getDurationandUnits(params.Ts);
    
    
end

a0 = 2^(1/nv);
scales = a0.^(1:na);
NbSc = numel(scales);

% Construct time series to analyze, pad if necessary
meanSIG = mean(x);
x = x - meanSIG;
NumExten = 0;

if params.pad
    %Pad the time series symmetrically
    np2 = nextpow2(nbSamp);
    NumExten = 2^np2-nbSamp;
    x = wextend('1d','symw',x,NumExten,'b');
end

%Record data length plus any extension
N = numel(x);

%Create frequency vector for CWT computation
omega = (1:fix(N/2));
omega = omega.*((2.*pi)/N);
omega = [0., omega, -omega(fix((N-1)/2):-1:1)];

% Compute FFT of the (padded) time series
xdft = fft(x);

% inserted sstwaveft
% function [wft,dwft] = sstwaveft(WAV,omega,scales,wavparam)
WAV = params.WAV;
wavparam = params.wavparam;

NbSc = numel(scales);
NbFrq = numel(omega);
wft = zeros(NbSc,NbFrq);

switch WAV
    case 'amor'
        
        cf = wavparam;
        
        for jj = 1:NbSc
            expnt = -(scales(jj).*omega - cf).^2/2.*(omega > 0);
            wft(jj,:) = exp(expnt).*(omega > 0);
        end
        
    case 'bump'
        
        mu = wavparam(1);
        sigma = wavparam(2);
        
        
        for jj = 1:NbSc
            w = (scales(jj)*omega-mu)./sigma;
            expnt = -1./(1-w.^2);
            daughter = exp(1)*exp(expnt).*(abs(w)<1-eps(1));
            daughter(isnan(daughter)) = 0;
            wft(jj,:) = daughter;
        end
end

%Compute derivative
omegaMatrix = repmat(omega,NbSc,1);
dwft = 1j*omegaMatrix.*wft;



% [psift,dpsift]  = sstwaveft(params.WAV,omega,scales,params.wavparam);
psift = wft;
dpsift = dwft;

%Obtain CWT coefficients and derivative
cwtcfs = ifft(repmat(xdft,NbSc,1).*psift,[],2);
dcwtcfs = ifft(repmat(xdft,NbSc,1).*dpsift,[],2);

%Remove padding if any
cwtcfs = cwtcfs(:,NumExten+1:end-NumExten);
dcwtcfs = dcwtcfs(:,NumExten+1:end-NumExten);

%Compute the phase transform
phasetf = imag(dcwtcfs./cwtcfs)./(2*pi);


% Threshold for synchrosqueezing
phasetf(abs(phasetf)<params.thr) = NaN;

% Create frequency vector for output
log2Nyquist = log2(1/(2*dt));
log2Fund = log2(1/(nbSamp*dt));
freq = 2.^linspace(log2Fund,log2Nyquist,na);



% inserted sstalgo function
% function Tx = sstalgo(cwtcfs,phasetf,gamma)

M = size(cwtcfs,1);
N = size(cwtcfs,2);
log2Fund = log2(1/N);
log2Nyquist = log2(1/2);
iRow = real(1 + floor(M/(log2Nyquist-log2Fund)*(log2(phasetf)-log2Fund)));
idxphasetf = find(iRow>0 & iRow<=M & ~isnan(iRow));
% changed gamma to params.thr
idxcwtcfs = find(abs(cwtcfs)>params.thr);
idx = intersect(idxphasetf,idxcwtcfs);
iCol = repmat(1:N,M,1);
Tx1 = accumarray([iRow(idx) iCol(idx)],cwtcfs(idx),size(cwtcfs));




Tx = 1/nv*Tx1;


if (nargout == 0)
    
    % inserted plotsst
% plotsst(Tx,freq,dt,params.engunitflag,params.normalizedfreq,Units);
% function plotsst(Tx,F,dt,engunitflag,isfreqnormalized,Units)
freq = F;
engunitflag = params.engunitflag;
isfreqnormalized = params.normalizedfreq;

if ~isempty(Units)
    freqUnits = Units(1:end-1);    
end

t = 0:dt:(size(Tx,2)*dt)-dt;
if engunitflag && isfreqnormalized
    frequnitstrs = wavelet.internal.wgetfrequnitstrs;
    freqlbl = frequnitstrs{1};
    xlbl = 'Samples';
elseif engunitflag && ~isfreqnormalized
    [F,~,uf] = engunits(F,'unicode');
    freqlbl = wavelet.internal.wgetfreqlbl([uf 'Hz']);
    [t,~,ut] = engunits(t,'unicode','time');
    xlbl = [getString(message('Wavelet:getfrequnitstrs:Time')) ' (' ut ')'];
    
else
    freqlbl = getString(message('Wavelet:synchrosqueezed:FreqLabel'));
    freqlbl = ...
        [freqlbl '/' freqUnits ')'];
    xlbl = getString(message('Wavelet:synchrosqueezed:Time'));
    xlbl = [xlbl ' (' Units ')'];
end


h = pcolor(t,F,abs(Tx));
h.EdgeColor = 'none';
shading interp;
ylabel(freqlbl); xlabel(xlbl);
title(getString(message('Wavelet:synchrosqueezed:SynchrosqueezedTitle')));
else
    sst = Tx;
    f = freq;
end





%-------------------------------------------------------------------------
function params = parseinputs(nbSamp,waveParameter,varargin)
% Set defaults.
params.fs = [];
params.dt = 1;
params.Ts = [];
params.sampinterval = false;
params.engunitflag = true;
params.WAV = 'amor';
params.wavparam = waveParameter;
params.thr = 1e-8;
params.nv = 32;
params.noct = floor(log2(nbSamp))-1;
params.pad = false;
params.normalizedfreq = true;

[varargin{:}] = convertStringsToChars(varargin{:});
% Error out if there are any calendar duration objects
tfcalendarDuration = cellfun(@iscalendarduration,varargin);
if any(tfcalendarDuration)
    error(message('Wavelet:FunctionInput:CalendarDurationSupport'));
end

tfsampinterval = cellfun(@isduration,varargin);

if (any(tfsampinterval) && nnz(tfsampinterval) == 1)
    params.sampinterval = true;
    params.Ts = varargin{tfsampinterval>0};
    if (numel(params.Ts) ~= 1 ) || params.Ts <= 0 || isempty(params.Ts)
        error(message('Wavelet:FunctionInput:PositiveScalarDuration'));
    end
    
    params.engunitflag = false;
    params.normalizedfreq = false;
    varargin(tfsampinterval) = [];
end

%Look for Name-Value pairs
numvoices = find(strncmpi('voicesperoctave',varargin,1));

if any(numvoices)
    params.nv = varargin{numvoices+1};
    %validate the value is logical
    validateattributes(params.nv,{'numeric'},{'positive','scalar',...
        'even','>=',10,'<=',48},'wsst','VoicesPerOctave');
    varargin(numvoices:numvoices+1) = [];
    if isempty(varargin)
        return;
    end
end


extendsignal = find(strncmpi('extendsignal',varargin,1));

if any(extendsignal)
    params.pad = varargin{extendsignal+1};
    
    if ~isequal(params.pad,logical(params.pad))
        error(message('Wavelet:FunctionInput:Logical'));
    end
    varargin(extendsignal:extendsignal+1) = [];
    if isempty(varargin)
        return;
    end
end


% Only scalar left must be sampling frequency or sampling interval
% Only scalar left must be sampling frequency
tfsampfreq = cellfun(@(x) (isscalar(x) && isnumeric(x)),varargin);

if (any(tfsampfreq) && (nnz(tfsampfreq) == 1) && ~params.sampinterval)
    params.fs = varargin{tfsampfreq};
    validateattributes(params.fs,{'numeric'},{'positive'},'wsst','Fs');
    params.normalizedfreq = false;
    params.engunits = true;
elseif any(tfsampfreq) && params.sampinterval
    error(message('Wavelet:FunctionInput:SamplingIntervalOrDuration'));
elseif nnz(tfsampfreq)>1
    error(message('Wavelet:FunctionInput:Invalid_ScalNum'));
end

%Only char variable left must be wavelet
tfwav = cellfun(@(x)ischar(x),varargin);
if (nnz(tfwav) == 1)
    params.WAV = varargin{tfwav>0};
    params.WAV = validatestring(params.WAV,{'bump','amor'},'wsst','WAV');
elseif nnz(tfwav)>1
    error(message('Wavelet:FunctionInput:InvalidChar'));
    
end

if strncmpi(params.WAV,'bump',1)
    params.wavparam = [5 1];
end

