% RSST Final Analysis User code
% last edit 10/5
% needs rsst.m

function RSSTAnalysis(plotTitle, currentData, currentFreq, varargin)
    VersionNum = 1005.2022; %monthday.year

%% optional input parameters
    p = inputParser;
    p.FunctionName = 'RSSTAnalysis';
    default_yMin = 0.3;
    default_yMax = 160;
    default_waveParameter = 4;
    default_currentTimeData = [];
    default_listBoxAlgorithm = 'SSWT';
    default_switchUnitType = 'Potential';
    default_switchYScale = 'Log';
    default_useMs = 'Off';
    default_filterAC = 'Off';
    default_removeCOI = 'On';
    default_switchSqHz = 'Off';
    default_boost = 1;
    default_editFieldSensitivity = 1;
    default_savePlotFile = '';
    checkOnOff = @(x) any(validatestring(x,{'On','Off'}));
    checkAlgorithm = @(x) any(validatestring(x,{'SSWT','CWT'}));
    checkUnitType = @(x) any(validatestring(x,{'Power','Potential'}));
    checkYScale = @(x) any(validatestring(x,{'Log','Linear'}));
    addParameter(p,'algorithm', default_listBoxAlgorithm, checkAlgorithm);
    addParameter(p,'unitType', default_switchUnitType, checkUnitType);
    addParameter(p, 'yScale', default_switchYScale, checkYScale);
    addParameter(p, 'sqHz', default_switchSqHz, checkOnOff);
    addParameter(p, 'removeCOI', default_removeCOI, checkOnOff);
    addParameter(p, 'filterAC', default_filterAC, checkOnOff);
    addParameter(p, 'useMs', default_useMs, checkOnOff);
    addParameter(p, 'waveParameter', default_waveParameter, @isnumeric)
    addParameter(p, 'currentTimeData', default_currentTimeData, @isnumeric)
    addParameter(p, 'sensitivity', default_editFieldSensitivity, @isnumeric)
    addParameter(p, 'yMin', default_yMin, @isnumeric)
    addParameter(p, 'yMax', default_yMax, @isnumeric)
    addParameter(p, 'boost', default_boost, @isnumeric)
    addParameter(p, 'savePlotFile', default_savePlotFile, @isstring)
    parse(p,varargin{:})

    currentTimeData = p.Results.currentTimeData;
    if isempty(currentTimeData)
        currentTimeData = linspace(0, length(currentData)-1, length(currentData))/currentFreq;
    end
    yMin = p.Results.yMin;
    yMax = p.Results.yMax;
    listBoxAlgorithm = p.Results.algorithm;
    switchUnitType = p.Results.unitType;
    waveParameter = p.Results.waveParameter;
    switchYScale  = p.Results.yScale;
    useMs = p.Results.useMs;
    filterAC = p.Results.filterAC;
    removeCOI = p.Results.removeCOI;
    switchSqHz = p.Results.sqHz;
    editFieldSensitivity = p.Results.sensitivity;
    boost = p.Results.boost;
    savePlotFile = p.Results.savePlotFile;



    passBands = [0.3, currentFreq/2-1]; %automatic pass-through frequencies from sampling rate

%% plots raw waveform graph
    rsstFigure = figure;
    wavesPlot = subplot(2, 2, 2);

    % filtering
    y = currentData;
    y = y - mean(y);
    y = double(y);
    filterCoefficients = {};
    [filterCoefficients{1}, filterCoefficients{2}] = butter(3, [passBands(1), passBands(2)]/(0.5*currentFreq));
    y2 =filtfilt(filterCoefficients{1},filterCoefficients{2},y);
    
    %filters out 60Hz AC current noise
    if strcmpi(filterAC,'On')
        desFilt = designfilt('bandstopiir','FilterOrder',2, 'HalfPowerFrequency1',55,'HalfPowerFrequency2',65, 'DesignMethod','butter','SampleRate',currentFreq);
        y2 = filtfilt(desFilt,y2);
    end

    wavesY = y;
    wavesX = currentTimeData;

    plot(wavesPlot, wavesX, wavesY, '-k');

    ylim([-max(abs(y)), max(abs(y))])

%% plot SSWT graph
    rsstPlot = subplot(2, 2, 4);

    % sswt calculations
    dataRealLength = length(y2);

    % zeros added to remove coi
    if strcmpi(removeCOI, 'On')
        disp('Imprecision from COI still exists, removal is only for aesthetics')
        zerosLength = ceil(dataRealLength/3);
        addedFront = zeros(1,zerosLength).';% + mean(currentData(1:10));
        addedEnd = zeros(1,zerosLength).';%+  mean(currentData(dataRealLength-10:dataRealLength));
        y3 = cat(1,addedFront,y2);
        y2 = cat(1,y3,addedEnd);
    end

    if strcmpi(listBoxAlgorithm,'SSWT')
        [cfs,frequencies] = rsst(y2,waveParameter,currentFreq);
    elseif strcmpi(listBoxAlgorithm,'CWT')
        % the legacy method allows for changing the wave parameter of CWT
        fb = round(2*((waveParameter)/(2*pi*1.5))^2, 2);
        scales = logspace(0, log10(currentFreq), 100);
        wavestr = ['cmor' num2str(fb) '-1.5'];
        [cfs,frequencies] = wavelet.internal.cwt(y2,scales,wavestr,1/currentFreq);
        %[cfs,frequencies] = cwt(y2,currentFreq);
    end

    if strcmpi(removeCOI, 'On')
        cfs = cfs(:,1+zerosLength:(dataRealLength)+zerosLength);
    end
    c=real(cfs);
    if strcmpi(switchUnitType, 'Power')
        c = abs(cfs).^2;
    end

    % increase sensitivity for higher frequency data
    if strcmpi(switchSqHz,'On')
        for i = 1:height(c)
            % multiplying factor ln(f+3)
            %m = log(frequencies(i)+3);
            % multiplying factor = root_2(f)
            m = nthroot(frequencies(i),2);

            c(i,:) = c(i,:)*m;

        end
    end


    
    sensitivityMag = 100.^editFieldSensitivity;
    % maximum magnitude of c
    Vmax=max(max(abs(real(c))));
    % c1 is each value of c as a proportion of the max magnitude
    % changes color to a log scale
    c1=c/Vmax;
    negIndex=find(c1<-1/sensitivityMag);
    posIndex=find(c1> 1/sensitivityMag);
    zerIndex=find(-1/sensitivityMag <= c1 & c1 <= 1/sensitivityMag);
    c1(negIndex)=-log10(sensitivityMag*abs(c1(negIndex)));
    c1(posIndex)=log10(sensitivityMag*abs(c1(posIndex)));
    c1(zerIndex)=0;
    analysisData = c1/(2*editFieldSensitivity);

    % records data
    analysisTimeData = currentTimeData;
    analysisFrequencies = frequencies;

    % set color scale
    pc = pcolor(rsstPlot, analysisTimeData,analysisFrequencies,analysisData);%*max(abs(currentData))
    grid on; shading interp;
    set(pc, 'EdgeColor', 'none');
    % set value
    Amax = max(analysisData, [], 'all');
    Amin = min(analysisData, [], 'all');
    if strcmp(switchUnitType, "Potential")
        ncol=2048;
        deleteCols = 2;
        shiftCol = -round(0.5*ncol*(Amax+Amin)/(abs(Amax)+abs(Amin)));
        ColorInit=jet(ncol);
        jetEdit=[ColorInit([1:ncol/2-deleteCols+shiftCol],:);repmat([1 1 1], [10,1]);ColorInit([ncol/2+deleteCols+shiftCol:ncol],:)];
        colormap(rsstPlot, jetEdit);
    else
        ncol=2048;
        deleteCols = 2;
        shiftCol = 0;
        ColorInit=jet(ncol);
        jetEdit=[ColorInit([1:ncol/2-deleteCols+shiftCol],:);ColorInit([ncol/2+deleteCols+shiftCol:ncol],:)];
        colormap(rsstPlot, jetEdit);
        caxis([0,1]);
    end
    Cmax = max(caxis);
    Cmin = min(caxis);

    % labels colorbar
    cb = colorbar();
    if strcmp(switchUnitType, "Potential")
        cb.TickLabels = round([-1, -0.3, -0.1, -0.03, 0, 0.03, 0.1, 0.3, 1],3);
        factor = 2*editFieldSensitivity;
        cb.Ticks = [Amin, Amin+log10(3)/(factor), Amin+log10(10)/(factor), Amin+log10(30)/(factor), 0, Amax-log10(30)/(factor), Amax-log10(10)/(factor),Amax-log10(3)/(factor), Amax];
    end

%% plots frequency x magnitude
    freqPlot = subplot(2, 2, [1,3]);
    
    freqData = mean(abs(c'));
    i = find(frequencies>yMin);
    j = find(frequencies<yMax);
    valid_freqData = freqData(i(1):j(end));
    fMax = max(valid_freqData,[],'all');
    freqData = freqData/fMax;

    L = length(y2);
    n = (2^nextpow2(L))*16;  % change 16 to higher power of 2 for more frequencies in FFT
    Y = fft(y2, n);
    P2 = abs(Y/n);
    P1 = P2(1:n/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = currentFreq*(0:(n/2))/n;       
    i = find(f>yMin);
    j = find(f<yMax);
    valid_P1 = P1(i(1):j(end));
    fMax_ = max(valid_P1,[],'all');
    P1 = P1/fMax_;

    % smooth from roughly the middle of frequency ticks   
    step_value = yMin;
    steps = nextpow2(yMax/step_value)-1;
    valid_length = j(end) - i(1);
    step_smooth = valid_length/(steps*20);  % Change 20 for different smoothing degree. smaller=smoother
    P = sgolayfilt(P1, steps, 101);
    for i = 1:steps
        smooth_start_value = step_value*2;
        smooth_start_index = round(n/currentFreq*step_value + n/currentFreq*(smooth_start_value-step_value)/2); 
        P = cat(1, P(1:smooth_start_index - 1), sgolayfilt(P1(smooth_start_index:end), steps-i, 101+2*floor((step_smooth*i)/2)));
        step_value = smooth_start_value;
    end

%     disp(class(P))
%     disp(size(P))
    
    % Fit FFT
    P_ = P(f>5 & f<yMax);
    f_ = f(f>5 & f<yMax);
    fit_func = fit(f_',P_,'power1');
    fit_cf = coeffvalues(fit_func);
%     disp(fit_func)
    y_ = f_';
    x_ = fit_cf(1)*(f_'.^fit_cf(2));
    
    % Fit RSST
    freqData_ = freqData(frequencies>5 & frequencies<yMax);
    frequencies_ = frequencies(frequencies>5 & frequencies<yMax);
    fit_func2 = fit(frequencies_',freqData_','power1');
    fit_cf2 = coeffvalues(fit_func2);
%     disp(fit_func2)
    y2_ = frequencies_';
    x2_ = fit_cf2(1)*(frequencies_'.^fit_cf2(2));
    
    ft = loglog(P, f, freqData, frequencies, x_, y_, x2_, y2_); % blue is FFT, black is rsst
    

    ft(1).Color = [0.15 0.5 0.75];
    ft(2).Color =  [0 0 0];
    ft(2).LineWidth = 1.2;
    ft(3).Color = [0.85 0.325 0.1];
    ft(3).LineWidth = 1.5;
    ft(3).LineStyle = '--';
    ft(4).Color = [0.5 0.18 0.56];
    ft(4).LineWidth = 1.5;
    ft(4).LineStyle = '--';
    ylim([yMin yMax])
    xlim([0.01, 1.5])
    xticks([0.01 0.1 1])
    xticklabels([0.01 0.1 1])
    
     legend(ft,{'FFT','RSST',sprintf('FFT \\alpha= %0.2f',-fit_cf(2)),sprintf('RSST \\alpha= %0.2f',-fit_cf2(2))},'Position',[0.1 0.75 0.06 0.04 ],'FontSize',8)


%     ft = loglog(P, f, freqData, frequencies); % blue is FFT, black is rsst
% 
%     ft(1).Color = [0.15 0.5 0.75];
%     ft(2).Color =  [0 0 0];
%     ft(2).LineWidth = 1.2;
%     ylim(freqPlot, [yMin yMax])
%     xlim(freqPlot, [0.01, 1.5])
%     xticks(freqPlot, [0.01 0.1 1])
%     xticklabels(freqPlot, [0.01 0.1 1])



%% Labeling
    % sets x axis
    xticks(rsstPlot, 'auto');

    % sets y axis
    if strcmpi(switchYScale, 'Linear')
        rsstPlot.YScale = 'linear';
        yticks(rsstPlot,'auto');
    else % log
        rsstPlot.YScale = 'log';
        ytickArray = [0.3, 0.625];
        while ytickArray(end) < yMax/2
            ytickArray(end+1) = ytickArray(end)*2;
        end
        freqPlot.YTick = ytickArray;
    end
    set(rsstPlot, 'ytick', [])

    % boost sensitivty of the data
    caxis(rsstPlot, [Cmin/boost, Cmax/boost]);



%% Formatting changes
    set(wavesPlot, 'XAxisLocation', 'top');
    xticks(wavesPlot, [0, 1]);

    wavesPlot.YAxis.Limits = wavesPlot.YAxis.Limits * 1.05;
    xlim(wavesPlot, xlim(rsstPlot))
    ylabel(wavesPlot, 'Potential (\muV)');
    ylabel(freqPlot, 'Frequency (Hz)');
    if strcmpi(useMs, 'On')
        xLabels = rsstPlot.XTick * 1000;
        rsstPlot.XTickLabel = xLabels;
        xlabel(rsstPlot, 'Time (ms)');
    else
        xlabel(rsstPlot, 'Time (s)');
    end
    ylim(rsstPlot, [yMin, yMax])

    set(rsstFigure,'color','w');
    box(rsstPlot,'on');
    box(wavesPlot,'on');
    wavesPlot.XTick = [];
    rsstPlot.XAxis.TickLength = [0 0];
    rsstPlot.YAxis.TickLength = [0 0];
    rsstPlot.YRuler.Axle.LineWidth = 1.01;
    grid(rsstPlot,"on");

    title(wavesPlot,plotTitle, 'FontSize', 20);
    cb.FontSize = 15;
    wavesPlot.FontSize = 15;
    rsstPlot.FontSize = 15;
    freqPlot.FontSize = 15;


    plotDescription = "F/SDf(waveParameter)=" + num2str(waveParameter) +...
        "," + "SamplingRate="+ num2str(currentFreq) +...
        ",Time=" + num2str(wavesX(end)) +...
        ",Yscale:"+ convertCharsToStrings(switchYScale) +...
        "," + convertCharsToStrings(listBoxAlgorithm) +...
        "," + convertCharsToStrings(switchUnitType) +...
        ",sqrt(Hz)" + convertCharsToStrings(switchSqHz) +...
        ",filterAC:" + convertCharsToStrings(filterAC) +...
        ",removeCOI:" + convertCharsToStrings(removeCOI) +...
        ",sensitivity:" + num2str(editFieldSensitivity) +...
        ",boost:" + num2str(boost) +...
        ",Version:" + num2str(VersionNum);
    annotation('textbox', [0, 1.0, 0, 0], 'string', plotDescription);

    cb.Position = [0.9 0.1 0.01 0.6];
    rsstPlot.InnerPosition = [0.25 0.1 0.65 0.6];
    wavesPlot.InnerPosition = [0.25 0.701 0.65 0.2];
    rsstFigure.Position = [100 100 1200 700];
    freqPlot.InnerPosition = [0.1 0.1 0.15 0.6];



    if ~strcmp(savePlotFile,'')
        savePlotFile = extractBefore(savePlotFile,'.');
        savePlotFile = append(savePlotFile, '.png');
        saveas(rsstFigure, savePlotFile);
        disp('figure saved as');
        disp(savePlotFile)
    end
    disp("finished plotting")

end