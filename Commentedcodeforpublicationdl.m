%%% Instructions for User
%Jlab from (http://www.jmlilly.net/jmlsoft.html) and 
%CalcScaleForCWT.m from Tanya Leise need to be added to the path of matlab running the code for the code to fully function. 
%Attempts have been made to be as clear as concise in the annotation as
%possible but if there is any questions feel free to reach out and we will happily attempt to clear up any confusion. 
%Please email Jpriggle@uchicago.edu 
%% create M/F logical (0=M)
tic
Sensor = xlsread('\Users\jrigg\Box\Wavelet\JP wavelets\Leslie Wavelet\ListOfSexesForJP.xlsx'); %where ever the channel sex decoder is stored
SensorLog = logical(Sensor(:,2)); %  converts the decoder to a logical
%% ultradian wavelet analysis, all channels: spectrogram adapted from wavelet_circ_leise.m 
load('\Users\jrigg\Box\Wavelet\JP wavelets\Leslie Wavelet\wavelet procedures\ultradDataGamma5beta9.mat')%load the data matrix
chNames = fieldnames(ultrad); %create an array of the field names
numChs = size(chNames,1); % number of channels
DL = {'longDL', 'intDL', 'shortDL'}; %three daylengths
numDL = length(DL); %size
week = {'W06' 'W08' 'W10'}; %week array
Longs = [30 6.5 4.26 2.13 1.07]; % Upper limit of analysis window in hours
Shorts = [.5 .5 2.13 1.07 .53]; % Lower limit
wname = ["CRFull", "URFull", "Long", "Medium", "Short"];
nvoices = [34 27.5 100 101 99]; %give tau sizes of ~100 (except the largest

NstepsPerHr=60; % sampling rate in minutes 
Nsteps=length(ultrad.ch183.longDL.W06.Cntmin);  % total number of samples
t=(0:Nsteps-1)'/NstepsPerHr; % creates timept vector
T=t(Nsteps)-t(1); % time at last time point [reqiuired for CalcScaleForCWT
for window = 1:length(Longs) %iterates through the various analysis windows
    longestperiod=Longs(window);shortestperiod=Shorts(window);  % range of periods to use in AWT in hours 
    gamma=3;beta=10; % parameter values for Morse wavelet function 
    nvoice = nvoices(window); % adjust to get 100 tau values depending on window size 
    [fs,tau, qscaleArray] = CalcScaleForCWT(shortestperiod,longestperiod,T,NstepsPerHr,nvoice);%approximates the periods with each wavelet scale
    [row,~,~] = find(tau == tau(tau >= shortestperiod & tau <= longestperiod)'); %allows you to set the true window you want; 
    %and insures that you are looking just at shortest to longest period. 
    %CalcScaleForCWT just approximates 
    
    vars{window, 1} = qscaleArray(row:row(end,1),1); %adjust each array appropriately so that it usese the exact window specified
    vars{window, 2} = fs(row:row(end,1),1); %adjust each array appropriately
    vars{window, 3} = tau(row:row(end,1),1); %adjust each array appropriately
    vars{window, 4} = longestperiod;
    vars{window, 5} = shortestperiod;
    vars{window, 6} = t;
end
clear shortestperiod longestperiod T NstepsPerHr nvoice nvoices
%% Exclusion of Poor Channels (Exclude channels with missing data or issues)
for channelnumber = 1:numChs %iterates over each channel in the study
    ch = chNames{channelnumber,1}; %stores the channel we are currently examining
    for daylength = 1:(numDL) %iterates over each photoperiod of the channel
        dl = DL{1, daylength};
        if strcmp('longDL',dl) == 1 %if the photoperiod is in long
            for weeknum = 1
                wk = week{1,weeknum};
                if strcmp('ch220',ch) == 1 || strcmp('ch224',ch) == 1 %these channels were broken or missing data
                    rubbish.(dl).(wk).(ch) = ultrad.(ch).(dl).(wk); 
                else 
                    ultradian.(ch).(dl).(wk) = ultrad.(ch).(dl).(wk);
                end   
            end
            
        elseif strcmp('intDL',dl) == 1 %if the photoperiod is intermediate
            for weeknum = 1:2
                wk = week{1,weeknum};
                if strcmp('W06',wk) == 1
                    if strcmp('ch183',ch) == 1 || strcmp('ch199',ch) == 1 || strcmp('ch210',ch) == 1 || strcmp('ch217',ch) == 1 || strcmp('ch224',ch) == 1
                        rubbish.(dl).(wk).(ch) = ultrad.(ch).(dl).(wk);
                    else
                        ultradian.(ch).(dl).(wk) = ultrad.(ch).(dl).(wk);
                    end
                elseif strcmp('W08',wk) == 1
                    if strcmp('ch183',ch) == 1 || strcmp('ch199',ch) == 1 || strcmp('ch210',ch) == 1 || strcmp('ch217',ch) == 1 || strcmp('ch224',ch) == 1 ||  strcmp('ch280',ch) == 1
                        rubbish.(dl).(wk).(ch) = ultrad.(ch).(dl).(wk); 
                    else 
                        ultradian.(ch).(dl).(wk) = ultrad.(ch).(dl).(wk);
                    end
                end
            end
        elseif strcmp('shortDL',dl) == 1 %if the photoperiod is in long
            for weeknum = 1:3  
                wk = week{1,weeknum};
                if strcmp('W06',wk) == 1 %These files contain 11 days instead of 1`0
                    %The computer failed in to collect data for almost 24
                    %hours. Rather than lose the power by running it on 9
                    %days off data, we recorded an additional 11th data.
                    %The following code excludes the missing day and
                    %concatenates.
                 
                    idx1 = isnan(ultrad.(ch).(dl).(wk).Cntmin) == false;
                    ultrad.(ch).(dl).(wk).Cntmin =  ultrad.(ch).(dl).(wk).Cntmin(idx1);
                    ultrad.(ch).(dl).(wk).Cntmin = ultrad.(ch).(dl).(wk).Cntmin(1:14400,:);
                    ultrad.(ch).(dl).(wk).LDphase = ultrad.(ch).(dl).(wk).LDphase(idx1);                
                    ultrad.(ch).(dl).(wk).LDphase = ultrad.(ch).(dl).(wk).LDphase(1:14400,:);
                    ultradian.(ch).(dl).(wk) = ultrad.(ch).(dl).(wk);
                else
                    ultradian.(ch).(dl).(wk) = ultrad.(ch).(dl).(wk);
                end
            end
        end
    end
end 
clear rubbish TimeStep idx1
%% Dealing with Outliers and converting NaN to zeroes 
for channelnumber = 1:numChs %iterate through channels
    ch = chNames{channelnumber,1};
    DL = fieldnames(ultradian.(ch));
    numDL = length(DL);
    for daylength = 1:numDL % iterates over photoperiod/day lengths, 
        dl = DL{daylength,1};
        weeknum = length(fieldnames(ultradian.(ch).(dl)));
        for wks = 1:weeknum % iterates over each week used for analysis
            wk = week{1,wks};
            Cntmin.(ch).(dl).(wk) = ultradian.(ch).(dl).(wk).Cntmin;
            Cntmin.(ch).(dl).(wk) = fillmissing(Cntmin.(ch).(dl).(wk),'movmean',10);
            Cntmin.(ch).(dl).(wk) = fillmissing(Cntmin.(ch).(dl).(wk),'constant',0);
            [row,~] = find(Cntmin.(ch).(dl).(wk) > 4* std(Cntmin.(ch).(dl).(wk)));
            Cntmin.(ch).(dl).(wk)(row,1) = 4* std(Cntmin.(ch).(dl).(wk));
            ultradian.(ch).(dl).(wk).Cntmin = Cntmin.(ch).(dl).(wk);
        end
    end
end
         

clear ultrad row
%% Wavelet/Normalization/Ridge/Sorting/Generating Power curves
lites = ["D", "L", "F"]; %Light phase: either dark phase, light phase, or unparsed by phase
count2 = 0;
malecount = 0;
femalecount = 0;
for channelnumber = 1:numChs %iterate through channels, day lengths, and week
    ch = chNames{channelnumber,1};
    DL = fieldnames(ultradian.(ch));
    numDL = length(DL);
    if Sensor(channelnumber,2) == 0 %identifies males
        sex = "male";
        malecount = malecount + 1; %NEED TO update this, so that if its excluded don't add the count
    elseif Sensor(channelnumber,2) == 1 %identifies females
        sex = "female";
        femalecount = femalecount + 1;
    end
    for daylength = 1:numDL % iterates over photoperiod/day lengths, 
        dl = DL{daylength,1};
        weeknum = length(fieldnames(ultradian.(ch).(dl)));
        for wks = 1:weeknum
            wk = week{1,wks};
            for window = 2:length(Longs)%set to 2 to ignore CRFull window as that takes up too much memory 
                windowname = wname(1,window);
               
                edgecut = (vars{window, 4} * 1.5) * 60; %this is the length of the edge that needs to be cut off on either side in minutes. 
                light_idx = mod(ultradian.(ch).(dl).(wk).LDphase,2) == 1; %seperates out the lightphase vs. the darkphase 
                dark_idx = mod(ultradian.(ch).(dl).(wk).LDphase,2) == 0; 
                Cntmin_D = ultradian.(ch).(dl).(wk).Cntmin(dark_idx);
                Cntmin_L = ultradian.(ch).(dl).(wk).Cntmin(light_idx);
                Cntmin_F = ultradian.(ch).(dl).(wk).Cntmin; 
                fulltotal = length(vars{window, 2}) * (Nsteps); %number of cells in the full wavelet matrix
                halftotal = fulltotal/2; %number of cells in either the light or dark phase wavelet matrix
                total = {halftotal, halftotal, fulltotal}; %THIS IS NOT A MISTAKE WE ARE PUTTING EVERYTHING ON THE SAME SCALE AS INTDL
                Cnt = {Cntmin_D, Cntmin_L, Cntmin_F};
                clearvars light_idx dark_idx
                for lite = 1:length(lites)
                    lit = lites(1,lite);
                    totaler = total{1,lite};
                    cnter = Cnt{1,lite};
                    fs = vars{window, 2};
                    tempcwt = wavetrans(cnter,{1,gamma,beta,fs,'bandpass'},'periodic'); %runs the wavelet with periodic boundary conditions on each of the three phases light, dark or full
                    clearvars cnter
                    tempnorm =  (abs(tempcwt))./(sum(abs(tempcwt),'all')/totaler); %normalize
                    tempnorm = tempnorm(ceil(edgecut):1:ceil(end-edgecut),:); %trim edges
                   
                    count = 0;
                    tau = vars{window, 3};
               
                    for ridge = 1:length(tempnorm)
                        tempridgepwr(ridge,1) = max(tempnorm(ridge,:));
                        [~,col,~] = find(tempnorm(ridge,:) == tempridgepwr(ridge,1));  
                        tempridge.(ch).(dl).(wk).(lit)(ridge,1) = tau(col,1);
                    end
                    %assemble components for ridge table 
                    count2 = count2 + 1; 
                    ridgeperiod(count2,1)= mean(tempridge.(ch).(dl).(wk).(lit));
                   % clearvars tempridge
                    light(count2,1) = convertCharsToStrings(lit);
                    windows(count2,1) = windowname;
                    weeks(count2,1) = convertCharsToStrings(wk);
                    daylengths(count2,1) = convertCharsToStrings(dl);
                    chans(count2,1) = convertCharsToStrings(ch); 
                    sexs(count2,1) = sex;
                    
                    %make composite of the curves
                    if window == 2 && strcmp(sex, "male") == 1 %just want the UR Full window
                       grpmat.male.(dl).(wk).(lit)(malecount,:) =  mean(tempnorm);
                    elseif window == 2 && strcmp(sex, "female") == 1 
                        grpmat.female.(dl).(wk).(lit)(femalecount,:) =  mean(tempnorm);
                    else 
                        continue
                    end
                end
               clearvars Cntmin_F Cntmin_L Cntmin_D tempcwt tempnorm
            end
        end
    end
end
%writes the data to a table and stores that table
ridge_table = table(ridgeperiod,light,windows,weeks,daylengths,chans,sexs); 
location3 = 'C:\Users\jrigg\Box\Data\Wavelet2021\ridgetabledl.csv';
writetable(ridge_table, location3);
clear ridgeperiod light windows weeks daylengths chans sexs
%% Get Rid of zeros in Grpmat
% It is very difficult to keep track of the few exceptions where counts
% were included. So here we will remove them.
 for daylength = 2%1:numDL
    dl = DL{daylength,1};
    weeknum = length(fieldnames(grpmat.male.(dl)));
    for wks = 1:weeknum
        wk = week{1,wks};
        for lite = 1:length(lites)
            lit = lites(1,lite);
            [row,~] = find(grpmat.male.(dl).(wk).(lit) == 0);
            grpmat.male.(dl).(wk).(lit)(row,:) =[];
            [row,~] = find(grpmat.female.(dl).(wk).(lit) == 0); 
            grpmat.female.(dl).(wk).(lit)(row, :) = [];
        end
    end
 end
%% Bootstrap
DL = {'longDL', 'intDL', 'shortDL'}; %three daylengths
numDL = length(DL); %size
weeknumDL = [1, 2, 3];
sextype = {'male', 'female'};
for stype = 1:length(sextype)
    sex = sextype{1,stype};
    for daylength = 1:numDL
        dl = DL{1,daylength};
        weeknum = weeknumDL(1,daylength);
        for wks = 1:weeknum
            wk = week{1,wks};
            for lite = 1:length(lites)
                lit = lites(1,lite);
                [bootCIs,bootstat] = bootci(2000,@nanmean, grpmat.(sex).(dl).(wk).(lit));
                bootCIsanalysis.(sex).(dl).(wk).(lit) = bootCIs;
                bootstatsanalysis.(sex).(dl).(wk).(lit) = bootstat;
            end
        end    
    end
end
                
                
%% Plot and Save power period curves 
for stype = 1:length(sextype) %iterates over the males and females 
    sex = sextype{1,stype};
    for daylength = 1:numDL %iterates over daylength
        dl = DL{1,daylength};
        weeknum = weeknumDL(1,daylength);
        for wks = 1:weeknum %iterates over the analysis week
            wk = week{1,wks};
            for lite = 1:length(lites)
                lit = lites(1,lite);
                tau = vars{2, 3}; %using only window 2
                windowname = convertStringsToChars(wname(1,2));
                figure
                plot(tau, bootCIsanalysis.(sex).(dl).(wk).(lit)','k')
                hold on
                plot(tau, mean(bootstatsanalysis.(sex).(dl).(wk).(lit),1),'k--')
                title(['The ', windowname, ' ', convertStringsToChars(sex), ' ', convertStringsToChars(lit), ' ', convertStringsToChars(dl), ' ', convertStringsToChars(wk),  ' Powercurve'])
                axis([0 7 0 2])
                locationfig = 'C:\Users\jrigg\Box\Data\Wavelet2021\';
                fname = ([windowname, sex, convertStringsToChars(dl), convertStringsToChars(wk), convertStringsToChars(lit),'powerperiodcurve','.eps']);
                saveas(gcf, fullfile(locationfig, fname), 'epsc')
            end
        end
    end
end

%% Lomb-Scargle Periodogram
counter =1;
for channelnumber = 1:numChs %iterate through channels, day lengths, and week
    ch = chNames{channelnumber,1};
    if Sensor(channelnumber,2) == 0
        sex = "male";
    elseif Sensor(channelnumber,2) == 1
        sex = "female";
    end
    DL = fieldnames(ultradian.(ch));
    numDL = length(DL);
    for daylength = 1:numDL
        dl = DL{daylength,1};
        weeknum = length(fieldnames(ultradian.(ch).(dl)));
        for wks = 1:weeknum
            wk = week{1,wks};
            x = ultradian.(ch).(dl).(wk).Cntmin;%specifies the time series being used 
            [pxx,f]= plomb(x,t); % run the lombs scargle 
            cirpwr = max(pxx(37:44,1)); % look for the maximum in the circadian range 
            [row,~,~] = find(pxx == cirpwr); %store variables of interest to export
            cirper = 1/f(row,1);
            maxpwr = max(pxx);
            clear row
            [row,~,~] = find(pxx == maxpwr); 
            maxper = 1/f(row,1);
            lsp_table(counter,1:8) = table(sex, convertCharsToStrings(ch), convertCharsToStrings(wk), convertCharsToStrings(dl), cirpwr, cirper, maxpwr, maxper); %store as a table
            counter = counter + 1; %to properly store this
        end
    end
end
%% Complete the Modified Discrete Wavelet Transform (Energy calculated with code adapted from Matlab Modwt documentation page)
decoder2 = {'2-4 min'; '4-8 min'; '8-16 min'; '16-32 min'; '32-64 min'; '64-128 min'; '128-256 min'; '256-512 min'; '512-1024 min'; '1024-2048 min'; '1024-2048 min';};
Levels = {'D1';'D2';'D3';'D4';'D5';'D6';'D7';'D8';'D9';'D10';'A10'};
counter = 1;
edgecutter = 2048 *1.5; %maximum period accessed for dealing with boundary effects 
for channelnumber = 1:numChs %iterate through channels, day lengths, and week
    ch = chNames{channelnumber,1};
    if Sensor(channelnumber,2) == 0
        sex = "male";
    elseif Sensor(channelnumber,2) == 1
        sex = "female";
    end
    DL = fieldnames(ultradian.(ch));
    numDL = length(DL);
    for daylength = 1:numDL
        dl = DL{daylength,1};
        weeknum = length(fieldnames(ultradian.(ch).(dl)));
        for wks = 1:weeknum
            wk = week{1,wks}; 
            premdwt = smoothdata(ultradian.(ch).(dl).(wk).Cntmin, 'movmean', 30);%takes the 30 minute moving average
            mdwt = modwt(premdwt,'sym6',length(Levels));%performs the modified discrete wavelet with periodic boundary extension (default)
            mdwt = mdwt(:,ceil(edgecutter):1:ceil(end-edgecutter));%trims for edge effects 
            energy_by_scales1 = sum(mdwt.^2,2);%calculates the energy explained by each scale
            energy_by_scales1 = energy_by_scales1/sum(energy_by_scales1, 'all');%normalizes it
            energy_by_scales2(counter,1:length(Levels)) = (energy_by_scales1(1:length(Levels))/sum(energy_by_scales1(1:length(Levels))));
            channels(counter,1) = convertCharsToStrings(ch);
            lengths(counter,1) = convertCharsToStrings(dl);
            sexes(counter,1) = sex;
            weeks2(counter,1) = convertCharsToStrings(wk);
            levelers(counter,1:length(Levels)) = Levels;
            decoders(counter,1:length(decoder2)) = decoder2;
            counter = counter + 1;
        end
    end
end

energy_table = table(weeks2, sexes, lengths, levelers, channels, decoders, energy_by_scales2);
%% Write the LSP, Ridge, and Energy to table
location = 'C:\Users\jrigg\Box\Data\Wavelet2021\energytabledl.csv';
writetable(energy_table, location);
location2 = 'C:\Users\jrigg\Box\Data\Wavelet2021\lsptabledl.csv';
writetable(lsp_table, location2);
toc