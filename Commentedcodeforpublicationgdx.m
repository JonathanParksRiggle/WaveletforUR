%% Instructions for User
%Jlab from (http://www.jmlilly.net/jmlsoft.html) and 
%CalcScaleForCWT.m from Tanya Leise need to be added to the path of matlab running the code for the code to fully function. 
%Attempts have been made to be as clear as concise in the annotation as
%possible but if there is any questions feel free to reach out and we will happily attempt to clear up any confusion. 

%% Parameter set up (adapted from wavelet_circ_leise.m)
load('\Users\jrigg\Box\Wavelet\JP wavelets\FourCore Project\gdx_struct.mat') %loads in the data from where it is stored
chnames = ["C001" "C002" "C004" "C005" "C008" "C009" "C010" "C011" "C012" "C013" "C014" "C015" "C016" "C017" "C018" "C019" "C020" "C021" "C022" "C023" "C024" "C025" "C026" "C027" "C028" "C057" "C058" "C059" "C061" "C062" "C064" "C065" "C066" "C067" "C068" "C069" "C070" "C072" "C073" "C074" "C075" "C076" "C077" "C078" "C079" "C080" "C081" "C082" "C083" "C084"]; % channel numbers (refering to each of our PIR sensors) that will be used 
Longs = [30 6.5 4.26 2.13 1.07]; 
Shorts = [.5 .5 2.13 1.07 .53];
wname = ["CRFull", "URFull", "Long", "Medium", "Short"];
nvoices = [34 27.5 100 101 99]; %give tau sizes of ~100 (except the largest
for lv = 1:length(chnames)
    ch = chnames(1,lv);
    NstepsPerHr=60; % sampling rate in minutes 
    Nsteps=length(gdx_struct.(ch).WheelTurns);  % total number of samples
    t=(0:Nsteps-1)'/NstepsPerHr; % creates timept vector
    T=t(Nsteps)-t(1); % time at last time point [reqiuired for CalcScaleForCWT
    for window = 1:length(Longs) 
        longestperiod=Longs(window);shortestperiod=Shorts(window);  % range of periods to use in AWT in hours 
        gamma=3;beta=10; % parameter values for Morse wavelet function 
        nvoice = nvoices(window); % adjust to get 100 tau values depending on window size 
        [fs,tau, qscaleArray] = CalcScaleForCWT(shortestperiod,longestperiod,T,NstepsPerHr,nvoice);%approximates the periods with each wavelet scale
        [row,~,~] = find(tau == tau(tau >= shortestperiod & tau <= longestperiod)'); %allows you to set
        %the true window you want insures that you are looking just at shortest to
        %longest period in your approximate was used for the windows for wavelet
        %ridge 
        vars{lv,1}{window, 1} = qscaleArray(row:row(end,1),1); %adjust each array appropriately
        vars{lv,1}{window, 2} = fs(row:row(end,1),1); %adjust each array appropriately
        vars{lv,1}{window, 3} = tau(row:row(end,1),1); %adjust each array appropriately
        vars{lv,1}{window, 4} = longestperiod;
        vars{lv,1}{window, 5} = shortestperiod;
        vars{lv,1}{window, 6} = t;
    end
end
File1 = ('\Users\jrigg\Box\Wavelet\JP wavelets\FourCore Project\gdxsexandsurgerydecoder.xlsx');
decoder = xlsread(File1);% this allows us to decode the sex and genotype of the animals associated with each channel.first column is surgery(sham = 0 gdx =1) second column is sex (female=0 male =1)
%% Dealing with Outliers and NaN
numchs = length(chnames);
for lv = 1:numchs %runs through channel and activity count within and if it is an outlier or an NaN replaces it with an interpolated value from the average of the count preceding and following it
    ch = chnames(1,lv);
    WheelTurns.(ch) = gdx_struct.(ch).WheelTurns;
    WheelTurns.(ch) = fillmissing(WheelTurns.(ch), 'movmean',10);
    [row,~] = find(WheelTurns.(ch) > 4* std(WheelTurns.(ch)));
    WheelTurns.(ch)(row,1) = 4* std(WheelTurns.(ch));
    gdx_struct.(ch).WheelTurns = WheelTurns.(ch);
end

%% Wavelet, Normalization, and Ridge Analyis
lites = ["D", "L", "F"];
for lv =1:numchs 
    ch = chnames(1,lv);
    for window = 1:length(Longs) 
        fulltotal = length( vars{lv,1}{window, 2}) * (Nsteps-1); %number of cells in the full wavelet matrix
        halftotal = fulltotal/2; %number of cells in either the light or dark phase wavelet matrix
        total = {halftotal, halftotal, fulltotal}; 
        edgecut = (vars{lv,1}{window, 4} * 1.5) * 60; %this is the length of the edge that needs to be cut off on either side in minutes. 
        light_idx = mod(gdx_struct.(ch).DayPhase,2) == 1; %seperates out the lightphase vs. the darkphase 
        dark_idx = mod(gdx_struct.(ch).DayPhase,2) == 0; 
        WheelTurns_D = gdx_struct.(ch).WheelTurns(dark_idx);
        WheelTurns_L = gdx_struct.(ch).WheelTurns(light_idx);
        WheelTurns_F = gdx_struct.(ch).WheelTurns;
        Wheel = {WheelTurns_D, WheelTurns_L, WheelTurns_F};
        for lite = 1:length(lites)
            lit = lites(1,lite);
            totaler = total{1,lite};
            wheeler = Wheel{1,lite};
            fs = vars{lv,1}{window, 2};
            gdx_struct.(ch).cwt.(lit){window,1} = wavetrans(wheeler,{1,gamma,beta,fs,'bandpass'},'periodic'); %runs the wavelet with periodic boundary conditions on each of the three phases light, dark or full
            gdx_struct.(ch).norm.(lit){window,1} =  (abs(gdx_struct.(ch).cwt.(lit){window,1}))./(sum(abs(gdx_struct.(ch).cwt.(lit){window,1}),'all')/totaler); %normalize
            gdx_struct.(ch).norm.(lit){window,1} = gdx_struct.(ch).norm.(lit){window,1}(ceil(edgecut):1:ceil(end-edgecut),:); %trim edges
            [~,~,JR,]= ridgewalk(gdx_struct.(ch).cwt.(lit){window,1}(ceil(edgecut):1:ceil(end-edgecut),:),fs);
            count = 0;
            tau = vars{lv,1}{window, 3};
            for ridge = 1:length(JR)
                if isnan(JR(ridge,1)) == true 
                    continue
                elseif isnan(JR(ridge,1)) == false 
                    count = count + 1;
                    gdx_struct.(ch).ridgeperiod.(lit){window,1}(count ,1) = tau(JR(ridge,1),1);
                end
            end
        end
    end
end 

%% Sorting into Norm CWT, Ridge, and Cells of Avg Power by Scale 
count = 0;
for lite = 1:length(lites)
    lit = lites(1,lite);
    count00 = 0;
    count01 = 0;
    count10 = 0;
    count11 = 0;
    for chs =1:numchs %sort through each channel and sort it into what sex and surgery
        ch = chnames(1,chs);
        if decoder(chs,1) == 1 && decoder(chs,2) == 1%first column is surgery(sham = 0 gdx =1) second column is sex (female=0 male =1)
            count11 = count11 + 1;
            for window = 1:length(Longs) 
                windowname = wname(1,window);
                count = count + 1;
                analysis.(lit).GM.Norm.(ch){window,1} = gdx_struct.(ch).norm.(lit){window,1};
                Ridgeperiod(count,:) =  mean(gdx_struct.(ch).ridgeperiod.(lit){window,1});
                analysis.(lit).GM.PbS.(windowname)(count11,:) = mean(gdx_struct.(ch).norm.(lit){window,1},1);
                analysis.F.GM.timeseries.(ch) = gdx_struct.(ch).WheelTurns;
                lights(count,1) = lit;
                chan(count,1) = ch;
                group(count,1) = "GM";
                windownames(count,1) = windowname; 
            end
        elseif decoder(chs,1) == 1 && decoder(chs,2) == 0
            count10 = count10 + 1;
            analysis.F.GF.timeseries.(ch) = gdx_struct.(ch).WheelTurns;
            for window = 1:length(Longs) 
                windowname = wname(1,window);
                count = count + 1;
                analysis.(lit).GM.Norm.(ch){window,2} = ch;
                analysis.(lit).GF.Norm.(ch){window,1} = gdx_struct.(ch).norm.(lit){window,1};
                Ridgeperiod(count,:) =  mean(gdx_struct.(ch).ridgeperiod.(lit){window,1});
                analysis.(lit).GF.PbS.(windowname)(count10,:) = mean(gdx_struct.(ch).norm.(lit){window,1},1);
                lights(count,1) = lit;
                chan(count,1) = ch;
                group(count,1) = "GF";
                windownames(count,1) = windowname; 
            end
        elseif decoder(chs,1) == 0 && decoder (chs,2) == 1
            count01 = count01 + 1;
            analysis.F.SM.timeseries.(ch) = gdx_struct.(ch).WheelTurns;
            for window = 1:length(Longs) 
                windowname = wname(1,window);
                count = count + 1;
                analysis.(lit).GM.Norm.(ch){window,2} = ch;
                analysis.(lit).SM.Norm.(ch){window,1} = gdx_struct.(ch).norm.(lit){window,1};
                Ridgeperiod(count,:) =  mean(gdx_struct.(ch).ridgeperiod.(lit){window,1});
                analysis.(lit).SM.PbS.(windowname)(count01,:) = mean(gdx_struct.(ch).norm.(lit){window,1},1);
                lights(count,1) = lit;
                chan(count,1) = ch;
                group(count,1) = "SM";
                windownames(count,1) = windowname;
            end
        elseif decoder(chs,1) == 0 && decoder(chs,2) == 0
            count00 = count00 + 1;
            analysis.F.SF.timeseries.(ch) = gdx_struct.(ch).WheelTurns;
            for window = 1:length(Longs) 
                windowname = wname(1,window);
                count = count + 1;
                analysis.(lit).GM.Norm.(ch){window,2} = ch;
                analysis.(lit).SF.Norm.(ch){window,1} = gdx_struct.(ch).norm.(lit){window,1};
                Ridgeperiod(count,:) =  mean(gdx_struct.(ch).ridgeperiod.(lit){window,1});
                analysis.(lit).SF.PbS.(windowname)(count00,:) = mean(gdx_struct.(ch).norm.(lit){window,1},1);
                lights(count,1) = lit;
                chan(count,1) = ch;
                group(count,1) = "SF";
                windownames(count,1) = windowname; 
            end
        elseif decoder(chs,1) == 2 %exclude animals whose surgical status could not be verified post mortem
           continue 
        end
    end
end


%% Write Ridge values
ridge_table = table(Ridgeperiod,lights,group,windownames,chan);%makes a table writes these values to a csv
location = 'C:\Users\jrigg\Box\Data\Wavelet2021\ridgetablegdx.csv';
writetable(ridge_table, location);


%% 95% Confidence Intervals
 grpfields = fieldnames(analysis.D)';
for grper = 1:length(grpfields)
    grp = grpfields{1,grper};
    for lite = 1:length(lites)
        lit = lites(1,lite);
        for window = 1:length(Longs) 
            windowname = wname(1,window);
            [bootCIs,bootstat] = bootci(2000,@mean,analysis.(lit).(grp).PbS.(windowname)); %boothstrap the mean period power curve to generate 95% confidence intervals 
            analysis.(lit).(grp).CIs{window,1}  = bootCIs; %store the confidence intervals 
            analysis.(lit).(grp).PC{window,1} = bootstat; % store the boot strapped mean the confidence intervals surroud 
        end
    end
end
%% Plot Group Curves
for grper = 1:length(grpfields)
    grp = grpfields{1,grper};
    for lite = 1:length(lites)
        lit = lites(1,lite);
        for window = 2 %Select just UR FULL
            tau = vars{1,1}{window, 3};
            windowname = convertStringsToChars(wname(1,window));
            figure
            plot(tau, analysis.(lit).(grp).CIs{window,1}','k')
            hold on
            plot(tau, nanmean(analysis.(lit).(grp).PC{window,1},1),'k--')
            title(['The ', windowname, ' ', convertStringsToChars(grp), ' ', convertStringsToChars(lit), ' Powercurve'])
            axis([0 7 0 2])
            locationfig = 'C:\Users\jrigg\Box\Data\Wavelet2021\';
            fname = ([windowname,convertStringsToChars(grp),convertStringsToChars(lit),'powerperiodcurve','.eps']);
            saveas(gcf, fullfile(locationfig, fname), 'epsc')
        end
    end
end
    
%% Lomb-Scargle Periodogram
counter = 1;
for grper = 1:length(grpfields)
    grp = grpfields{1,grper};
    chnamesbygrp = fieldnames(analysis.F.(grp).timeseries);
    for chs =1:length(chnamesbygrp)
        ch = chnamesbygrp{chs,1};
        x = analysis.F.(grp).timeseries.(ch);%specifies the time series being used 
        [pxx,f]= plomb(x,t); % run the lombs scargle 
        cirpwr = max(pxx(37:44,1)); % look for the maximum in the circadian range 
        [row,~,~] = find(pxx == cirpwr); %store variables of interest to export
        cirper = 1/f(row,1);
        maxpwr = max(pxx);
        clear row
        [row,~,~] = find(pxx == maxpwr); 
        maxper = 1/f(row,1);
        lsp_table(counter,1:6) = table(convertCharsToStrings(ch), convertCharsToStrings(grp), cirpwr, cirper, maxpwr, maxper); %store as a table
        counter = counter + 1; %to properly store this
    end
end

%% Modified Discrete Wavelet Transform and Energy Calculation (code adapted from matlab documention of Modwt function)
decoder2 = {'2-4 min'; '4-8 min'; '8-16 min'; '16-32 min'; '32-64 min'; '64-128 min'; '128-256 min'; '256-512 min'; '512-1024 min'; '1024-2048 min'; '1024-2048 min';};
Levels = {'D1';'D2';'D3';'D4';'D5';'D6';'D7';'D8';'D9';'D10';'A10'};
counter = 1;
edgecutter = 2048 *1.5; %maximum period accessed dealing with boundary effects
for grper = 1:length(grpfields)
    grp = grpfields{1,grper};
    chnamesbygrp = fieldnames(analysis.F.(grp).timeseries);
    for chs =1:length(chnamesbygrp)
        ch = chnamesbygrp{chs,1};
        premdwt = smoothdata(analysis.F.(grp).timeseries.(ch), 'movmean', 30);%takes the 30 minute moving average
        mdwt = modwt(premdwt,'sym6',length(Levels));%performs the modified discrete wavelet with periodic boundary extension (default)
        mdwt = mdwt(:,ceil(edgecutter):1:ceil(end-edgecutter));%trims for edge effects 
        energy_by_scales1 = sum(mdwt.^2,2);%calculates the energy explained by each scale
        energy_by_scales1 = energy_by_scales1/sum(energy_by_scales1, 'all');%normalizes it
        energy_by_scales2(counter,1:length(Levels)) = (energy_by_scales1(1:length(Levels))/sum(energy_by_scales1(1:length(Levels))));
        channels(counter,1) = convertCharsToStrings(ch);
        groups(counter,1) = convertCharsToStrings(grp);
        levelers(counter,1:length(Levels)) = Levels;
        decoders(counter,1:length(decoder2)) = decoder2;
        counter = counter + 1;
    end
end

energy_table = table(groups, levelers, channels, decoders, energy_by_scales2);
%% Save the Energy and LSP Data
location = 'C:\Users\jrigg\Box\Data\Wavelet2021\energytablegdx.csv';
location2 = 'C:\Users\jrigg\Box\Data\Wavelet2021\lsptablegdx.csv';
writetable(energy_table, location);
writetable(lsp_table,location2);

