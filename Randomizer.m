%% Parameter set up
load('\Users\jrigg\Box\Wavelet\JP wavelets\FourCore Project\gdx_struct.mat') 
chnames = {'C001' 'C002' 'C004' 'C005' 'C008' 'C009' 'C010' 'C011' 'C012' 'C013' 'C014' 'C015' 'C016' 'C017' 'C018' 'C019' 'C020' 'C021' 'C022' 'C023' 'C024' 'C025' 'C026' 'C027' 'C028' 'C057' 'C058' 'C059' 'C061' 'C062' 'C064' 'C065' 'C066' 'C067' 'C068' 'C069' 'C070' 'C072' 'C073' 'C074' 'C075' 'C076' 'C077' 'C078' 'C079' 'C080' 'C081' 'C082' 'C083' 'C084'};
NstepsPerHr=60; % sampling rate
TimeStep=1/NstepsPerHr;
Nsteps=14400;  % total number of samples
t=(0:Nsteps-1)'/NstepsPerHr; % creates timept vector
T=t(Nsteps)-t(1);      % time at last time point [not sure why this is needed]
longestperiod=6.5;shortestperiod=.5;  % range of periods to use in AWT
%6.5 AND .5 FOR FULL 5.5 and 3 for short frequency and 3 and 1 for short
gamma=3;beta=10; % parameter values for Morse wavelet function (Need to be played around with usually best to hold gamma at 3) can use morsewave from jlab to plot the wavelet function itself
nvoice=64;
[fs,tau, qscaleArray] = CalcScaleForCWT(shortestperiod,longestperiod,T,NstepsPerHr,nvoice);
%[row,~,~] = find(tau == tau(tau >= 3.5 & tau <= 4.5)'); %allows you to set the true window you want
%qscaleArray = qscaleArray(row:row(end,1),1);
%fs = fs(row:row(end,1),1);
%tau = tau(row:row(end,1),1); 
numchs = length(chnames);
File1 = '\Users\jrigg\Box\Wavelet\JP wavelets\FourCore Project\gdxsexandsurgerydecoder.xlsx';
decoder = xlsread(File1);%first column is surgery(sham = 0 gdx =1) second column is sex (female=0 male =1)
%% Dealing with Outliers and NaN
NaNnum = 0;
Internum = 0;
for i = 1:numchs
    ch = chnames{1,i};
    WheelTurns.(ch) = gdx_struct.(ch).WheelTurns;
    for cnts = 1:14400
        if isnan(WheelTurns.(ch)(cnts,1)) == 1
            if cnts == 1
                WheelTurns.(ch)(cnts,1) = (WheelTurns.(ch)((cnts + 1),1) + WheelTurns.(ch)((cnts + 2),1))/2;
                NaNnum = NaNnum + 1;
            elseif cnts == 14400
                WheelTurns.(ch)(cnts,1) = (WheelTurns.(ch)((cnts - 1),1) + WheelTurns.(ch)((cnts - 2),1))/2;
                NaNnum = NaNnum + 1;
            else
                WheelTurns.(ch)(cnts,1) = (WheelTurns.(ch)((cnts + 1),1) + WheelTurns.(ch)((cnts -1),1))/2;
                NaNnum = NaNnum + 1;
            end
        elseif WheelTurns.(ch)(cnts,1) > 4 * std(WheelTurns.(ch))
            if cnts == 1
               WheelTurns.(ch)(cnts,1) = (WheelTurns.(ch)((cnts + 1),1) + WheelTurns.(ch)((cnts + 2),1))/2;
               Internum = Internum + 1;
            elseif cnts == 14400
                WheelTurns.(ch)(cnts,1) = (WheelTurns.(ch)((cnts - 1),1) + WheelTurns.(ch)((cnts - 2),1))/2;
                Internum = Internum + 1;
            else
               WheelTurns.(ch)(cnts,1) = (WheelTurns.(ch)((cnts + 1),1) + WheelTurns.(ch)((cnts - 1),1))/2;
               Internum = Internum + 1;
            end
        else
            continue 
        end
    end
   gdx_struct.(ch).WheelTurns = WheelTurns.(ch);
end
%% Randomizer to test the system only. 
startL = 1;
stopL = 720;
startD = 721;
stopD = 1440;
for i = 1:numchs
    ch = chnames{1,i};
    temp = gdx_struct.(ch).WheelTurns;
    day1L = temp(startL:stopL);
    day1L = day1L(randperm(length(day1L)));
    day1D = temp(startD:stopD);
    day1D = day1D(randperm(length(day1D)));
    day2L = temp((startL + 1440):(stopL + 1440));
    day2L = day2L(randperm(length(day2L)));
    day2D = temp((startD + 1440):(stopD + 1440));
    day2D = day2D(randperm(length(day2D)));
    day3L = temp((startL + 1440*2):(stopL + 1440*2));
    day3L = day3L(randperm(length(day3L)));
    day3D = temp((startD + 1440*2):(stopD + 1440*2));
    day3D = day3D(randperm(length(day3D)));
    day4L = temp((startL + 1440*3):(stopL + 1440*3));
    day4L = day4L(randperm(length(day4L)));
    day4D = temp((startD + 1440*3):(stopD + 1440*3));
    day4D = day4D(randperm(length(day4D)));
    day5L = temp((startL + 1440*4):(stopL + 1440*4));
    day5L = day5L(randperm(length(day5L)));
    day5D = temp((startD + 1440*4):(stopD + 1440*4));
    day5D = day5D(randperm(length(day5D)));
    day6L = temp((startL + 1440*5):(stopL + 1440*5));
    day6L = day6L(randperm(length(day6L)));
    day6D = temp((startD + 1440*5):(stopD + 1440*5));
    day6D = day6D(randperm(length(day6D)));
    day7L = temp((startL + 1440*6):(stopL + 1440*6));
    day7L = day7L(randperm(length(day7L)));
    day7D = temp((startD + 1440*6):(stopD + 1440*6));
    day7D = day7D(randperm(length(day7D)));
    day8L = temp((startL + 1440*7):(stopL + 1440*7));
    day8L = day8L(randperm(length(day8L)));
    day8D = temp((startD + 1440*7):(stopD + 1440*7));
    day8D = day8D(randperm(length(day8D)));
    day9L = temp((startL + 1440*8):(stopL + 1440*8));
    day9L = day9L(randperm(length(day9L)));
    day9D = temp((startD + 1440*8):(stopD + 1440*8));
    day9D = day9D(randperm(length(day9D)));
    day10L = temp((startL + 1440*9):(stopL + 1440*9));
    day10L = day10L(randperm(length(day10L)));
    day10D = temp((startD + 1440*9):(stopD + 1440*9));
    day10D = day10D(randperm(length(day10D)));
    gdx_struct.(ch).WheelTurns = [day1L day1D day2L day2D day3L day3D day4L day4D day5L day5D day6L day6D day7L day7D day8L day8D day9L day9D day10L day10D];
    gdx_struct.(ch).WheelTurns =  gdx_struct.(ch).WheelTurns';
    clear  temp day1L day1D day2L day2D day3L day3D day4L day4D day5L day5D day6L day6D day7L day7D day8L day8D day9L day9D day10L day10D
end
%%
edgecut = (longestperiod * 1.5) * 60;
for i =1:numchs 
    ch = chnames{1,i};
    light_idx = mod(gdx_struct.(ch).DayPhase,2) == 1;
    dark_idx = mod(gdx_struct.(ch).DayPhase,2) == 0; 
    WheelTurns_D{i,1} = gdx_struct.(ch).WheelTurns(dark_idx);
    WheelTurns_L{i,1} = gdx_struct.(ch).WheelTurns(light_idx);
    WheelTurns_F{i,1} = gdx_struct.(ch).WheelTurns;
    gdx_struct.(ch).cwt.D = wavetrans(WheelTurns_D{i,1},{1,gamma,beta,fs,'bandpass'},'periodic');
    gdx_struct.(ch).cwt.L = wavetrans(WheelTurns_L{i,1},{1,gamma,beta,fs,'bandpass'},'periodic');
    gdx_struct.(ch).cwt.F = wavetrans(WheelTurns_F{i,1}, {1,gamma,beta,fs,'bandpass'},'periodic');
end 
%% Sorting and Normalization(avg each animal CWT and then divide CWT by it)
fulltotal = length(qscaleArray) * (Nsteps-1);
halftotal = fulltotal/2;
edgecut = (longestperiod * 1.5) * 60;
for i =1:numchs
     ch = chnames{1,i};
     if decoder(i,1) == 1
         if decoder(i,2) == 1
             
             analysis.gdx.male.L.(ch) = gdx_struct.(ch).cwt.L;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.L),'all')/halftotal;
             analysis.norm.GML.(ch) = abs(analysis.gdx.male.L.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.GML.(ch) = analysis.norm.GML.(ch)(edgecut:1:end-edgecut,:);
             
             analysis.gdx.male.D.(ch) = gdx_struct.(ch).cwt.D;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.D),'all')/halftotal;
             analysis.norm.GMD.(ch) = abs(analysis.gdx.male.D.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.GMD.(ch) = analysis.norm.GMD.(ch)(edgecut:1:end-edgecut,:);
             
             analysis.gdx.male.F.(ch) = gdx_struct.(ch).cwt.F;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.F),'all')/fulltotal;
             analysis.norm.GMF.(ch) = abs(analysis.gdx.male.F.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.GMF.(ch) = analysis.norm.GMF.(ch)(edgecut:1:end-edgecut,:);
             analysis.timeseries.GMF.(ch) = gdx_struct.(ch).WheelTurns;
        
         elseif decoder(i,2) == 0
             
             analysis.gdx.female.L.(ch) = gdx_struct.(ch).cwt.L;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.L),'all')/halftotal;
             analysis.norm.GFL.(ch) = abs(analysis.gdx.female.L.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.GFL.(ch) = analysis.norm.GFL.(ch)(edgecut:1:end-edgecut,:);
             
             analysis.gdx.female.D.(ch) = gdx_struct.(ch).cwt.D;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.D),'all')/halftotal;
             analysis.norm.GFD.(ch) = abs(analysis.gdx.female.D.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.GFD.(ch) = analysis.norm.GFD.(ch)(edgecut:1:end-edgecut,:);
             
             analysis.gdx.female.F.(ch) = gdx_struct.(ch).cwt.F;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.F),'all')/fulltotal;
             analysis.norm.GFF.(ch) = abs(analysis.gdx.female.F.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.GFF.(ch) = analysis.norm.GFF.(ch)(edgecut:1:end-edgecut,:);
             analysis.timeseries.GFF.(ch) = gdx_struct.(ch).WheelTurns;
             
         end
     elseif decoder(i,1) == 0
        if decoder (i,2) == 1
             analysis.sham.male.L.(ch) = gdx_struct.(ch).cwt.L;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.L),'all')/halftotal;
             analysis.norm.SML.(ch) = abs(analysis.sham.male.L.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.SML.(ch) = analysis.norm.SML.(ch)(edgecut:1:end-edgecut,:);
             
             analysis.sham.male.D.(ch) = gdx_struct.(ch).cwt.D;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.D),'all')/halftotal;
             analysis.norm.SMD.(ch) = abs(analysis.sham.male.D.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.SMD.(ch) = analysis.norm.SMD.(ch)(edgecut:1:end-edgecut,:);
             
             analysis.sham.male.F.(ch) = gdx_struct.(ch).cwt.F;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.F),'all')/fulltotal;
             analysis.norm.SMF.(ch) = abs(analysis.sham.male.F.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.SMF.(ch) = analysis.norm.SMF.(ch)(edgecut:1:end-edgecut,:); 
             analysis.timeseries.SMF.(ch) = gdx_struct.(ch).WheelTurns;
             
        elseif decoder(i,2) == 0 
             analysis.sham.female.L.(ch) = gdx_struct.(ch).cwt.L;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.L), 'all')/halftotal;
             analysis.norm.SFL.(ch) = abs(analysis.sham.female.L.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.SFL.(ch) = analysis.norm.SFL.(ch)(edgecut:1:end-edgecut,:);
             
             analysis.sham.female.D.(ch) = gdx_struct.(ch).cwt.D;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.D), 'all')/halftotal;
             analysis.norm.SFD.(ch) = abs(analysis.sham.female.D.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.SFD.(ch) =  analysis.norm.SFD.(ch)(edgecut:1:end-edgecut,:);
             
             analysis.sham.female.F.(ch) = gdx_struct.(ch).cwt.F;
             analysis.avgtotalpower.(ch) = sum(abs(gdx_struct.(ch).cwt.F),'all')/fulltotal;
             analysis.norm.SFF.(ch) = abs(analysis.sham.female.F.(ch))/ analysis.avgtotalpower.(ch);
             analysis.norm.SFF.(ch) = analysis.norm.SFF.(ch)(edgecut:1:end-edgecut,:);
             analysis.timeseries.SFF.(ch) = gdx_struct.(ch).WheelTurns;
        
        end
     elseif decoder(i,1) == 2 
         analysis.reject.L.(ch) = gdx_struct.(ch).cwt.L;
         analysis.reject.D.(ch) = gdx_struct.(ch).cwt.D;
         analysis.reject.F.(ch) = gdx_struct.(ch).cwt.F;
     end
end
%% Grouping Powercurves for each Animal 
chnamesGM = fieldnames(analysis.norm.GMF); %channels in each of the four groups, won't matter depending on light
chnamesSM = fieldnames(analysis.norm.SMF);
chnamesGF = fieldnames(analysis.norm.GFF);
chnamesSF = fieldnames(analysis.norm.SFF);
grpnamesGM = {'GMF', 'GMD', 'GML'}; %allows for iteration through the different light groupings
grpnamesSM = {'SMF', 'SMD', 'SML'};
grpnamesGF = {'GFF', 'GFD', 'GFL'};
grpnamesSF = {'SFF', 'SFD', 'SFL'};

    for j = 2:3
        grp = grpnamesGM{1,j};
        for n = 1:size(chnamesGM,1)
            ch = chnamesGM{n,1};
            analysis.powercurves.(grp){n,1} = mean(abs(analysis.norm.(grp).(ch))); 
        end
    end
    for k = 2:3
        grp = grpnamesSM{1,k};
        for o = 1:size(chnamesSM,1)
            ch = chnamesSM{o,1};
            analysis.powercurves.(grp){o,1} = mean(abs(analysis.norm.(grp).(ch)));
        end
    end
    for l = 2:3
        grp = grpnamesGF{1,l};
        for p = 1:size(chnamesGF,1)
            ch = chnamesGF{p,1};
            analysis.powercurves.(grp){p,1} = mean(abs(analysis.norm.(grp).(ch)));
        end
    end
    for m = 2:3
        grp = grpnamesSF{1,m};
        for q = 1:size(chnamesSF,1)
            ch = chnamesSF{q,1};
            analysis.powercurves.(grp){q,1} = mean(abs(analysis.norm.(grp).(ch)));
        end
    end
    


%% Make Matrices of Power Curves and Ridge Frequencies/Powers
time = t(edgecut:1:(end - edgecut),:);
dtime = t(dark_idx);
dtime = dtime(edgecut:1:end - edgecut,:);
ltime = t(light_idx);
ltime = ltime(edgecut:1:end - edgecut,:);
grpfields = {'GFF' 'GFL' 'GFD' 'GMF' 'GML' 'GMD' 'SFF' 'SFL' 'SFD' 'SMF' 'SML' 'SMD'};

for group = 1:length(grpfields)
   grp = grpfields{1,group};
   if strcmp('GFF',grp) ~= 1 && strcmp('GMF',grp) ~= 1 && strcmp('SMF',grp) ~= 1 && strcmp('SFF',grp) ~= 1
    analysis.grpmat.powercurve.(grp) = cell2mat(analysis.powercurves.(grp));
   else
       continue 
   end
end
%% 95% Confidence Intervals
   for i = 1:length(grpfields)
        grp = grpfields{1,i};
       if strcmp('GFF',grp) ~= 1 && strcmp('GMF',grp) ~= 1 && strcmp('SMF',grp) ~= 1 && strcmp('SFF',grp) ~= 1
        [bootCIs,bootstat] = bootci(2000,@nanmean,analysis.grpmat.powercurve.(grp)); 
        analysis.PC.CIs.(grp)  = bootCIs;
        analysis.PC.bootstat.(grp) = bootstat;
       else
        continue
       end
   end

%% Plot Group Curves and Histograms
for i = 1:length(grpfields)
    grp = grpfields{1,i};
    if strcmp('GFF',grp) == 1||strcmp('GMF',grp) == 1||strcmp('SMF',grp) == 1||strcmp('SFF',grp) == 1
        continue 
    elseif strcmp('GFD',grp) == 1 ||strcmp('GMD',grp) == 1 ||strcmp('SMD',grp) == 1||strcmp('SFD',grp) == 1
        tim = dtime;  
        figure
        plot(tau, analysis.PC.CIs.(grp)','k')
        hold on
        plot(tau, nanmean(analysis.PC.bootstat.(grp),1),'k--')
        title([grp,' Powercurve'])
        axis([0 7 0.5 2])
    elseif strcmp('GFL',grp) == 1||strcmp('GML',grp) == 1||strcmp('SML',grp) == 1||strcmp('SFL',grp) == 1
        tim = ltime;
        figure
        plot(tau, analysis.PC.CIs.(grp)','k')
        hold on
        plot(tau, nanmean(analysis.PC.bootstat.(grp),1),'k--')
        title([grp,' Powercurve'])
        axis([0 7 0.5 2])
    end
end