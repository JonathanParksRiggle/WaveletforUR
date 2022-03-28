locationa = uigetdir; % Identify where to search for files
a = dir([locationa, '/*.csv']); %Store the name of all .xls files as a vector
filename = {a(:).name}.'; %extract the intact file names
files = cell(length(a),1);
for ii = length(a):-1:1 
      % Create the full file name and partial filename
      fullname = [locationa filesep a(ii).name];
      temp = strsplit(a(ii).name,'_'); 
      temp = strcat('Channel', temp); 
      channels{ii,1} = temp{1};
      % Read in the data
      files{ii} =   csvread(fullname, 4, 0);    
end
datafiles = cell2struct(files, channels); 
clear files filename temp fullname
per2_exemplars = datafiles;
clear datafiles channels
%% Set up Variables
chNames = fieldnames(per2_exemplars); %create an array of the field names
numChs = size(chNames,1); % number of channels
Longs = [30 30 6.5 4.26 2.13 1.07]; 
Shorts = [.5 20 .5 2.13 1.07 .53];
wname = ["Full", "CRFull", "URFull", "Long", "Medium", "Short" ];
nvoices = [85 172 27.5 100 101 99]; %give tau sizes of ~100 except thet first window which is ~500 for better resolution
NstepsPerHr=60; % sampling rate in minutes 
Nsteps=length(per2_exemplars.(chNames{1,1})(:,1));  % total number of samples
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
    vars{window, 1} = qscaleArray(row:row(end,1),1); %adjust each array appropriately
    vars{window, 2} = fs(row:row(end,1),1); %adjust each array appropriately
    vars{window, 3} = tau(row:row(end,1),1); %adjust each array appropriately
    vars{window, 4} = longestperiod;
    vars{window, 5} = shortestperiod;
    vars{window, 6} = t;
end
clear shortestperiod longestperiod T NstepsPerHr nvoice nvoices

%% Data Cleaning
for chs = 1:numChs
    ch = chNames{chs, 1};
    Cntmin.(ch) = per2_exemplars.(ch)(:,4);
    Cntmin.(ch) = fillmissing(Cntmin.(ch),'movmean',10);
    Cntmin.(ch) = fillmissing(Cntmin.(ch),'constant',0);
    [row,~] = find(Cntmin.(ch) > 4* std(Cntmin.(ch)));
    Cntmin.(ch)(row,1) = 4* std(Cntmin.(ch));
    per2_exemplars.(ch)(:,4)= Cntmin.(ch);
end
%% Wavelet, Normalization
for chs = 1:numChs
    ch = chNames{chs, 1};
    for window = 1
        windowname = wname(1, window);
        fulltotal = length(vars{window, 2}) * (Nsteps);
        edgecut = (vars{window, 4} * 1.5) * 60;
        fs = vars{window, 2};
        cnts = per2_exemplars.(ch)(:,4);
        analysis.(ch).(windowname).fullcwt = wavetrans(cnts,{1,gamma,beta,fs,'bandpass'},'periodic');
        
    end
end
%% PLOT SCALOGRAMS
for chs = 1:numChs
    ch = chNames{chs, 1};
    for window = 1
        windowname = wname(1, window);
        figure ('visible', 'off')
        h = pcolor((1:length(analysis.(ch).(windowname).fullcwt))/1440, vars{window,3}, abs(analysis.(ch).(windowname).fullcwt)');
        colormap jet
        set(h, 'EdgeColor', 'none')
        set(gca,'TickDir','out');
        filename = strcat('C:\Users\jrigg\Box\Prendergast Lab Box\Wavelet Summer 2021\Figure 3 Data\', ch, 'and', windowname, 'jet');
        saveas(gcf, filename , 'tiffn')
        close gcf
    end
end

%%  SELECT REGIONS OF INTERESTS
% Starts on 9/30 1201
%10/20-10/29 DD 1 1221-1230
%11/19-11/28 ARR 1251-1260
%12/12-12/21 2:2 1274-1283
%1/16 -1/25 DD 1309-1318
chNames = fieldnames(per2_exemplars); %create an array of the field names
numChs = size(chNames,1); % number of channels
Longs = [30 30 6.5 4.26 2.13 1.07]; 
Shorts = [.5 20 .5 2.13 1.07 .53];
wname = ["Full", "CRFull", "URFull", "Long", "Medium", "Short" ];
nvoices = [85 172 27.5 100 101 99]; %give tau sizes of ~100 except thet first window which is ~500 for better resolution
NstepsPerHr=60; % sampling rate in minutes 
Nsteps=14401;  % total number of samples
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
    vars{window, 1} = qscaleArray(row:row(end,1),1); %adjust each array appropriately
    vars{window, 2} = fs(row:row(end,1),1); %adjust each array appropriately
    vars{window, 3} = tau(row:row(end,1),1); %adjust each array appropriately
    vars{window, 4} = longestperiod;
    vars{window, 5} = shortestperiod;
    vars{window, 6} = t;
end
clear shortestperiod longestperiod T NstepsPerHr nvoice nvoices



count2= 0;
intlong = [1230 1260 1283 1318];
intshort = [1221 1251 1274 1309];
intname = ["DD1", "ARR", "TwoTwo", "DD2"];
for interval = 1:length(intlong)
    intL = intlong(interval);
    intS = intshort(interval);
    intN = intname(interval);
    [rowL,~] = find(per2_exemplars.(ch)(:,1) == intL+1);
    [rowS,~] = find(per2_exemplars.(ch)(:,1) == intS);
    for chs = 1:numChs
        ch = chNames{chs, 1};
        if strcmp(ch,'ChannelBF') || strcmp(ch,'ChannelPP') == 1 || strcmp(ch,'ChannelLL') == 1 || strcmp(ch,'ChannelDC') == 1 || strcmp(ch,'ChannelDL') == 1
            geno = "Per2";
        elseif strcmp(ch,'ChannelRU') == 1 || strcmp(ch,'ChannelRV') == 1 || strcmp(ch,'ChannelUE') == 1 || strcmp(ch,'ChannelUX') == 1 || strcmp(ch,'ChannelBQ') == 1 || strcmp(ch,'ChannelXA') == 1 || strcmp(ch,'ChannelQV') == 1
            geno = "WT";
        end           
        for window = 2:length(Longs) 
            windowname = wname(1, window);
            fulltotal = length(vars{window, 2}) * (Nsteps);
            edgecut = (vars{window, 4} * 1.5) * 60;
            fs = vars{window, 2};
            cnts = per2_exemplars.(ch)(rowS:rowL,4);
            intervalanalysis.(intN).(geno).(ch).(windowname).fullcwt = wavetrans(cnts,{1,gamma,beta,fs,'bandpass'},'periodic');
            intervalanalysis.(intN).(geno).(ch).(windowname).norm =  (abs(intervalanalysis.(intN).(geno).(ch).(windowname).fullcwt))./(sum(abs(intervalanalysis.(intN).(geno).(ch).(windowname).fullcwt),'all')/fulltotal); %normalize
            intervalanalysis.(intN).(geno).(ch).(windowname).norm =   intervalanalysis.(intN).(geno).(ch).(windowname).norm(ceil(edgecut):1:ceil(end-edgecut),:);
            intervalanalysis.(intN).(geno).(ch).(windowname).PC = mean(intervalanalysis.(intN).(geno).(ch).(windowname).norm);
            %[~,~,JR,] = ridgewalk(intervalanalysis.(intN).(geno).(ch).(windowname).fullcwt(ceil(edgecut):1:ceil(end-edgecut),:),fs);
            clearvars intervalanalysis.(intN).(ch).(windowname).fullcwt
            count = 0;
            tau = vars{window, 3};
            for ridge = 1:length(intervalanalysis.(intN).(geno).(ch).(windowname).norm)     
                 tempridgepwr(ridge,1) = max(intervalanalysis.(intN).(geno).(ch).(windowname).norm(ridge,:));
                 [~,col,~] = find(intervalanalysis.(intN).(geno).(ch).(windowname).norm(ridge,:) == tempridgepwr(ridge,1));  
                 tempridge(ridge,1) = tau(col,1);
            end
            count2 = count2 + 1; 
            ridgeperiod(count2,1)= mean(tempridge);
            clearvars tempridge
            intervals(count2,1) = intN;
            windows(count2,1) = windowname;
            chans(count2,1) = convertCharsToStrings(ch); 
            genos(count2,1) = geno;
        end
    end
end
location3 = 'C:\Users\jrigg\Box\Prendergast Lab Box\Wavelet Summer 2021\Figure 3 Data\ridgetableper2.csv';
ridgetable = table(genos,windows,intervals,chans,ridgeperiod);
writetable(ridgetable, location3); 
%%
winlenlong = [30.5 7];
intaxis = [2 2 3 2];
winlenshort = [19.5 0];
for interval = 1:length(intlong)
    intL = intlong(interval);
    intS = intshort(interval);
    intN = intname(interval);
    intA = intaxis(interval);
    [rowL,~] = find(per2_exemplars.(ch)(:,1) == intL+1);
    [rowS,~] = find(per2_exemplars.(ch)(:,1) == intS);
    for chs = 1:numChs
        ch = chNames{chs, 1};
        for window = 2:3
            windowname = wname(1, window);
            wll = winlenlong(window-1);
            wls = winlenshort(window-1);
            if strcmp(ch,'ChannelBF') || strcmp(ch,'ChannelPP') == 1 || strcmp(ch,'ChannelLL') == 1 || strcmp(ch,'ChannelDC') == 1 || strcmp(ch,'ChannelDL') == 1
                geno = "Per2";
                figure 
                plot(vars{window, 3}, intervalanalysis.(intN).(geno).(ch).(windowname).PC)
                title([convertStringsToChars(windowname), ' ', convertStringsToChars(intN), ' ', convertStringsToChars(geno), ' ' convertStringsToChars(ch),  ' Powercurve'])
                axis([wls wll 0 intA])
                filename = strcat('C:\Users\jrigg\Box\Prendergast Lab Box\Wavelet Summer 2021\Figure 3 Data\', ch, windowname, geno, intN, 'powercurve');
                saveas(gcf, filename , 'eps')
                close gcf
            elseif strcmp(ch,'ChannelRU') == 1 || strcmp(ch,'ChannelRV') == 1 || strcmp(ch,'ChannelUE') == 1 || strcmp(ch,'ChannelUX') == 1 || strcmp(ch,'ChannelBQ') == 1 || strcmp(ch,'ChannelXA') == 1 || strcmp(ch,'ChannelQV') == 1
                geno = "WT";
                figure 
                plot(vars{window, 3}, intervalanalysis.(intN).(geno).(ch).(windowname).PC)
                title([convertStringsToChars(windowname), ' ', convertStringsToChars(intN), ' ', convertStringsToChars(geno), ' ' convertStringsToChars(ch),  ' Powercurve'])
                axis([wls wll 0 intA])
                filename = strcat('C:\Users\jrigg\Box\Prendergast Lab Box\Wavelet Summer 2021\Figure 3 Data\', ch, windowname, geno, intN, 'powercurve');
                saveas(gcf, filename , 'eps')
                close gcf
            end      
        end
    end
end

