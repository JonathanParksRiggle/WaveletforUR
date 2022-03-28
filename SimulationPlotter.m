load('~/justgeneratedtimeseries.mat') 
files = {'timeseries.mat', 'timeseries2.mat', 'timeseriesrand2.mat'};
fnames = {'timeseries', 'timeseries2', 'timeseriesrand2'};
%%
URtable = table(data.ur1', data.ur2', data.ur3', data.ur4');
 filename = ['~/urtable', '.csv'];
writetable(URtable, filename);
% premake variable so they will exist outside parfor
%% 
beta = 12;
for tlen = 1 :length(files)
    filenamer = files{1,tlen};
    filename = fnames{1,tlen};
    filenamer2 = strcat('/mnt/ide0/share/prendergastlab/', num2str(beta), filenamer);
    load(filenamer2)
    printer.(filename).PC = PCsreps;
    printer.(filename).Ridge = Ridgereps;
    printer.(filename).Norm = Normsreps;
    printer.(filename).Ridge2 = Ridgereps2;
    clearvars PCsreps Ridgereps Ridgereps2 Normsreps
 end
        


%% 4 real and 4 simulated activity records (one low-sh, 2 mediums, and one high)
    % (generate as TimeSeries and as Actograms)(ensure that 2 L/D transitions are clearly marked)
%3 simulated activity records plotted as TimeSeries
       %   unmodulated
      %    modulated
     %     noUR
  %……..with UR overlays for all of them

%And the spike - revover values for each of them. L and D phases.
 
f9 = figure;
sgtitle([num2str(data.ur1{1,24}), ' & ', num2str(data.ur2{1,24}), ' modulated by ', num2str(data.ur3{1,24}), ' & ', num2str(data.ur4{1,24})]);
for rep = 1:5
    subplot(5,1,rep)
    plot(simdata.timeseries2{rep,24})
    axis([0 15000 0 30])
    hold on
    %Write to clocklab
    file(1,1) = convertCharsToStrings(num2str(rep));
    file(2,1) = "01-01-2021";
    file(3,1) = "18:00";
    file(4,1) = "4";
    file(5,1) = "1";
    file(6,1) = "1";
    file(7,1) = "i";
    file(8:length(simdata.timeseries2{rep,24})+7,1) = simdata.timeseries2{rep,24};
    filename = ['~/simactogrammod' num2str(rep) '.txt'];
    writematrix(file, filename)
    
end
f8 = figure;
sgtitle([num2str(data.ur1{1,24}), ' & ', num2str(data.ur2{1,24}), ' unmodulated']);
for rep = 1:5
    subplot(5,1,rep)
    plot(simdata.timeseries{rep,24})
    axis([0 15000 0 30])
    hold on
    
    %Write to clocklab
    file(1,1) = convertCharsToStrings(num2str(rep));
    file(2,1) = "01-01-2021";
    file(3,1) = "18:00";
    file(4,1) = "4";
    file(5,1) = "1";
    file(6,1) = "1";
    file(7,1) = "i";
    file(8:length(simdata.timeseries2{rep,24})+7,1) = simdata.timeseries2{rep,24};
    filename = ['~/simactogramunmod' num2str(rep) '.txt'];
    writematrix(file, filename)
end
        %%
% tlen -> sim -> rep -> lite -> variable
 Lurs = num2str(mean(printer.timeseries.Ridge2{1,24}.light.Medium));
Lurl = num2str(mean(printer.timeseries.Ridge2{1,24}.light.Long));
Durs = num2str(mean(printer.timeseries.Ridge2{1,24}.dark.Medium));
Durl = num2str(mean(printer.timeseries.Ridge2{1,24}.dark.Long));

Lurs2 = num2str(mean(printer.timeseries2.Ridge2{1,24}.light.Medium));
Lurl2 = num2str(mean(printer.timeseries2.Ridge2{1,24}.light.Long));
Durs2 = num2str(mean(printer.timeseries2.Ridge2{1,24}.dark.Medium));
Durl2 = num2str(mean(printer.timeseries2.Ridge2{1,24}.dark.Long));

Lurs3 = num2str(mean(printer.timeseriesrand2.Ridge2{1,24}.light.Medium));
Lurl3 = num2str(mean(printer.timeseriesrand2.Ridge2{1,24}.light.Long));
Durs3 = num2str(mean(printer.timeseriesrand2.Ridge2{1,24}.dark.Medium));
Durl3 = num2str(mean(printer.timeseriesrand2.Ridge2{1,24}.dark.Long));


ts=60;% replicate the 1 minute sampling
T=14400*60; %total time in seconds
ti = 1:ts:T;% create a time stream        
ult  = square(2*pi*(1/((data.ur1{1,24})*3600))*ti)';
ult2  =  square(2*pi*(1/((data.ur2{1,24})*3600))*ti)';
ult3 = square(2*pi*(1/((data.ur3{1,24})*3600))*ti)';
ult4 = square(2*pi*(1/((data.ur4{1,24})*3600))*ti)';
f7 = figure;
subplot(3,1,1)
plot(simdata.timeseries{1,24}(1,721:2160),'b')
hold on
subplot(3,1,2)
plot(ult(721:2160,1) * 10, 'r') 
axis([0 1500 0 20])
subplot(3,1,3)
plot(ult2(721:2160,1) * 10, 'y')
axis([0 1500 0 20])
sgtitle([num2str(data.ur1{1,24}), ' & ', num2str(data.ur2{1,24}), ' unmodulated']);
xlabel(['Recovered Circadian on: ', Durl, ' & ', Durs, ' And ', 'Recovered Circadian off: ', Lurl, ' & ', Lurs]) 
fname = 'unmodtimeseriesexamplefig2';
filenamer1 = '/home/jpriggle/acropolis_home/';
saveas(f7, fullfile(filenamer1, fname), 'epsc')

f6 = figure; 
subplot(5,1,1)
plot(simdata.timeseries2{1,24}(1,721:2161),'b')
axis([0 1500 0 20])
hold on
subplot(5,1,2)
plot(ult(721:2160,1) * 10, 'r')
axis([0 1500 0 20])
subplot(5,1,3)
plot(ult2(721:2160,1) * 10, 'y')
axis([0 1500 0 20])
subplot(5,1,4)
plot(ult3(721:2160,1) * 10, 'm')
axis([0 1500 0 20])
subplot(5,1,5)
plot(ult4(721:2160,1) * 10, 'c')
axis([0 1500 0 20])
sgtitle([num2str(data.ur1{1,24}), ' & ', num2str(data.ur2{1,24}), ' modulated by ', num2str(data.ur3{1,24}), ' & ', num2str(data.ur4{1,24})]);
xlabel(['Recovered Circadian on: ', Durl2, ' & ', Durs2, ' And ', 'Recovered Circadian off: ', Lurl2, ' & ', Lurs2 ])  %make this an easy variable
fname = 'modtimeseriesexamplefig2';
saveas(f6, fullfile(filenamer1, fname), 'epsc')

f5 = figure; 
subplot(5,1,1)
plot(simdata.timeseries2rand{1,24}(1,721:2160),'b')
hold 
subplot(5,1,2)
plot(ult(721:2160,1) * 10, 'r')
axis([0 1500 0 30])
subplot(5,1,3)
plot(ult2(721:2160,1) * 10, 'y')
axis([0 1500 0 30])
subplot(5,1,4)
plot(ult3(721:2160,1) * 10, 'm')
axis([0 1500 0 30])
subplot(5,1,5)
plot(ult4(721:2160,1) * 10, 'c')
axis([0 1500 0 30])
sgtitle(['No URs:', num2str(data.ur1{1,24}), ' & ', num2str(data.ur2{1,24}),' modulated by ', num2str(data.ur3{1,24}), ' & ', num2str(data.ur4{1,24}), 'randomized']);
xlabel(['Recovered Circadian on: ', Durl3, ' & ', Durs3, ' And ', 'Recovered Circadian off: ', Lurl3, ' & ', Lurs3 ]);
fname = 'noURtimeseriesexamplefig2';
saveas(f5, fullfile(filenamer1, fname), 'epsc')

%%
%ADD  A criterion to individual ridgepoints if they are at edge of windows
%tau START OR TAU END
%T
%or perhaps too far from mean
litenames = fieldnames(printer.timeseries.Ridge2{1,1});
 reps = 10;
 files = {'timeseries', 'timeseries2','timeseriesrand2'};
 ur = {'ur1', 'ur2', 'ur3', 'ur4'};
for tlen = 1:length(files)
    time = files{1,tlen};
    countl = 0;
    countd = 0;
    countf = 0;
    countfd = 0;
    countfl = 0;
    counts = {countl,countd, countf, countfd, countfl};
    for sim = 1:100
        ur1 = data.ur1{1,sim};
        ur2 = data.ur2{1,sim};
        ur3 = data.ur3{1,sim};
        ur4 = data.ur4{1,sim};
        for rep = 1:reps
            for lite = 1:length(litenames)
                lites = litenames{lite,1};
                count = counts{1,lite};
                if tlen == 1
                 
                   count = count + 1;
                   recovery.(time).(lites).real(count,:) = ur1;
                   recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Long);
                   recovery.(time).(lites).rm(count,:) =  mean(printer.(time).Ridge2{rep,sim}.(lites).Long);
                   count = count + 1;
                   recovery.(time).(lites).real(count,:) = ur2; 
                   recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Medium);
                   recovery.(time).(lites).rm(count,:) = mean(printer.(time).Ridge2{rep,sim}.(lites).Medium);
                   counts{1,lite} = count;
                elseif tlen > 1
                    if lite == 1 || lite == 5
                       count = count + 1;
                       recovery.(time).(lites).real(count,:) = ur3;
                       recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Long);
                       recovery.(time).(lites).rm(count,:) =  mean(printer.(time).Ridge2{rep,sim}.(lites).Long);
                       count = count + 1;
                       recovery.(time).(lites).real(count,:) = ur4; 
                       recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Medium);
                       recovery.(time).(lites).rm(count,:) = mean(printer.(time).Ridge2{rep,sim}.(lites).Medium);
                       counts{1,lite} = count;
                    elseif lite == 2 || lite == 4
                       count = count + 1;
                       recovery.(time).(lites).real(count,:) = ur1;
                       recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Long);
                       recovery.(time).(lites).rm(count,:) =  mean(printer.(time).Ridge2{rep,sim}.(lites).Long);
                       count = count + 1;
                       recovery.(time).(lites).real(count,:) = ur2;
                       recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Medium);
                       recovery.(time).(lites).rm(count,:) = mean(printer.(time).Ridge2{rep,sim}.(lites).Medium);
                       counts{1,lite} = count;
                    elseif lite == 3
                       count = count + 1;
                       recovery.(time).(lites).real(count,:) = ur1;
                       recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Long);
                       recovery.(time).(lites).rm(count,:) =  mean(printer.(time).Ridge2{rep,sim}.(lites).Long);
                       count = count + 1;
                       recovery.(time).(lites).real(count,:) = ur2;
                       recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Medium);
                       recovery.(time).(lites).rm(count,:) = mean(printer.(time).Ridge2{rep,sim}.(lites).Medium);
                       count = count + 1;
                       recovery.(time).(lites).real(count,:) = ur3;
                       recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Long);
                       recovery.(time).(lites).rm(count,:) =  mean(printer.(time).Ridge2{rep,sim}.(lites).Long);
                       count = count + 1;
                       recovery.(time).(lites).real(count,:) = ur4; 
                       recovery.(time).(lites).rw(count,:) =  nanmean(printer.(time).Ridge{rep,sim}.(lites).Medium);
                       recovery.(time).(lites).rm(count,:) = mean(printer.(time).Ridge2{rep,sim}.(lites).Medium);
                       counts{1,lite} = count;
                    end
                end
            end
        end
    end
end
%%
for tlen = 1:length(files)
    time = files{1,tlen};
    for lite = 1:length(litenames)
        lites = litenames{lite,1};
        filename = ['~/recovery' time, lites, '.csv'];
        real = recovery.(time).(lites).real;
        rw = recovery.(time).(lites).rw;
        rm = recovery.(time).(lites).rm;
        rec = table(real, rw, rm);
        writetable(rec, filename);
        figure
        scatter(real, rw)
        title(['real vs ridgewalk-' time lites])
        axis([0 5 0 5])
        figure
        scatter(real, rm)
         title(['real vs ridgem-' time lites])
         axis([0 5 0 5])
    end
end

%%
litenames = fieldnames(printer.timeseries.Ridge{1,1});
 reps = 10;
 files = {'timeseries', 'timeseries2', 'timeseriesrand2'};
for tlen = 1:length(files)
    time = files{1,tlen};
    for sim = 25
        for rep = 1:reps
            for lite = 1:length(litenames)
                lites = litenames{lite,1};
                powercurves.(time).(lites){rep,1}= printer.(time).PC{rep,sim}.(lites).URFull(1,:);
                
            end
        end
    end
end

%%
tau = vars{1,2};
for tlen = 1:length(files)
    time = files{1,tlen};
    
    for lite = 1:length(litenames)
         lites = litenames{lite,1};
         powercurves.(time).(lites) = cell2mat(powercurves.(time).(lites));
         [bootCIs,bootstat] = bootci(2000,@mean,powercurves.(time).(lites));
         boot.(time).(lites).bootCI = bootCIs;
         boot.(time).(lites).bootstat = bootstat;
         f11 = figure;
         plot(tau, bootCIs')
         hold on
         plot(tau, mean(bootstat))
         title([time, ' ', lites])
         axis([0 7 0 0.035])
         fname = (['~/powercurves', time, lites, num2str(sim)]);
         saveas(f11, fname, 'epsc')
      
    end
end

%%
for tlen = 1:length(files)
    time = files{1,tlen};
    for lite = 1:length(litenames)
        lites = litenames{lite,1};
        vname = 'URFull';
        figure ('visible', 'off')
        h = pcolor((1:length(printer.(time).Norm{1,24}.(lites).(vname)))/1440, vars{1,2}, printer.(time).Norm{1,24}.(lites).(vname)' );
        colormap hot
        set(h, 'EdgeColor', 'none')
        set(gca,'TickDir','out');
        filename = (['~/scalograms', time, lites, vname]);
        saveas(gcf, filename , 'tiffn')
    end
end
