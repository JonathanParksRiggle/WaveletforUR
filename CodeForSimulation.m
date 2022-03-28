%load('\Users\jrigg\Box\Wavelet\JP wavelets\Leslie Wavelet\wavelet procedures\ultradDataGamma5beta9.mat');
load('home/jpriggle/acropolis_home/Ultrad.mat') 
ts=60;% replicate the 1 minute sampling
T=14400*60; %total time in seconds
ti = 1:ts:T;% create a time stream
f = 1/(24*3600); %once a day in hz
numberofsimulations= 100;
ur = {'ur1', 'ur2', 'ur3', 'ur4'}; 
for sim = 1:numberofsimulations
    for urs = 1:length(ur)
        urnum = ur{1,urs};
        if strcmp(urnum, 'ur1') > 0 || strcmp(urnum, 'ur3') > 0
            data.(urnum){1,sim} = 2.13 + (4.26-2.13) .* rand(1,1); %generate a random number within range to serve as UR period
        elseif strcmp(urnum, 'ur2') > 0 || strcmp(urnum, 'ur4') > 0 %convert to hz
            data.(urnum){1,sim} = 1.07 + (2.13-1.07) .* rand(1,1); %generate a random number within range to serve as UR period
        end
    end      
end



%% 
reps = 10;
timeseries = zeros(1,14400); 
timeseries2 = zeros(1,14400); 
rtimeseries = zeros(1,14400); 
rtimeseries2 = zeros(1,14400);

for sim = 1:numberofsimulations
    
    circmodulatornoisegrp = rand;
        if circmodulatornoisegrp >= .50
           circgrphigh= 1;
           circgrplow = .6;
        elseif circmodulatornoisegrp <= .50 
            circgrphigh =.6;
            circgrplow = .2;
        end
    
    
    for rep = 1:reps
        URwiggle1 = (-.008 * data.ur1{1,sim}) + ((.008 * data.ur1{1,sim})-(-.008 * data.ur1{1,sim})) .*rand(1,1); %scale each of these ranges appropriately to the URs
        URwiggle2 = (-.008 * data.ur2{1,sim}) + ((.008 * data.ur2{1,sim})-(-.008 * data.ur2{1,sim})) .*rand(1,1) ; 
        URwiggle3 = (-.008 * data.ur3{1,sim}) + ((.008 * data.ur3{1,sim})-(-.008 * data.ur3{1,sim})) .*rand(1,1);
        URwiggle4 = (-.008 * data.ur4{1,sim}) + ((.008 * data.ur4{1,sim})-(-.008 * data.ur4{1,sim})) .*rand(1,1);
        circa  = square(2*pi*(1/((24)*3600))*ti); %once a day in hz
        ult  = square(2*pi*(1/((data.ur1{1,sim} + URwiggle1)*3600))*ti);
        ult2  =  square(2*pi*(1/((data.ur2{1,sim} + URwiggle2)*3600))*ti);
        ult3 = square(2*pi*(1/((data.ur3{1,sim} + URwiggle3)*3600))*ti);
        ult4 = square(2*pi*(1/((data.ur4{1,sim} + URwiggle4)*3600))*ti);
        frand = 1/(rand*3600); %high frequency noise
        noise = square(2*pi*frand*ti);
        
        for timelength = 1:length(circa)
            if (noise(1,timelength) == 1) %start with noise
                rtimeseries(1,timelength) = rand;
                rtimeseries2(1,timelength) = rand;
                
            elseif (noise(1,timelength) == -1)
                rtimeseries(1,timelength) = 0;
                rtimeseries2(1,timelength) = 0;

            end
        end
        
       
        circmodulator = circgrplow + (circgrphigh-circgrplow) .* rand(1,1);
            
    
        % now signal
        for timelength = 1:length(circa)
            if (circa(1,timelength) == -1)%Circadian off
                if (ult(1,timelength) == -1)  && (ult2(1,timelength) == -1) 
                     random = rand;
                   if (random > .01) 
                      timeseries(1,timelength) = 0 ;
                   elseif (random <= .01)
                      timeseries(1,timelength) = randi([0 30])  ;
                   end
                elseif (ult(1,timelength) == 1)  && (ult2(1,timelength) == -1) 
                     random = rand;
                   if (random > .20) && (random < .99)
                      timeseries(1,timelength) = 0 ;                  
                   elseif (random <= .20) 
                      timeseries(1,timelength) = randi([0 20]) ;
                   elseif (random >= .99)
                      timeseries(1,timelength) = randi([0 30])  ;
                   end
                elseif (ult(1,timelength) == -1)  && (ult2(1,timelength) == 1) 
                     random = rand;
                    if (random > .20) && (random < .99)
                        timeseries(1,timelength) = 0 ;
                    elseif (random <= .20)
                        timeseries(1,timelength) = randi([0 20]) ;
                    elseif (random >= .99)
                        timeseries(1,timelength) = randi([0 30])  ;
                    end
                elseif (ult(1,timelength) == 1)  && (ult2(1,timelength) == 1)
                    random = rand;
                    if (random >= .3334) && (random < .99)
                      timeseries(1,timelength) = 0 ;
                    elseif (random < .3334) && (random > .20)
                      timeseries(1,timelength) = randi([0 30]) ; 
                    elseif (random <= .20)
                      timeseries(1,timelength) = randi([0 20]) ;
                    elseif (random >= .99)
                      timeseries(1,timelength) = randi([0 30])  ;
                    end 
                end
                if (ult3(1,timelength) == -1)  && (ult4(1,timelength) == -1) 
                   random = rand;
                   if (random <= .01) 
                      timeseries2(1,timelength) = randi([0 30]) ;
                   elseif (random > .01) 
                      timeseries2(1,timelength) = 0 ;
                   end
                elseif (ult3(1,timelength) == -1)  && (ult4(1,timelength) == 1) 
                     random = rand;
                   if (random > .20) && (random < .99)
                      timeseries2(1,timelength) = 0 ;
                   elseif (random <= .20)
                      timeseries2(1,timelength) = randi([0 20])  ;
                   elseif (random >= .99)
                      timeseries2(1,timelength) = randi([0 30])  ;
                   end
                elseif (ult3(1,timelength) == 1)  && (ult4(1,timelength) == -1) 
                     random = rand;
                   if (random > .20) && (random < .99)
                      timeseries2(1,timelength) = 0 ;
                   elseif (random <= .20) 
                      timeseries2(1,timelength) = randi([0 20]) ;  
                   elseif (random >= .99)
                      timeseries2(1,timelength) = randi([0 30])  ;
                   end
                elseif (ult3(1,timelength) == 1)  && (ult4(1,timelength) == 1)
                    random = rand;
                    if (random >= .3334) && (random < .99)
                      timeseries2(1,timelength) = 0 ;
                    elseif (random < .3334) && (random > .20)
                      timeseries2(1,timelength) = randi([0 30])  ; 
                   elseif (random <= .20)
                      timeseries2(1,timelength) = randi([0 20]) ;
                   elseif (random >= .99)
                      timeseries2(1,timelength) = randi([0 30])  ;
                   end
                end
            elseif circa(1,timelength) == 1 
                if (ult(1,timelength) == -1)  && (ult2(1,timelength) == -1) 
                   random = rand;
                   if (random > .50) && (random <.99)
                      timeseries(1,timelength) = 0 ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                   elseif (random <= .50) 
                      timeseries(1,timelength) = randi([0 20]) ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                   elseif (random >= .99)
                      timeseries(1,timelength) = randi([0 30])  ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                   end
                elseif (ult(1,timelength) == -1)  && (ult2(1,timelength) == 1) 
                     random = rand;
                   if (random > .70) && (random < .99)
                      timeseries(1,timelength) = 0 ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                   elseif (random <= .70)
                      timeseries(1,timelength) = randi([0 20]);
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                   elseif (random >= .99)
                      timeseries(1,timelength) = randi([0 30])  ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                   end
                elseif (ult(1,timelength) == 1)  && (ult2(1,timelength) == -1) 
                     random = rand;
                   if (random > .70) && (random < .99)
                      timeseries(1,timelength) = 0 ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                   elseif (random <= .70) 
                      timeseries(1,timelength) = randi([0 20]);
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                   elseif (random >= .99)
                      timeseries(1,timelength) = randi([0 30])  ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                   end
                elseif (ult(1,timelength) == 1)  && (ult2(1,timelength) == 1)
                    random = rand;
                    if (random >= .90) && (random < .99)
                      timeseries(1,timelength) = 0 ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                    elseif (random >= .72) && (random < .90)
                      timeseries(1,timelength) = randi([0 30]);  
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                    elseif (random < .72) && (random >= .18)
                      timeseries(1,timelength) = randi([0 20]) ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                    elseif (random < .18)
                      timeseries(1,timelength) = randi([0 10]) ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                    elseif (random >= .99)
                      timeseries(1,timelength) = randi([0 30])  ;
                      timeseries2(1,timelength) = timeseries(1,timelength) ;
                    end
                end
            end
        end
        timeseries = timeseries - rtimeseries;
        timeseries2= timeseries2 - rtimeseries;
        simdata.timeseries{rep,sim} = timeseries * circmodulator;
        simdata.timeseries2{rep,sim} = timeseries2 * circmodulator;
    end
end
%% get rid of negatives

for sim = 1:numberofsimulations
    for rep = 1:reps
      simdata.timeseries{rep,sim}(1,simdata.timeseries{rep,sim} < 0) = 0;
      simdata.timeseries2{rep,sim}(1,simdata.timeseries2{rep,sim} < 0) = 0;
    end
end
%%
f9 = figure;
sgtitle('modulated');
for rep = 1:reps/2
    
    subplot(reps/2,1,rep)
    plot(simdata.timeseries2{rep,5})
    axis([0 15000 0 30])
    hold on
    file(1,1) = convertCharsToStrings(num2str(rep));
    file(2,1) = "01-01-2021";
    file(3,1) = "18:00";
    file(4,1) = "4";
    file(5,1) = "1";
    file(6,1) = "1";
    file(7,1) = "i";
    file(8:length(simdata.timeseries2{rep,1})+7,1) = simdata.timeseries2{rep,1};
    filename = ['~/simactogrammodtest' num2str(rep) '.txt'];
    writematrix(file, filename)
end

f8 = figure;
sgtitle('unmodulated');
for rep = 1:reps/2
   
    subplot(reps/2,1,rep)
    plot(simdata.timeseries{rep,5})
    axis([0 15000 0 30])
    hold on
    
     file(1,1) = convertCharsToStrings(num2str(rep));
    file(2,1) = "01-01-2021";
    file(3,1) = "18:00";
    file(4,1) = "4";
    file(5,1) = "1";
    file(6,1) = "1";
    file(7,1) = "i";
    file(8:length(simdata.timeseries2{rep,1})+7,1) = simdata.timeseries2{rep,1};
    filename = ['~/simactogramunmodtest' num2str(rep) '.txt'];
    writematrix(file, filename)
end



filenamer1 = '/home/jpriggle/acropolis_home/';
fname = 'simunmodtime';
saveas(f8, fullfile(filenamer1, fname), 'epsc')

filenamer1 = '/home/jpriggle/acropolis_home/';
fname = 'simmodtime';
saveas(f9, fullfile(filenamer1, fname), 'epsc')
%%
%sex = readmatrix('/home/jpriggle/acropolis_home/JP wavelets/Leslie Wavelet/ListOfSexesForJP');
%chnames = fieldnames(ultrad);

%for chs = 1: length(chnames)
  %  ch = chnames{chs,1};
  %  if strcmp('ch183',ch) == 1 || strcmp('ch199',ch) == 1 || strcmp('ch210',ch) == 1 || strcmp('ch217',ch) == 1 || strcmp('ch224',ch) == 1 ||  strcmp('ch280',ch) == 1
    %    continue
  %  else
   %     newUR.(ch) = ultrad.(ch).intDL.W08.Cntmin;
   %     if sex(chs,2) == 0
   %         sexUR.(ch) = 'male';
   %     elseif sex(chs,2) == 1
   %         sexUR.(ch) = 'female';
     %   end
    %end
%end
%%
% 
% newchnames = fieldnames(newUR);
% count1 = 1;
% count2 = 1;
% count3 = 1;
% count4 = 1;
% count5 = 1;
% count6 = 1;
% count7 = 1;
% f1= figure;
% f2 = figure;
% f3 = figure;
% f4 = figure;
% f5 = figure;
% f6 = figure;
% f7 = figure;
% for chs = 1:length(newchnames)
%     ch = newchnames{chs,1};
%     if chs <= 5
%         figure(f1)
%         subplot(5, 1, count1)
%         plot(newUR.(ch))
%         axis([0 15000 0 30])
%         count1 = 1 + count1;
%         title([ch, ' ', sexUR.(ch)])
%         sgtitle('first')
%       elseif chs >= 5 && chs <= 10
%         figure(f2)
%         subplot(5, 1, count2)
%         plot(newUR.(ch))
%         axis([0 15000 0 30])
%         count2 = 1 + count2;
%         title([ch, ' ', sexUR.(ch)])
%         sgtitle('Second')
%     elseif chs >= 10 && chs <= 15 
%         figure(f3)
%         subplot(5, 1, count3)
%         plot(newUR.(ch))
%         axis([0 15000 0 30])
%         count3 = 1 + count3;
%         title([ch, ' ', sexUR.(ch)])
%         sgtitle('third')
%     elseif chs >= 15 && chs <= 20
%         figure(f4)
%         subplot(5, 1, count4)
%         plot(newUR.(ch))
%         axis([0 15000 0 30])
%         count4 = 1 + count4;
%         title([ch, ' ', sexUR.(ch)])
%         sgtitle('four')
%     elseif chs >= 20 && chs <= 25
%         figure(f5)
%         subplot(5, 1, count5)
%         plot(newUR.(ch))
%         axis([0 15000 0 30])
%         count5 = 1 + count5;
%         title([ch, ' ', sexUR.(ch)])
%         sgtitle('five')
%     elseif chs >= 25 && chs <= 30
%         figure(f6)
%         subplot(5, 1, count6)
%         plot(newUR.(ch))
%         axis([0 15000 0 30])
%         count6 = 1 + count6;
%         title([ch, ' ', sexUR.(ch)])
%         sgtitle('six')
%     elseif chs >= 30  && chs <= length(newchnames)
%         figure(f7)
%         subplot(5, 1, count7)
%         plot(newUR.(ch))
%         axis([0 15000 0 30])
%         count7 = 1 + count7;
%         title([ch, ' ', sexUR.(ch)])
%         sgtitle('seven')
%     end
% end
% filenamer1 = '/home/jpriggle/acropolis_home/';
% fname = 'realbatch1';
% saveas(f1, fullfile(filenamer1, fname), 'epsc')
% fname = 'realbatch2';
% saveas(f2, fullfile(filenamer1, fname), 'epsc')
% fname = 'realbatch3';
% saveas(f3, fullfile(filenamer1, fname), 'epsc')
% fname = 'realbatch4';
% saveas(f4, fullfile(filenamer1, fname), 'epsc')
% fname = 'realbatch5';
% saveas(f5, fullfile(filenamer1, fname), 'epsc')
% fname = 'realbatch6';
% saveas(f6, fullfile(filenamer1, fname), 'epsc')
% fname = 'realbatch7';
% saveas(f7, fullfile(filenamer1, fname), 'epsc')
%%
for sim = 1:25 % look at a 1/4 of them 
    figure
    for rep = 1:reps/2
        subplot(reps/2,1,rep)
        plot(simdata.timeseries{rep,sim})
        axis([0 15000 0 30])
        hold on
    end
end
%% create the timeseries with random noise 

timeseries = {'timeseries' 'timeseries2'};
for time = 1:length(timeseries) %just want the half of the time series with actual ultradian signal
    timeser = timeseries{1,time};  
    for sim = 1:numberofsimulations
        for rep = 1:(reps)
            startL = 1;
            stopL = 720;
            startD = 721;
            stopD = 1440;
            temp = simdata.(timeser){rep,sim}; %create my composite timeseries  
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
            if strcmp(timeser,'timeseries') > 0
                simdata.timeseries1rand{rep,sim} = [day1L day1D day2L day2D day3L day3D day4L day4D day5L day5D day6L day6D day7L day7D day8L day8D day9L day9D day10L day10D];
            elseif strcmp(timeser,'timeseries2') > 0
                simdata.timeseries2rand{rep,sim} = [day1L day1D day2L day2D day3L day3D day4L day4D day5L day5D day6L day6D day7L day7D day8L day8D day9L day9D day10L day10D];
            end  
            clear  temp day1L day1D day2L day2D day3L day3D day4L day4D day5L day5D day6L day6D day7L day7D day8L day8D day9L day9D day10L day10D
        end
    end
end
filename = '~/justgeneratedtimeseries.mat' ;
save(filename, 'simdata', 'data');

