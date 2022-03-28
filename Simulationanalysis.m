
%% Seperate light and dark, perform wavelet, normalize, perform ridge analysis, and group period power curves and ffts 
tic
load('~/justgeneratedtimeseries.mat') 

NstepsPerHr=60; % sampling rate
Nsteps=24*10*NstepsPerHr+1;  % total number of samples
t=(0:Nsteps-1)'/NstepsPerHr; % creates timept vector
T=t(Nsteps)-t(1);            % time at last time point [not sure why this is needed]
% change following for dark period of long days to 6.5 longestperiod
longestperiods = ([6.5, 1.07, 2.13, 4.26]);
shortestperiods = ([.5, .53, 1.07, 2.13]);

nvoices = ([27, 99, 101, 100]);
empty = {{}, {}, {}, {}} ;
vars = empty;
edge = empty;
for window = 1:length(longestperiods)
    longestperiod=longestperiods(1, window);shortestperiod=shortestperiods(1,window);  % range of periods to use in AWT
    gamma=3;beta=4;
    P = sqrt(beta*gamma);
    nvoice=nvoices(1,window); 
    [fs,tau,qscaleArray] = CalcScaleForCWT(shortestperiod,longestperiod,T,NstepsPerHr,nvoice);
    IntDLScale = length(qscaleArray) * length(tau); %size of the matrix
    [row,~,~] = find(tau == tau(tau >= shortestperiod & tau <= longestperiod)'); %allows you to set the true window you want
    qscaleArray = qscaleArray(row:end,1);
    fs = fs(row:end,1);
    tau = tau(row:end,1); 
    clear col row
    vars{window, 1} = qscaleArray;
    vars{window, 2} = tau;
    vars{window, 3} = fs;
    vars{window,4} = IntDLScale;
    edge{window, 1} = round((longestperiod * 1.5) * 60, 0);
end
%% Splitting the data before
lites = {'L' 'D' 'F'};
timeseries = {'timeseries' 'timeseries2' 'timeseries2rand'};
ts=60;% replicate the 1 minute sampling
T=14400*60; %total time in seconds
ti = 1:ts:T;% create a time stream
circa  = square(2*pi*(1/((24)*3600))*ti); %once a day in hz generated the same way as when simulating timeseries
idx = circa;
clear ultrad circa
numberofsimulations =100;
reps = 10;
counts = {'LCounts' 'DCounts' 'FCounts'}; 
 
for time = 1:length(timeseries) 
    timeser = timeseries{1,time};
    for sim = 1:numberofsimulations
        for rep = 1:reps
            for lite = 1:length(counts)
                count = counts{1,lite};
                if strcmp(count,'LCounts') == 1
                    light.(timeser){rep,sim} = simdata.(timeser){rep,sim}(idx == -1);
                elseif strcmp(count,'DCounts') == 1
                    dark.(timeser){rep,sim} = simdata.(timeser){rep,sim}(idx == 1);
                elseif strcmp(count,'FCounts') == 1
                    full.(timeser){rep,sim} = simdata.(timeser){rep,sim}; 
                end 
            end
        end
    end
end
save('countsimdata.mat', 'dark', 'light', 'full', 'edge', 'idx')
clearvars -except idx vars timeseries edge reps numberofsimulations gamma beta lites P

%
 %% Runs the wavelet and ridge on the time series
files = {'timeseries.mat', 'timeseries2.mat', 'timeseriesrand2.mat'};
% premake variable so they will exist outside parfor

for tlen = 1 :length(files)
 
    PCs = {'lPC', 'dPC', 'fPC', 'fPCDark', 'fPCLight'}; 
    Norms = {'lNorm', 'dNorm', 'fnorm',  'fNormDark', 'fNormLight'}; 
    Ridges = {'lRidge', 'dRidge', 'fRidge',  'fRidgeDark', 'fRidgeLight'}; 
    empty = {{}, {}, {}, {}, {}} ;
    PC = cell2struct(empty,PCs,2);
    Norm = cell2struct(empty,Norms,2);
    PCsreps = cell(length(numberofsimulations),length(reps));
    Normsreps = cell(length(numberofsimulations),length(reps));
    tempPC = cell(numberofsimulations);
    tempNorm = cell(numberofsimulations);
    ridgeperiod = cell(numberofsimulations);
    ridgeperiod2 = cell(numberofsimulations);
    tempPCs = cell(numberofsimulations);
    tempNorms = cell(numberofsimulations);
    tempRidges = cell(numberofsimulations);
    tempRidges2 = cell(numberofsimulations);
    filenamer = files{1,tlen}; 
    timeser = timeseries{1,tlen};
    load('countsimdata.mat')
    tempdata{1,1} = light.(timeser);
    tempdata{2,1}  = dark.(timeser);
    tempdata{3,1} = full.(timeser);

    parfor sim = 1:numberofsimulations
      
        tempdats = tempdata;
         count = 0; 
        for rep = 1:reps
            nbedge = edge;
            nbidx = idx;
            for lite = 1:length(PCs)
                count = count + 1; 
               
                lights = {'light', 'dark', 'full', 'fullDark', 'fullLight'}; 
                lites = lights{1,lite};
                if lite <= 3 %presplit or unsplit
                    tempdat = tempdats{lite,1};
                    columnofinterest = tempdat{rep, sim};
                    tempdat = {}; 
                    for variables = 1:size(vars,1)
                        vnames = ["URFull", "Short", "Medium", "Long" ];
                        vname = vnames{1,variables};
                        tau = vars{variables, 2};
                        IntDLScale = vars{variables,4};
                        edgecut = edge{variables,1};
                        temp= wavetrans(columnofinterest', {1,gamma,beta,vars{variables, 3},'bandpass'},'periodic');
                        tempavg= sum(sum(abs(temp)))/IntDLScale; %switched from abs here
                        temp2 = abs(temp)./tempavg;
                        tempavg = [];
                        tempsize = size(vars{variables, 2},1);
                        temp3 = temp2(edgecut:(end-edgecut), 1:tempsize);
                        temp4 = temp(edgecut:(end-edgecut), 1:tempsize);
                        tempPC{1,sim}.(vname) = [];
                        tempPC{1,sim}.(vname)(:, 1:tempsize)  = mean(abs(temp3));
                        if sim == 24
                            tempNorm{1,sim}.(vname) = [];
                            tempNorm{1,sim}.(vname) = temp3;
                        end
                        [~,~,JR,]= ridgewalk(temp4,vars{variables, 3});
                        count = 0;
                        ridgeperiod2{1,sim}.(vname) =[];
                        for ridge2 = 1:length(temp3)
                            ridgeperiod2{1,sim}.(vname)(ridge2,1) = max(temp3(ridge2,:));
                            [~,col,~] = find(temp3(ridge2,:) == ridgeperiod2{1,sim}.(vname)(ridge2,1));  
                            ridgeperiod2{1,sim}.(vname)(ridge2,1) = tau(col,1);
                        end
                        ridgeperiod{1,sim}.(vname) = [];
                        for ridge = 1:length(JR)
                            if isnan(JR(ridge,1)) == true 
                                count = count + 1;
                                ridgeperiod{1,sim}.(vname)(1, count) = JR(ridge,1);
                            elseif isnan(JR(ridge,1)) == false 
                                count = count + 1;
                                ridgeperiod{1,sim}.(vname)(1, count) = tau(JR(ridge,1),1);
                            end
                        end
                        temp = {}; temp2 = {}; temp3 = {}; JR = {};
                    end
                elseif lite == 4 %split after wavelet
                    tempdat = tempdats{3,1};
                    columnofinterest = tempdat{rep, sim};
                    tempdat = {}; 
                    for variables = 1:size(vars,1)
                        vnames = ["URFull", "Short", "Medium", "Long" ];
                        vname = vnames{1,variables};
                        tau = vars{variables, 2};
                        IntDLScale = vars{variables,4};
                        edgecut = edge{variables,1};
                        temp= wavetrans(columnofinterest', {1,gamma,beta,vars{variables, 3},'bandpass'},'periodic');
                        tempavg= sum(sum(abs(temp)))/IntDLScale;
                        temp2 = (abs(temp)./tempavg);
                        tempavg = [];
                        tempsize = size(vars{variables, 2},1);
                        tempidx = temp2(idx == 1, 1:length(tau));
                        tempidx2 = temp(idx == 1, 1:length(tau));
                        tempidx = tempidx(edgecut:(end-edgecut), 1:tempsize);
                        tempidx2 = tempidx2(edgecut:(end-edgecut), 1:tempsize);
                        tempPC{1,sim}.(vname) = [];
                        tempPC{1,sim}.(vname)(:, 1:tempsize)  = mean(abs(tempidx));
                        if sim == 24   
                            tempNorm{1,sim}.(vname) = [];
                            tempNorm{1,sim}.(vname)  = tempidx;
                        end
                        [~,~,JR,]= ridgewalk(tempidx2,vars{variables, 3});
                        count = 0;
                        ridgeperiod2{1,sim}.(vname) =[];
                         for ridge2 = 1:length(tempidx)
                            ridgeperiod2{1,sim}.(vname)(ridge2,1) = max(tempidx(ridge2,:));
                            [~,col,~] = find(tempidx(ridge2,:) == ridgeperiod2{1,sim}.(vname)(ridge2,1));  
                            ridgeperiod2{1,sim}.(vname)(ridge2,1) = tau(col,1);
                        end
                        for ridge = 1:length(JR)
                            if isnan(JR(ridge,1)) == true 
                                count = count + 1;
                                ridgeperiod{1,sim}.(vname)(1, count) = JR(ridge,1);
                            elseif isnan(JR(ridge,1)) == false 
                                count = count + 1;
                                ridgeperiod{1,sim}.(vname)(1, count) = tau(JR(ridge,1),1);
                            end
                        end
                        temp = {}; temp2 = {}; tempidx = {}; JR = {};
                    end
                elseif lite == 5 
                    tempdat = tempdats{3,1};
                    columnofinterest = tempdat{rep, sim};
                    tempdat = {}; 
                    for variables = 1:size(vars,1)
                        vnames = ["URFull", "Short", "Medium", "Long" ];
                        vname = vnames{1,variables};
                        tau = vars{variables, 2};
                        IntDLScale = vars{variables,4};
                        edgecut = edge{variables,1};
                        temp= wavetrans(columnofinterest', {1,gamma,beta,vars{variables, 3},'bandpass'},'periodic');
                        tempavg= sum(sum(abs(temp)))/IntDLScale;
                        temp2 = (abs(temp)./tempavg);
                        tempavg = [];
                        tempsize = size(vars{variables, 2},1);
                        tempidx = temp2(idx == -1, 1:length(tau));
                        tempidx2 = temp(idx == -1, 1:length(tau));
                        tempidx = tempidx(edgecut:(end-edgecut), 1:tempsize);
                        tempidx2 = tempidx2(edgecut:(end-edgecut), 1:tempsize);
                        tempPC{1,sim}.(vname) = [];
                        tempPC{1,sim}.(vname)(:, 1:tempsize)  = mean(abs(tempidx));
                        if sim == 24
                            tempNorm{1,sim}.(vname) = [];
                            tempNorm{1,sim}.(vname) = tempidx;
                        end
                        [~,~,JR,]= ridgewalk(tempidx2,vars{variables, 3});
                        count = 0;
                        ridgeperiod{1,sim}.(vname) = [];
                        ridgeperiod2{1,sim}.(vname) =[];
                         for ridge2 = 1:length(tempidx)
                            ridgeperiod2{1,sim}.(vname)(ridge2,1) = max(tempidx(ridge2,:));
                            [~,col,~] = find(tempidx(ridge2,:) == ridgeperiod2{1,sim}.(vname)(ridge2,1));  
                            ridgeperiod2{1,sim}.(vname)(ridge2,1) = tau(col,1);
                        end
                        for ridge = 1:length(JR)
                            if isnan(JR(ridge,1)) == true 
                                count = count + 1;
                                ridgeperiod{1,sim}.(vname)(1, count) = JR(ridge,1);
                            elseif isnan(JR(ridge,1)) == false 
                                count = count + 1;
                                ridgeperiod{1,sim}.(vname)(1, count) = tau(JR(ridge,1),1);
                            end
                        end
                        temp = []; temp2 = [] ; tempidx = [] ; JR = [];
                    end
                end
                tempPCs{1,sim}.(lites) = tempPC{1,sim};
                if sim == 24
                    tempNorms{1,sim}.(lites) = tempNorm{1,sim};
                end
                tempRidges{1,sim}.(lites) = ridgeperiod{1,sim};
                tempRidges2{1,sim}.(lites) = ridgeperiod2{1,sim};
                tempPC{1,sim} = []; tempNorm{1,sim} = []; ridgeperiod{1,sim} = []; ridgeperiod2{1,sim} = [];
            end
             PCsreps{rep,sim} = tempPCs{1,sim};
             if sim == 24  
                Normsreps{rep,sim} = tempNorms{1,sim};
             end
             Ridgereps{rep,sim} = tempRidges{1,sim};
             Ridgereps2{rep,sim} = tempRidges2{1,sim};
             tempPCs{1,sim} = []; tempNorms{1,sim} = []; tempRidges{1,sim} = []; tempRidges2{1,sim} = [];
        end      
     end
    filenamer = strcat('/mnt/ide0/share/prendergastlab/', num2str(beta),  filenamer);
    save(filenamer, 'vars', 'Normsreps', 'PCsreps', 'Ridgereps', 'Ridgereps2', '-v7.3');
    clearvars -except idx vars timeseries edge reps numberofsimulations gamma beta files 
end
toc