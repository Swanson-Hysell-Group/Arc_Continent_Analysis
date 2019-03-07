% Macdonald et al., Science, 2019
% Arc-continent collisions in the tropics set Earth�s climate state
% Authors: Francis A. Macdonald, Nicholas L. Swanson-Hysell, Yuem Park, Lorraine Lisiecki, Oliver Jagoutz
% The code calculates statistics comparing ice extent with suture lengths for
% the past 520 Myr
% code written by Lorraine Lisiecki (lisiecki@geol.ucsb.edu)

clear all
close all

alldata=load('./code_output/ice_LIP_suture_lengths.txt');
suture=alldata(:,1:10);
time5=suture(:,1);
ice5=suture(:,2);  % Ice sampled every 5 Myr

%list of RGB colors to use in plotting
col={.5*[1 1 1],[72,21,103]/255,[51 99 141]/255,[32 163 135]/255,[200 200 0]/255};

rng(1) % Seed random numbers so results are reproducible

%First 2 columns of suture file are age and ice extent
ice=suture(:,1:2);
icenorm=ice(:,2)/max(ice(:,2));  % Scale ice extent from 0-1 (for plotting)


icelat_real=flipud(ice(:,2));  %Flip time series to start at 520 Ma
time=-1*ice(:,1);  %Use negative numbers for past times

N=length(time); %105 5-Myr time steps

%Create variables for each suture record, normalized by max/modern length
sut=suture(:,[3:6 10]); % Extract desired columns with suture data
suture_all=sut(:,1)/max(sut(:,1)); % all suture
suture10=sut(:,2)/max(sut(:,2));  % within 10 degree of equator
suture15=sut(:,3)/max(sut(:,3));  % within 15 degree of equator
suture20=sut(:,4)/max(sut(:,4));  % within 20 degree of equator
suture40=sut(:,5)/max(sut(:,5));  % MORE than 40 degree from equator

%Flip time series to start at 520 Ma
flipsutall=flipud(suture_all);
flipsut10=flipud(suture10);
flipsut15=flipud(suture15);
flipsut20=flipud(suture20);
flipsut40=flipud(suture40);

% Identify true(1)/false(0) whether sutures are present at each time step
% Only count sutures that are at least 20% as extensize as modern
binsutall=flipud(suture_all)'>0.2;
binsut10=flipud(suture10)'>0.2;
binsut15=flipud(suture15)'>0.2;
binsut20=flipud(suture20)'>0.2;
binsut40=flipud(suture40)'>0.2;

disp('_________________________________________________________________________________')
disp('Percent of record with sutures >20% of modern: globally, within 10,15,20 degrees, and >40')
100/length(suture)*[sum(binsutall) sum(binsut10) sum(binsut15) sum(binsut20) sum(binsut40)]

figure(1)
plot(time, icenorm,'b','LineWidth',2)
hold on
plot(time, suture_all,'Color',col{1})
plot(time, suture10,'Color',col{2})
plot(time, suture15,'Color',col{3})
plot(time, suture20,'Color',col{4})
plot(time, suture40,'Color',col{5})
hold off
legend('ice latitude','suture (total)','<10^{o}','<15^{o}','<20^{o}','>40^{o}','Location','EastOutside')
title('Suture and Ice extent comparison')
ylabel('Fraction of maximum value')
xlabel('Time (Myr)')
axis tight
grid on


% Correlation calculations for real data

% Correlation coefficients between ice extent and suture lengths
cc=corrcoef(icelat_real,flipsutall);
cc_realall=cc(1,2);
cc=corrcoef(icelat_real,flipsut10);
cc_real10=cc(1,2);
cc=corrcoef(icelat_real,flipsut15);
cc_real15=cc(1,2);
cc=corrcoef(icelat_real,flipsut20);
cc_real20=cc(1,2);
cc=corrcoef(icelat_real,flipsut40);
cc_real40=cc(1,2);

disp('Correlation between ice extent and suture length [Global, <10, <15, <20, >40]')
[cc_realall cc_real10 cc_real15 cc_real20 cc_real40]

%% Create suture records with random age errors
% to evaluate effect of age uncertainty on correlation estimates
% Age errors are simulated as random walks scaled to have a maximum of 5 or 10 Myr

cc_sut=[];
clear ssut %variable for sutures with simulated age uncertainty
for i=1:1000  %Generate ensemble of 1000 suture records with simulated age uncertainty
    rn=randi(3,102,1)-2; %generate numbers -1, 0 or 1
    rn=[0;0;0;rn];  % Assume no age error for 0-10 Ma
    ageunc=cumsum(rn);  % Create random walk
    ageunc=10*ageunc/max(abs(ageunc)); %scale to max error of 10 Myr (can also be set to 5 Myr)
    if ageunc(end)<0
        ageunc(end)=-1*ageunc(end); % make sure last date is >520 Ma
    end
    timeunc = time5 + ageunc; % Add age errors to time

    % Use splines to interpolate age-perturbed records to even 5-Myr sampling
    for j=1:5 % REPEAT FOR DIFFERENT LATITUDE BANDS
        simsut=spline(timeunc,sut(:,j),[0:5:520]'); % spline interpolation of suture length to 5 Myr increments
        ind=find(simsut<0);  % Identify negative suture lengths produced by spline
        simsut(ind)=0;       % Remove negative values
        simsut_flip=flipud(simsut);  % Flip time series to start at 520 Ma
        ssut{j}(:,i)=simsut_flip; % Store suture record perturbed with simulated age uncertainty

        cc=corrcoef(icelat_real,simsut_flip);  % correlation of ice extent and sutures with simulated age errors
        cc_sut(i,j)=cc(1,2);  %store iterations of suture-ice correlation with age uncertainty
    end
end

disp('Average correlation and std dev for age-perturbed sutures [Global, <10, <15, <20, >40]')
[mean(cc_sut); std(cc_sut)]


for j=1:5 % REPEAT FOR DIFFERENT LATITUDE BANDS
    if j>1
        ind=find(cc_sut(:,j)>cc_sut(:,1)); % find how correlations that are better than total suture length
        better_than_tot(j)=length(ind); % store number of better correlations
    end
    cc_sut_sort(:,j)=sort(cc_sut(:,j)); % sort iterations by correlation
    cclow(:,j)=cc_sut_sort(25,j); % find lower band 95 percent confidence interval
    ccup(:,j)=cc_sut_sort(975,j); % find upper band 95 percent confidence interval
end
disp('Fraction of correlations for lat band that are better than global total [<10, <15, <20, >40]')
better_than_tot/1000


%PLOT HISTOGRAMS of ice correlation to suture records with age uncertainty
figure(2)
subplot(511)
histogram(cc_sut(:,1),'FaceColor',col{1}) % total ice-suture correlation with age uncertainty
hold on
plot(cc_realall,1,'ko','MarkerFaceColor',col{1}) %observed correlation
hold off
legend('Total','Location','NorthWest')
xlim([.5 .7])
title('Age uncertainty simulations (10 Myr)')
%title('Correlation of sutures (w/ age uncertainty) & ice (original time samples)')

subplot(512)
histogram(cc_sut(:,2),'FaceColor',col{2}); % tropical (<10) ice-suture correlation with age uncertainty
hold on
plot(cc_real10,1,'ko','MarkerFaceColor',col{2}) %observed correlation
hold off
legend('<10^o','Location','NorthWest')
xlim([.5 .7])

subplot(513)
histogram(cc_sut(:,3),'FaceColor',col{3}); % tropical (<15) ice-suture correlation with age uncertainty
hold on
plot(cc_real15,1,'ko','MarkerFaceColor',col{3})  %observed correlation
hold off
legend('<15^o','Location','NorthWest')
xlim([.5 .7])

subplot(514)
histogram(cc_sut(:,4),'FaceColor',col{4}) % tropical (<20) ice-suture correlation with age uncertainty
hold on
plot(cc_real20,1,'ko','MarkerFaceColor',col{4}) % observed correlation
hold off
legend('<20^o','Location','NorthWest')
xlim([.5 .7])

%Comparison of total sutures, tropical (<20), and high-latitude (>40)
subplot(515)
h1=histogram(cc_sut(:,5),'FaceColor',col{5}); % hiagh-lat (>40) ice-suture correlation with age uncertainty
hold on
h2=histogram(cc_sut(:,1),'FaceColor',col{1}); % total ice-suture correlation with age uncertainty
h3=histogram(cc_sut(:,4),'FaceColor',col{4}); % tropical (<20) ice-suture correlation with age uncertainty
plot(cc_real20,1,'ko','MarkerFaceColor',col{4}) % observed correlation
plot(cc_realall,1,'ko','MarkerFaceColor',col{1}) % observed correlation
plot(cc_real40,1,'ko','MarkerFaceColor',col{5}) % observed correlation
hold off
legend('>40^o','total','<20^o','Location','NorthWest')
xlabel('Pearson Correlation')
xlim([-.3 .7])

% HISTOGRAMS TO COMPARE SUTURES WITH LIP AND ARC LENGTH
figure(3)
h2=histogram(cc_sut(:,1),[.1:.01:.72],'FaceColor',col{1}); % total ice-suture correlation with age uncertainty
hold on
h3=histogram(cc_sut(:,3),[.1:.01:.72],'FaceColor',col{3}); % tropical (<15) ice-suture correlation with age uncertainty
plot(cc_realall,1,'ko','MarkerFaceColor',col{1})
plot(cc_real15,1,'ko','MarkerFaceColor',col{3})
hold off
xlabel('Pearson Correlation')
ylabel('Frequency')
title('Effect of simulated 10 Myr age uncertainty')
xlim([-.3 .7])


%%  Identify glaciated and unglaciated intervals (OVERLAP CALCULATIONS)


% Identify true(1)/false(0) whether ice is present at each time step
% Only count ice extent at least 10 degrees from poles
binice=flipud(ice(:,2))'>10;
disp(['Percent of record with ice extent > 10^o is ' num2str(100*sum(binice)/length(ice))])
disp('')

% Identify/measure durations of glaciated and unglaciated intervals
% CAUTION: Code assumes that there are 4 identified glaciations

tempbinice=binice;  % true(1)/false(0) whether ice is present
tempice=flipud(ice(:,2))';  % latitude extent of ice
for i=1:4
    ind=min(find(tempbinice==1)); % Find first instance of ice
    no_ice_dur(i)=ind-1;          % Duration of first unglaciated interval
    tempbinice=tempbinice(ind:end);  %Remove identified portion from temp record
    tempice=tempice(ind:end);        %Remove identified portion from temp record
    if i<4
        ind=min(find(tempbinice==0));  % Find next time without ice
        ice_lat{i}=tempice(1:ind-1);   % ice extent within glaciated interval
        ice_dur(i)=ind-1;              % duration of glaciated interval
        tempbinice=tempbinice(ind:end); %Remove identified portion from temp record
        tempice=tempice(ind:end);
    else
        ice_dur(i)=length(tempbinice);  % duration of last glaciated interval
        ice_lat{i}=tempice(1:end);      % ice extent within last glaciated interval
    end
end
ice_lat_real=ice_lat;

disp('Number of 5-Myr time steps in each ice age (1st line) & between ice ages (2nd line)')
[ice_dur; no_ice_dur]

Ngap=sum(no_ice_dur);  % Total time steps without ice
pctice=sum(ice_dur)/105;  %Percent of time with ice present

% Calculate overlap between presence of ice and sutures
ovlpall=sum(binice.*binsutall);
ovlp10=sum(binice.*binsut10);
ovlp15=sum(binice.*binsut15);
ovlp20=sum(binice.*binsut20);
ovlp40=sum(binice.*binsut40);
disp('Ice-suture overlap as percent of time glaciated (ice extent > 10 degrees latitude)')
disp('[Global, <10, <15, <20, >40]')
[100/sum(binice)*[ovlpall ovlp10 ovlp15 ovlp20 ovlp40]]
disp('Ice-suture overlap as percent of time sutures >20% of modern [Global, <10, <15, <20, >40]')
[100*[ovlpall/sum(binsutall) ovlp10/sum(binsut10) ovlp15/sum(binsut15)...
    ovlp20/sum(binsut20) ovlp40/sum(binsut40)]]




%% Create simulated climate records by re-arranging ice intervals

% Simulations will be used to test the hypothesis that ice ages are
% NOT related to suture length. Therefore, climate simulations are
% designed to generate uniform probability for 4 ice ages occuring
% at anytime throughout the last 520 Myr

% clear memory for variables created within loop
all_simice=[]; all_gap=[]; all_end=[]; all_simicelat=[];
ovlp_sim=[]; cc_sim=[]; allstart=[];

ITER=10000; % Number of Monte carlo simulations (use 10,000 or 20,000)

for j=1:ITER
    ice_lat=ice_lat_real;  % create copy of ice latitude extent
    iceflip=randi(2,1,4);  % 50% chance of mirror-image for lat extent variation
    for i=1:4
        if iceflip(i)==1
            ice_lat{i}=fliplr(ice_lat_real{i});
        end
    end

    iceorder=randperm(4);  % pick random order for ice intervals

    % clear variables used to create climate simulation
    simi=[]; sim=[]; simlat=[];

    % Percent of time with ice determines chance of ice at start or end of record
    endstate=(rand(1,2) < pctice);
    all_end(j,:)=endstate; %save endstate results

    if sum(endstate)>=1 % If record starts OR ends in a glacial state

        % If starting with ice, 40 percent chance simulation starts midway through an ice age
        % In which case ice age wraps around (i.e., one portion at start of
        % simulation and the rest at end of simulation)
        if rand <.4 && ice_dur(iceorder(1))>1
            startdur=randi(ice_dur(iceorder(1))-1);  % Pick which time step in the glaciated interval

        % If not starting midway through glaciation, 50 percent chance ice age is at start of simulation
        elseif rand < .5
            startdur=0;  %starts in non-glacial state and ends with ice
        else
            startdur=ice_dur(iceorder(1)); %starts with full duration of an ice age
        end

        % Begin constructing climate simulation

        if startdur>0  % If simulation starts with an ice age, apend duration/latitude of ice
            simi=[simi startdur];
            sim=[sim ones(1,startdur)];
            simlat=[simlat ice_lat{iceorder(1)}(end-startdur+1:end)];
        end

        allstart(end+1,:)=[iceorder(1) startdur]; %save info about initial ice age

        rn=lognrnd(0,1,1,4);  %Generate log-normal wait times between ice ages
        rand_gap=round(Ngap*rn/sum(rn));  %Scale wait time sum to total non-glaciated timespan

        % Small adjustment if rounding doesn't produce correct sum
        rand_gap(2)=rand_gap(2)-(sum(rand_gap)-Ngap);
        while min(rand_gap)<=0 % repeat if adjustment produces a negative wait time
            rn=lognrnd(0,1,1,4);
            rand_gap=round(Ngap*rn/sum(rn));
            rand_gap(2)=rand_gap(2)-(sum(rand_gap)-Ngap);
        end

        for i=1:3  % Append subsequent non-glacial and glacial intervals
            simi=[simi rand_gap(i) ice_dur(iceorder(i+1))];
            sim=[sim zeros(1,rand_gap(i)) ones(1,ice_dur(iceorder(i+1)))];
            simlat=[simlat zeros(1,rand_gap(i)) ice_lat{iceorder(i+1)}];
        end

        % If starting with non-glacial state or starting midway through glacial interval
        % add glacial state to end of simulation
        %(remaining portion if starting midway through)
        if startdur < ice_dur(iceorder(1))
            simi=[simi rand_gap(4) ice_dur(iceorder(1))-startdur];
            sim=[sim zeros(1,rand_gap(4)) ones(1,ice_dur(iceorder(1))-startdur)];
            simlat=[simlat zeros(1,rand_gap(4)) ice_lat{iceorder(1)}(1:end-startdur)];

        else % otherwise end in non-glacial state
            simi=[simi rand_gap(4)];
            sim=[sim zeros(1,rand_gap(4))];
            simlat=[simlat zeros(1,rand_gap(4))];
        end


    else  % If simulation starts AND ends in non-glacial state

        rn=lognrnd(0,1,1,5); % create log-normal wait times for 5 non-glacial intervals
        rand_gap=round(Ngap*rn/sum(rn));
        % Small adjustment if rounding doesn't produce correct sum
        rand_gap(2)=rand_gap(2)-(sum(rand_gap)-Ngap);
        while min(rand_gap)<=0 % repeat if adjustment produces a negative wait time
            rn=lognrnd(0,1,1,5);
            rand_gap=round(Ngap*rn/sum(rn));
            rand_gap(2)=rand_gap(2)-(sum(rand_gap)-Ngap);
        end

        % Append successive non-glacial and glacial intervals
        for i=1:4
            simi=[simi rand_gap(i) ice_dur(iceorder(i))];
            sim=[sim zeros(1,rand_gap(i)) ones(1,ice_dur(iceorder(i)))];
            simlat=[simlat zeros(1,rand_gap(i)) ice_lat{iceorder(i)}];
        end
        % add last non-glacial state
        simi(end+1)=rand_gap(5);
        sim=[sim zeros(1,rand_gap(5))];
        simlat=[simlat zeros(1,rand_gap(5))];
    end
    % end construction of climate simulation

    % check that simulation has correct total length or print an error
    if length(simlat)~=105
        disp(['Error simulation length =' num2str(sum(simi))])
        disp(['Error simulation length =' num2str(length(simlat))])
        sim=sim(1:105);
    end

    % add simulation to matrix containing all simulations
    simice=sim;
    all_simice(j,:)=simice;
    all_simicelat(j,:)=simlat;

    % calculate and store overlap between simulated ice and real suture data
    ovlp_simall(j)=sum(simice.*binsutall);
    ovlp_sim10(j)=sum(simice.*binsut10);
    ovlp_sim15(j)=sum(simice.*binsut15);
    ovlp_sim20(j)=sum(simice.*binsut20);
    ovlp_sim40(j)=sum(simice.*binsut40);


    % calculate correlation between simulated ice and real suture data
    cc=corrcoef(simlat,flipsutall);
    cc_simall(j)=cc(1,2);
    cc=corrcoef(simlat,flipsut10);
    cc_sim10(j)=cc(1,2);
    cc=corrcoef(simlat,flipsut15);
    cc_sim15(j)=cc(1,2);
    cc=corrcoef(simlat,flipsut20);
    cc_sim20(j)=cc(1,2);
    cc=corrcoef(simlat,flipsut40);
    cc_sim40(j)=cc(1,2);
end
% End loop that creates all climate simulations


%% Calculate and print p-values for real result compared to simulations
% p < 0.05 means we reject the null hypothesis that ice is unrelated
% to a particular suture time series.

p_ovlpall=length(find(ovlp_simall>=ovlpall))/ITER;
p_ovlp10=length(find(ovlp_sim10>=ovlp10))/ITER;
p_ovlp15=length(find(ovlp_sim15>=ovlp15))/ITER;
p_ovlp20=length(find(ovlp_sim20>=ovlp20))/ITER;
p_ovlp40=length(find(ovlp_sim40>=ovlp40))/ITER;
disp('Null hypothesis test for random timing of glaciations')
disp('p values: overlap [Global, <10, <15, <20, >40]')
[p_ovlpall p_ovlp10 p_ovlp15 p_ovlp20 p_ovlp40]

p_corrall=length(find(cc_simall>=cc_realall))/ITER;
p_corr10=length(find(cc_sim10>=cc_real10))/ITER;
p_corr15=length(find(cc_sim15>=cc_real15))/ITER;
p_corr20=length(find(cc_sim20>=cc_real20))/ITER;
p_corr40=length(find(cc_sim40>=cc_real40))/ITER;
disp('p values: correlation [Global, <10, <15, <20, >40]')
[p_corrall p_corr10 p_corr15 p_corr20 p_corr40]


% Create figures to show results
figure(4)
% % Summary statistics of ice simulations
% % Flat lines of ~1 are best fit with evenly distributed probabilities
% % for ice through time (ie, null hypothesis)
% subplot(311)
% plot(sum(all_simice)/(sum(ice_dur)/105*ITER))
% hold on
% plot(sum(all_simicelat)/(mean(icelat_real)*ITER))
% hold off
% axis tight
% title(['Ice lat > 10, ITER=' num2str(ITER)])
% xlabel('Time')
% ylabel('Fraction of Expected')
% legend('Presence','Average Latitude')

% Colors in each horizontal line represent ice latitude vs time for each
% simulation (should look like random noise)
%subplot('Position',[.09 .05 .85 .58])
imagesc(all_simicelat)
colorbar
title('Simulated latitude of ice extent')
xlabel('Time (in 5-Myr steps)')
ylabel('Iteration')

%%
% PLOT OF NULL HYPOTHESIS TEST P VALUES
% Histograms illustrating p-values for suture overlaps and correlations

figure(5)
% OVERLAP HISTOGRAMS
% Display overlap as a percent of amount of time sutures are present
% because, e.g., 5-degree sutures are present/extensive for less total time
% than 10-degree sutures
subplot(5,2,1)
histogram(100*ovlp_simall/sum(binsutall),[0:5:100],'FaceColor',col{1})
hold on
plot(100*ovlpall/sum(binsutall),1,'ko','MarkerFaceColor',col{1})
hold off
axis tight
text(60,1000,['p = ' num2str(p_ovlpall)])
legend('total')
title('Overlap (Ice extent > 10^o, Sutures > 20% of max)')

subplot(5,2,3)
histogram(100*ovlp_sim10/sum(binsut10),[0:5:100],'FaceColor',col{2})
hold on
plot(100*ovlp10/sum(binsut10),1,'ko','MarkerFaceColor',col{2})
text(60,800,['p = ' num2str(p_ovlp10)])
hold off
axis tight
legend('<10')

subplot(5,2,5)
histogram(100*ovlp_sim15/sum(binsut15),[0:5:100],'FaceColor',col{3})
hold on
plot(100*ovlp15/sum(binsut15),1,'ko','MarkerFaceColor',col{3})
text(60,1000,['p = ' num2str(p_ovlp15)])
hold off
axis tight
legend('<15')

subplot(5,2,7)
histogram(100*ovlp_sim20/sum(binsut20),[0:5:100],'FaceColor',col{4})
hold on
plot(100*ovlp20/sum(binsut20),1,'ko','MarkerFaceColor',col{4})
text(60,800,['p = ' num2str(p_ovlp20)])
hold off
axis tight
legend('<20')

subplot(5,2,9)
histogram(100*ovlp_sim40/sum(binsut40),[0:5:100],'FaceColor',col{5})
hold on
plot(100*ovlp40/sum(binsut40),1,'ko','MarkerFaceColor',col{5})
text(60,1000,['p = ' num2str(p_ovlp40)])
hold off
axis tight
xlabel('%Time glaciated')
legend('>40')

% CORRELATION HISTOGRAMS
subplot(5,2,2)
histogram(cc_simall,'FaceColor',col{1})
hold on
plot(cc_realall,1,'ko','MarkerFaceColor',col{1})
text(.4,350,['p = ' num2str(p_corrall)])
hold off
axis tight
legend('total')
xlim([-.6 .8])
title('Correlation (bars = null hyp., o = data)')

subplot(5,2,4)
histogram(cc_sim10,'FaceColor',col{2})
hold on
plot(cc_real10,1,'ko','MarkerFaceColor',col{2})
hold off
axis tight
text(.4,350,['p = ' num2str(p_corr10)])
legend('<10')
xlim([-.6 .8])

subplot(5,2,6)
histogram(cc_sim15,'FaceColor',col{3})
hold on
plot(cc_real15,1,'ko','MarkerFaceColor',col{3})
hold off
axis tight
text(.4,300,['p = ' num2str(p_corr15)])
legend('<15')
xlim([-.6 .8])

subplot(5,2,8)
histogram(cc_sim20,'FaceColor',col{4})
hold on
plot(cc_real20,1,'ko','MarkerFaceColor',col{4})
hold off
axis tight
text(.4,300,['p = ' num2str(p_corr20)])
legend('<20')
xlim([-.6 .8])

subplot(5,2,10)
histogram(cc_sim40,'FaceColor',col{5})
hold on
plot(cc_real40,1,'ko','MarkerFaceColor',col{5})
hold off
axis tight
xlabel('Correlation')
text(.4,300,['p = ' num2str(p_corr40)])
legend('>40')
xlim([-.6 .8])

%%
%*****************************************************
%             LIP and ARC STATS
%*****************************************************

rng(1) % Seed random numbers so results are reproducible

liparc=alldata(:,[11:17]);

%Create variables for each LIP record
LIP_allD=liparc(:,1)/max(liparc(:,1)); % total_LIP_decay
LIP_allDB=liparc(:,2)/max(liparc(:,2));  % total_LIP_decay_burial
LIP_15D=liparc(:,3)/max(liparc(:,3));  % within_15_LIP_decay
LIP_15DB=liparc(:,4)/max(liparc(:,4));  % within_15_LIP_decay_burial
LIP_hilatD=liparc(:,5)/max(liparc(:,5));  % outside_15_LIP_decay
LIP_hilatDB=liparc(:,6)/max(liparc(:,6));  % outside_15_LIP_decay_burial

arc=alldata(:,17);
arcn=-1*(arc-max(arc))/(max(arc)-min(arc)); % multiply by -1 due to expected negative correlation, scale to a range of 0-1
liparc(:,7)=arcn;

%Flip time series to start at 520 Ma
fliplipalld=flipud(LIP_allD);
fliplipalldb=flipud(LIP_allDB);
fliplip15d=flipud(LIP_15D);
fliplip15db=flipud(LIP_15DB);
fliplip_hilatd=flipud(LIP_hilatD);
fliplip_hilatdb=flipud(LIP_hilatDB);
fliparc=flipud(arcn);

% Identify true(1)/false(0) whether LIPs/arcs are present at each time step
% Only count LIP that are at least 20% as extensize as max
binlipalld=flipud(LIP_allD)'>0.2;
binlipalldb=flipud(LIP_allDB)'>0.2;
binlip15d=flipud(LIP_15D)'>0.2;
binlip15db=flipud(LIP_15DB)'>0.2;
binlip_hilatd=flipud(LIP_hilatD)'>0.2;
binlip_hilatdb=flipud(LIP_hilatDB)'>0.2;
% Only count arcs that are more than 70% of range (max scaled to 0, min rescaled to 1)
binarc=flipud(arcn)'>0.7;

disp('_________________________________________________________________________________')
disp('Percent of record with LIP >20% of max: total (d, d+b), <15 (d, d+b), >15 (d, d+b), arc')
100/length(liparc)*[sum(binlipalld) sum(binlipalldb) sum(binlip15d) sum(binlip15db) sum(binlip_hilatd) sum(binlip_hilatdb) sum(binarc)]

figure(6)
plot(time, icenorm,'LineWidth',2)
hold on
plot(time, LIP_allD)
plot(time, LIP_allDB)
plot(time, LIP_15D)
plot(time, LIP_15DB)
plot(time, LIP_hilatD,'--')
plot(time, LIP_hilatDB,':')
plot(time, arcn,'-.')
hold off
legend('ice','LIP decay','LIP d+b','LIP <15 s','LIP <15 d+b','LIP >15 d','LIP >15 d+b','-1*arc','Location','EastOutside')
ylabel('Fraction of maximum value (or range)')
xlabel('Time (Myr)')
axis tight
grid on


% Correlation calculations for real data

% Correlation coefficients between ice extent and LIP lengths
cc=corrcoef(icelat_real,fliplipalld);
cc_real_lipalld=cc(1,2);
cc=corrcoef(icelat_real,fliplipalldb);
cc_real_lipalldb=cc(1,2);
cc=corrcoef(icelat_real,fliplip15d);
cc_real_lip15d=cc(1,2);
cc=corrcoef(icelat_real,fliplip15db);
cc_real_lip15db=cc(1,2);
cc=corrcoef(icelat_real,fliplip_hilatd);
cc_real_liphid=cc(1,2);
cc=corrcoef(icelat_real,fliplip_hilatdb);
cc_real_liphidb=cc(1,2);
cc=corrcoef(icelat_real,fliparc);
cc_realarc=cc(1,2);

disp('Correlation between ice extent and LIP length [Global (d, d+b), <15 (d, d+b), >15 (d, d+b), arc')
[cc_real_lipalld cc_real_lipalldb cc_real_lip15d cc_real_lip15db cc_real_liphid cc_real_liphidb cc_realarc]

% Calculate overlap between presence ice and LIPs
ovlp_lipalld=sum(binice.*binlipalld);
ovlp_lipalldb=sum(binice.*binlipalldb);
ovlp_lip15d=sum(binice.*binlip15d);
ovlp_lip15db=sum(binice.*binlip15db);
ovlp_liphid=sum(binice.*binlip_hilatd);
ovlp_liphidb=sum(binice.*binlip_hilatdb);
ovlparc=sum(binice.*binarc);

disp('Ice-prediction overlap as percent of time glaciated (ice>10) [Global (d, d+b), <15 (d, d+b) and >15 (d, d+b), -1*arc]')
100/sum(binice)*[ovlp_lipalld ovlp_lipalldb ovlp_lip15d ovlp_lip15db ovlp_liphid ovlp_liphidb ovlparc]
disp('Ice-prediction overlap as percent of time predicted by LIP/arc [Global (d, d+b), <15 (d, d+b) and >15 (d, d+b), -1*arc]')
100*[ovlp_lipalld/sum(binlipalld) ovlp_lipalldb/sum(binlipalldb) ovlp_lip15d/sum(binlip15d)...
    ovlp_lip15db/sum(binlip15db) ovlp_liphid/sum(binlip_hilatd) ovlp_liphidb/sum(binlip_hilatdb) ovlparc/sum(binarc)]


%% Create LIP records with random age errors
% to evaluate effect of age uncertainty on correlation estimates
% Age errors are simulated as random walks scaled to have a maximum of 10 Myr

cc_liparc=[];
clear sliparc %simulated LIP or arc record
for i=1:1000
    rn=randi(3,102,1)-2; %generate numbers -1, 0 or 1
    rn=[0;0;0;rn];  % Assume no age error for 0-10 Ma
    ageunc=cumsum(rn);  % Create random walk
    ageunc=10*ageunc/max(abs(ageunc)); %scale to max error of 10 Myr
    if ageunc(end)<0
        ageunc(end)=-1*ageunc(end); % make sure last date is >520 Ma
    end
    timeunc=time5+ageunc; %Add age errors to time

    % Use splines to interpolate age-perturbed records to even 5-Myr sampling
    for j=1:7
        simliparc=spline(timeunc,liparc(:,j),[0:5:520]'); %spline interpolation
        ind=find(simliparc<0);  % Identify negative LIP lengths produced by spline
        simliparc(ind)=0;       % Remove negative values
        simliparc_flip=flipud(simliparc);  %Flip time series to start at 520 Ma
        sliparc{j}(:,i)=simliparc_flip;

        cc=corrcoef(icelat_real,simliparc_flip);
        cc_liparc(i,j)=cc(1,2);
    end
end

for j=1:7
    cc_liparc_sort(:,j)=sort(cc_liparc(:,j));
    cclow(:,j)=cc_liparc_sort(25,j);
    ccup(:,j)=cc_liparc_sort(975,j);
end

%[Global (d, d+b), <15 (d, d+b) and >15 (d+b)]

figure(7)
subplot(311)
histogram(cc_liparc(:,1))
title('Correlation of LIP (w/ age uncertainty) & ice (original time samples)')
hold on
histogram(cc_liparc(:,2))
histogram(cc_liparc(:,7))
plot(cc_real_lipalld,1,'ko','MarkerFaceColor','b')
plot(cc_real_lipalldb,1,'ko','MarkerFaceColor','r')
plot(cc_realarc,1,'ko','MarkerFaceColor','y')
hold off
legend('total d','total d+b','-1*arc','Location','NorthEast')
xlim([-.25 .5])

subplot(312)
histogram(cc_liparc(:,3))
hold on
histogram(cc_liparc(:,4))
plot(cc_real_lip15d,1,'ko','MarkerFaceColor','b')
plot(cc_real_lip15db,1,'ko','MarkerFaceColor','r')
hold off
legend('<15 (d)','<15 (d+b)','Location','NorthEast')
xlim([-.25 .5])

subplot(313)
histogram(cc_liparc(:,5))
hold on
histogram(cc_liparc(:,6))
plot(cc_real_liphid,1,'ko','MarkerFaceColor','b')
plot(cc_real_liphidb,1,'ko','MarkerFaceColor','r')
hold off
legend('>15 (d)','>15 (d+b)','Location','NorthEast')
xlabel('Pearson Correlation')
xlim([-.25 .5])

% HISTOGRAMS TO COMPARE SUTURES WITH LIP AND ARC LENGTH
figure(3)
hold on
h4=histogram(cc_liparc(:,4),[0:.01:.72],'FaceColor',[1 .4 0]);
h5=histogram(cc_liparc(:,7),[0:.01:.72],'FaceColor',[.4 0 0]);
plot(cc_real_lip15db,1,'ko','MarkerFaceColor',[1 .4 0])
plot(cc_realarc,1,'ko','MarkerFaceColor',[.5 0 0])
hold off
legend([h2 h3 h4 h5],'suture total','suture <15^o','LIP <15^o (d+b)','-1*arc','Location','NorthWest')
xlim([0 .7])

disp('Average correlation and std dev for age-perturbed LIP [Global (d, d+b), <15 (d, d+b), >15 (d, d+b), -1*arc]')
[mean(cc_liparc); std(cc_liparc)]

%% NULL HYPOTHESIS TESTS for simulated random glaciations

% Loop through all previously generated climate simulations
for j=1:ITER
    % calculate overlap between simulated ice and real LIP data
    ovlp_simlipalld(j)=sum(all_simice(j,:).*binlipalld);
    ovlp_simlipalldb(j)=sum(all_simice(j,:).*binlipalldb);
    ovlp_simlip15d(j)=sum(all_simice(j,:).*binlip15d);
    ovlp_simlip15db(j)=sum(all_simice(j,:).*binlip15db);
    ovlp_simliphid(j)=sum(all_simice(j,:).*binlip_hilatd);
    ovlp_simliphidb(j)=sum(all_simice(j,:).*binlip_hilatdb);
    ovlp_simarc(j)=sum(all_simice(j,:).*binarc);


    % calculate correlation between simulated ice and real LIP data
    cc=corrcoef(all_simicelat(j,:),fliplipalld);
    cc_simlipalld(j)=cc(1,2);
    cc=corrcoef(all_simicelat(j,:),fliplipalldb);
    cc_simlipalldb(j)=cc(1,2);
    cc=corrcoef(all_simicelat(j,:),fliplip15d);
    cc_simlip15d(j)=cc(1,2);
    cc=corrcoef(all_simicelat(j,:),fliplip15db);
    cc_simlip15db(j)=cc(1,2);
    cc=corrcoef(all_simicelat(j,:),fliplip_hilatd);
    cc_simliphid(j)=cc(1,2);
    cc=corrcoef(all_simicelat(j,:),fliplip_hilatdb);
    cc_simliphidb(j)=cc(1,2);
    cc=corrcoef(all_simicelat(j,:),fliparc);
    cc_simarc(j)=cc(1,2);
end


%% Calculate and print p-values for real result compared to simulations
% p < 0.05 means we reject the null hypothesis that ice is unrelated
% to a particular LIP time series.

p_ovlp_lipalld=length(find(ovlp_simlipalld>=ovlp_lipalld))/ITER;
p_ovlp_lipalldb=length(find(ovlp_simlipalldb>=ovlp_lipalldb))/ITER;
p_ovlp_lip15d=length(find(ovlp_simlip15d>=ovlp_lip15d))/ITER;
p_ovlp_lip15db=length(find(ovlp_simlip15db>=ovlp_lip15db))/ITER;
p_ovlp_liphid=length(find(ovlp_simliphid>=ovlp_liphid))/ITER;
p_ovlp_liphidb=length(find(ovlp_simliphidb>=ovlp_liphidb))/ITER;
p_ovlparc=length(find(ovlp_simarc>=ovlparc))/ITER;
disp('p values: overlap [Global (d, d+b), <15 (d, d+b), >15 (d, d+b), -1*arc]')
[p_ovlp_lipalld p_ovlp_lipalldb p_ovlp_lip15d p_ovlp_lip15db p_ovlp_liphid p_ovlp_liphidb p_ovlparc]

p_corr_lipalld=length(find(cc_simlipalld>=cc_real_lipalld))/ITER;
p_corr_lipalldb=length(find(cc_simlipalldb>=cc_real_lipalldb))/ITER;
p_corr_lip15d=length(find(cc_simlip15d>=cc_real_lip15d))/ITER;
p_corr_lip15db=length(find(cc_simlip15db>=cc_real_lip15db))/ITER;
p_corr_liphid=length(find(cc_simliphidb>=cc_real_liphid))/ITER;
p_corr_liphidb=length(find(cc_simliphidb>=cc_real_liphidb))/ITER;
p_corrarc=length(find(cc_simarc>=cc_realarc))/ITER;
disp('p values: correlation [Global (d, d+b), <15 (d, d+b), >15 (d, d+b), -1*arc]')
[p_corr_lipalld p_corr_lipalldb p_corr_lip15d p_corr_lip15db p_corr_liphid p_corr_liphidb p_corrarc]


% Histograms illustrating p-values for overlaps and correlations
figure(8)
% Display overlap as a percent of amount of time LIPs are present
% because, e.g., 5-degree LIPs are present/extensive for less total time
% than 10-degree LIPs
subplot(421)
histogram(100*ovlp_simlipalld/sum(binlipalld),[0:5:100])
hold on
histogram(100*ovlp_simlipalldb/sum(binlipalldb),[0:5:100])
plot(100*ovlp_lipalld/sum(binlipalld),1,'ko','MarkerFaceColor','b')
plot(100*ovlp_lipalldb/sum(binlipalldb),1,'ko','MarkerFaceColor','r')
hold off
axis tight
%title('5-Myr ice sample, LIPs > 20% of modern')
title('Overlap: LIP > 20% max or Arc <30% range')
xlabel(['%Time glaciated: p_{d}=' num2str(p_ovlp_lipalld) ', p_{db}=' num2str(p_ovlp_lipalldb)])
legend('LIP total (d)','LIP total (d+b)')


subplot(423)
histogram(100*ovlp_simlip15d/sum(binlip15d),[0:5:100])
hold on
histogram(100*ovlp_simlip15db/sum(binlip15db),[0:5:100])
plot(100*ovlp_lip15d/sum(binlip15d),1,'ko','MarkerFaceColor','b')
plot(100*ovlp_lip15db/sum(binlip15db),1,'ko','MarkerFaceColor','r')
hold off
axis tight
xlabel(['%Time glaciated: p_{d}=' num2str(p_ovlp_lip15d) ', p_{db}=' num2str(p_ovlp_lip15db)])
legend('LIP <15^o (d)','LIP <15^o (d+b)')

subplot(425)
histogram(100*ovlp_simliphid/sum(binlip_hilatd),[0:5:100])
hold on
histogram(100*ovlp_simliphidb/sum(binlip_hilatdb),[0:5:100])
plot(100*ovlp_liphid/sum(binlip_hilatd),1,'ko','MarkerFaceColor','b')
plot(100*ovlp_liphidb/sum(binlip_hilatdb),1,'ko','MarkerFaceColor','r')
hold off
axis tight
xlabel(['%Time glaciated: p_{d}=' num2str(p_ovlp_liphid) ', p_{db}=' num2str(p_ovlp_liphidb)])
legend('LIP >15^o (d)','LIP >15^o (d+b)')

subplot(427)
histogram(100*ovlp_simarc/sum(binarc),[0:5:100])
hold on
plot(100*ovlparc/sum(binarc),1,'ko','MarkerFaceColor','b')
hold off
axis tight
xlabel(['%Time glaciated: p_{arc}=' num2str(p_ovlparc)])
legend('arc')



subplot(422)
histogram(cc_simlipalld)
hold on
histogram(cc_simlipalldb)
plot(cc_real_lipalld,1,'ko','MarkerFaceColor','b')
plot(cc_real_lipalldb,1,'ko','MarkerFaceColor','r')
hold off
axis tight
title('Histograms = simulations, Circles = data')
xlabel(['Correlation: p_{d}=' num2str(p_corr_lipalld) ', p_{db}=' num2str(p_corr_lipalldb)])
legend('LIP total (d)','LIP total (d+b)')
xlim([-.6 .8])

subplot(424)
histogram(cc_simlip15d)
hold on
histogram(cc_simlip15db)
plot(cc_real_lip15d,1,'ko','MarkerFaceColor','b')
plot(cc_real_lip15db,1,'ko','MarkerFaceColor','r')
hold off
axis tight
xlabel(['Correlation: p_{d}=' num2str(p_corr_lip15d) ', p_{db}=' num2str(p_corr_lip15db)])
legend('LIP <15^o (d)','LIP <15^o (d+b)')
xlim([-.6 .8])

subplot(426)
histogram(cc_simliphid)
hold on
histogram(cc_simliphidb)
plot(cc_real_liphid,1,'ko','MarkerFaceColor','b')
plot(cc_real_liphidb,1,'ko','MarkerFaceColor','r')
hold off
axis tight
xlabel(['Correlation: p_{d}=' num2str(p_corr_liphid) ', p_{db}=' num2str(p_corr_liphidb)])
legend('LIP >15^o (d)','LIP >15^o (d+b)')
xlim([-.6 .8])

subplot(428)
histogram(cc_simarc)
hold on
plot(cc_realarc,1,'ko','MarkerFaceColor','b')
hold off
axis tight
xlabel(['Correlation: p_{arc}=' num2str(p_corrarc)])
legend('-1*arc')
xlim([-.6 .8])
