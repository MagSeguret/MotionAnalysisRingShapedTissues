% MESSAGE TO USERS:
% This code works on .mp4 videos containing one single ring each.
% It was written on MacOS, for Windows the "/" may have to be changed to
% "\" (l.493)

% It requires the user to first select the folder of where the videos to
% analyze are located and it will analyze all the videos of the folder.

% To start analyzing not on the first video, change the number of the
% starting video at the start of the "for" loop on videos.

%For each video the sensitivity will be asked (number between 0 and 1 which
%depends on the brightness and contrast of the video), the user must choose
%it to be able to see the inner circle of the ring clearly. In the article we use a sensotivity around 0.3. 
%The user will be asked if he is satistied with the sensitivity, and will answer "n" for no or "y" for
%yes. In the "n" case, the user can enter a new sensitivity. In the "y"
%case, the video will be analyzed.
% If the videos can all be analyzed with the same sensitivity, the
% sensitivity can be written directly in the code.

% At the end of the analysis, an Excel file containing the parameters of
% the differents videos and the variables used for the calculation are all
% saved in a .mat file which can be re-opened to add some parameters or retrieve the contraction traces.


clear all, close all

% Name of the Excel file containing the analyzed parameters
ExcelName='ExportVal_Tissues_Motion_analysis';
FinalExcelName=strcat(ExcelName,'.xlsx');


% Tissue recordings location
disp('Please select the movies directory');
rootPath = uigetdir(path, 'Please select the movies directory');
SavePath=rootPath; %analysis will be saved in the folder of the videos, can be changed.

cd(rootPath)
mkdir motion_analysis
                                                                                                                        
videoPaths=struct();
videoPaths=dir('*.mp4');


millimeterPerPixel=1/910.2; %depends on the camera and the objective used
diameterMM=0.360; % theoretical diameter of the relaxed ring in mm: 0.360mm



for video =1:length(videoPaths)%loop on videos, 1 can be changed to the number of the video on which you want to start your analysis
   
        cd(rootPath)
        
        clear V
        clear AreaFrames
        clear RawData
        clear ExportVal
    


            fileName = videoPaths(video).name;
            fprintf(' >>> analyzing video:  %s \n', fileName)
            cap = [];
            cap = VideoReader(fileName); 


            h = cap.Height;
            w = cap.Width;
            n = cap.NumFrames;

            V = (zeros(h,w,n));
            vidFrame_g = [];

            for i = 1:n
                vidFrame = read(cap,i);
                vidFrame_g(:,:) = rgb2gray(vidFrame);
                tif_file = fullfile(rootPath,strcat(fileName,'_TimeLapse.tif'));
                V(:,:,i) = double(vidFrame_g)/255;
            end
            
            
            %BWsensitivity is between 0 and 1 and depends on the camera,
            %and brightness of the video
            
            %%%%%%%%%%%%
            str='n'; 
            while str=='n'

                promptParam = 'Please enter BWsensitivity: ';
                BWsensitivity = input(promptParam);
                

                pixeldots=2000;

                BW = im2bw(V(:,:,100),BWsensitivity);
                % imshow(BW)
                BW2 = bwareaopen(BW,pixeldots);
                imshow(BW2)
                prompt = 'Are you happy with the analysis? y/n [y]: ';
                str = input(prompt,'s');
                
            end
            
            %%%%%%%%%%%%
            
            % If contrast and brightness is the same between the videos of
            % the folder the part between %%%%%%%%%%%% can be removed and replaced with:
            
            % BWsensitivity = SENSIIVITY;
            % pixeldots=2000;
            % BW2 = bwareaopen(BW,pixeldots);
            % imshow(BW2)
            
            
            diameterPix=diameterMM/millimeterPerPixel;% theoretical diameter in pixels
            
            
            RawData=struct();
            RawData.Name=cap.Name;
            RawData.Duration=cap.Duration;
            RawData.FrameRate=cap.FrameRate;
            RawData.NumFrames=cap.NumFrames;
            RawData.Height=cap.Height;
            RawData.Width=cap.Width;

            %Numbers of frames to remove at the end
            Remove=100;

            DeltaT=1/RawData.FrameRate;
            AnalysisFrames=size(V,3)-Remove;
            AnalysisTime=RawData.Duration-DeltaT*Remove;
            time=[0:DeltaT:AnalysisTime];

            AreaFrames=zeros(AnalysisFrames,1);
            MinorLengthFrames=zeros(AnalysisFrames,1);
            MajorLengthFrames=zeros(AnalysisFrames,1);


            AnalysisTime=size(V,3)-Remove;

            for i=1:AnalysisTime
                BW = im2bw(V(:,:,i+Remove),BWsensitivity);% to remove at the end of the film im2bw(V(:,:,i),BWsensitivity) - to remove at the begining im2bw(V(:,:,i+Remove),BWsensitivity)
                BW2 = bwareaopen(BW,pixeldots);

                stats = regionprops('table',BW2,'Centroid',...
                'MajorAxisLength','MinorAxisLength','Eccentricity','Area');
            
                % Delete regions which are not the expected circle
                ToDelete1=stats.Eccentricity>0.6;%region wanted should be round
                ToDelete2=mean([stats.MajorAxisLength stats.MinorAxisLength],2)<0.75*diameterPix; 
                ToDelete3=mean([stats.MajorAxisLength stats.MinorAxisLength],2)>1.5*diameterPix;
                ToDeleteFinal=max(ToDelete1,max(ToDelete2,ToDelete3));
                Tnew=stats;
                Tnew(ToDeleteFinal,:)=[];

                centers = Tnew.Centroid;
                diameters = mean([Tnew.MajorAxisLength Tnew.MinorAxisLength],2);
                radii = diameters/2;

                imshow(BW2)
                h = viscircles(centers,radii);  

                AreaFrames(i)=Tnew.Area;
                MinorLengthFrames(i)=Tnew.MajorAxisLength;
                MajorLengthFrames(i)=Tnew.MinorAxisLength;

                     
            end


            RawData.AreaFrames=AreaFrames;
            RawData.MinorLengthFrames=MinorLengthFrames;
            RawData.MajorLengthFrames=MajorLengthFrames;
            RawData.Time=time;          
            
           
            %Butterworth filter
            Fs=1/DeltaT;%data sampling
            fc=5;%cut-off frequency
            % Filter lowpass 10Hz
            [b,a] = butter(2, fc/Fs,'low');
            RawData.FilteredArea= filtfilt(b,a,RawData.AreaFrames);

            %Calculation of the strain with filtered data
            Max=max(RawData.FilteredArea);
            RawData.Contraction=-(RawData.FilteredArea-Max);
            RawData.FilteredStrain=-(RawData.FilteredArea-Max)/Max;
            
            % Frequency with fft
            Y=fft(RawData.FilteredStrain);
            L=length(RawData.FilteredStrain);
            fftfreq = 1/DeltaT*[0:floor(L/2)]./L; %vecteur fréquences
            P2 = abs(Y/L);
            P1 = P2(1:floor(L/2)+1);
            P1(2:end-1) = 2*P1(2:end-1);
            [fftpeaks indexfftpeaks]=findpeaks(P1);
            [fftmax indexmax]=max(fftpeaks);
            freqfromfft=fftfreq(indexfftpeaks(indexmax));

            RawData.ContractionFreqFromFFT=freqfromfft;
            
            %Contraction characterization
            MH=0.6*max(RawData.Contraction);
            MP=1/freqfromfft/2*RawData.FrameRate;
            %determination of minima and maxima
            Index_maxima=[];
            Index_minima=[];
            Index_min_final=[];
            
            %the characterization is done on the amplitude of filtered Area
            [Contrac_max, Index_maxima]= findpeaks (RawData.Contraction,'MinPeakHeight',MH,'MinPeakDistance',MP);
            [Contrac_min, Index_minima]= findpeaks (-RawData.Contraction,'MinPeakDistance',MP);
            
            RawData.IndexMax=Index_maxima;
            RawData.ContractMaxStrain=RawData.FilteredStrain(Index_maxima);
            RawData.MeanContractMaxStrain=mean(RawData.ContractMaxStrain);
            
            
            % Remove the minima in excess 

       
            % we make sure that there is only one minimum between 2 max
            if Index_maxima(1)<Index_minima(1) %if the first area max is before the first area min
                for j=1:length(Index_maxima)-1
                    difference=zeros(length(Index_minima),1);
                    difference=Index_maxima(j+1)*ones(size(Index_minima))-Index_minima;
                    Mindiff=sort(difference);
                    k=1;
                    for k=1:length(Mindiff)
                        if Mindiff(k)>0
                            break
                        end
                    end
                    Indexmin=Index_maxima(j+1)-Mindiff(k);
                    Index_min_final=[Index_min_final Indexmin];                   

                end
            end

            if Index_maxima(1)>Index_minima(1) %if the first area min is before the first area max
                for j=1:length(Index_maxima)
                    difference=zeros(length(Index_minima),1);
                    difference=Index_maxima(j)*ones(size(Index_minima))-Index_minima;
                    Mindiff=sort(difference);
                    for k=1:length(Mindiff)
                        if Mindiff(k)>0
                            break
                        end
                    end
                    Indexmin=Index_maxima(j)-Mindiff(k);
                    Index_min_final=[Index_min_final Indexmin];                   

                end
            end
            
            RawData.IndexMin=Index_min_final;
            RawData.ContractMinStrain=RawData.FilteredStrain(Index_min_final);
            RawData.MeanContractMinStrain=mean(RawData.ContractMaxStrain);
            

%             figure,            
%             plot(time,RawData.AreaFramesMilli)
%             axis tight
%             hold on
%             plot(time,RawData.FilteredArea)
%             
            

            

            %Rise time
            RawData.RiseTimeStrain=risetime(RawData.FilteredStrain,RawData.Time);
            RawData.MeanRiseTimeStrain=mean(RawData.RiseTimeStrain);  

            %Fall time
            RawData.FallTimeStrain=falltime(RawData.FilteredStrain,RawData.Time);
            RawData.MeanFallTimeStrain=mean(RawData.FallTimeStrain);
            
            

            %Show data
            figure, plot(RawData.Time,RawData.FilteredStrain), hold on
            plot(RawData.Time(Index_maxima),RawData.ContractMaxStrain,'*'),hold on
            plot(RawData.Time(RawData.IndexMin),RawData.ContractMinStrain,'*')
            
            

            %% Calculation of parameters of the signal


        rel90={};
        Index90={};


        length_index=min(length(RawData.IndexMax),length(RawData.IndexMin));
        rel90{1,i}=zeros(length_index,1);
        Index90{1,i}=zeros(length_index,1);    




        epsilon90=(max(RawData.FilteredStrain)-min(RawData.FilteredStrain))/25;%margin(tolerance) for determination of relaxation 90%
        length_index=min(length(RawData.IndexMax),length(RawData.IndexMin));

        RawData.Amplitude=zeros(1,length_index-1); %initialisation of amplitude
        RawData.ContracTime=zeros(1,length_index-1);%initialisation of time of contraction
        RawData.ContracSpeed=zeros(1,length_index-1);%initialisation of speed of contraction
        RawData.RelaxTime=zeros(1,length_index-1);%initialisation of time of relaxation
        RawData.RelaxSpeed=zeros(1,length_index-1);%initialisation of speed of relaxation
        RawData.IndexRel90=zeros(1,length_index-1);%initialisation of index of rel90
        RawData.RestTime=zeros(1,length_index-1);

        if RawData.IndexMax(1)<RawData.IndexMin(1) %if the first contraction max is before the first contraction min

            %About contraction

            for j=1:length_index-1%we start at the second max to have the amplitude of contraction
                RawData.AmplitudeStrain(j)=RawData.ContractMaxStrain(j+1)-RawData.ContractMinStrain(j);%difference between a max and the previous min of area

                Index_beg=RawData.IndexMin(j);%index of the beginning of strain
                time_beg=RawData.Time(Index_beg);
                Index_end=RawData.IndexMax(j+1);%index of the end of strain
                time_end=RawData.Time(Index_end);
                RawData.ContracTime(j)=time_end-time_beg;%time of strain

                RawData.ContractSpeedStrain(j)=RawData.AmplitudeStrain(j)/RawData.ContracTime(j); %speed of contraction=amplitude/time
            end


            %About relaxation
            for j=1:length_index
                rel90{1,i}(j)=(RawData.ContractMaxStrain(j)-RawData.ContractMinStrain(j))*0.1;%amplitude of relaxation at 90%
                for k=RawData.IndexMax(j):RawData.IndexMin(j)-1%we look for the index of rel90(i,j) between the min and max of contraction
                    if RawData.FilteredStrain(k)<rel90{1,i}(j)+RawData.ContractMinStrain(j)+epsilon90 && RawData.FilteredStrain(k)>rel90{1,i}(j)+RawData.ContractMinStrain(j)-epsilon90%we look for the index of the value of rel90(i,j)
                        Index90{1,i}(j)=k;
                    end
                end
                if Index90{1,i}(j)~=0
                    RawData.RelaxTime(j)=RawData.Time(Index90{1,i}(j))-RawData.Time(RawData.IndexMax(j));%time of relaxation=time(rel90(i,j))-time(max_contraction)
                    RawData.RelaxSpeedStrain(j)=(RawData.ContractMaxStrain(j)-rel90{1,i}(j))/RawData.RelaxTime(j);
                    %speed of relaxation to 90% 
                    RawData.IndexRel90(j)=Index90{1,i}(j);
                    if Index90{1,i}(j)<RawData.IndexMin(j)
                        RawData.RestTime(j)=RawData.Time(RawData.IndexMin(j))-RawData.Time(Index90{1,i}(j));
                    else
                        RawData.RestTime(j)=RawData.Time(RawData.IndexMin(j+1))-RawData.Time(Index90{1,i}(j));
                    end

                end
            end

        else %if the first contraction min is before the first contraction max

            %About contraction
            for j=1:length_index % we start at the first min to have the amplitude of contraction
                RawData.AmplitudeStrain(j)=RawData.ContractMaxStrain(j)-RawData.ContractMinStrain(j);% difference between a max and the next min of area
                Index_beg=RawData.IndexMin(j);%index of the beginning of contraction
                time_beg=RawData.Time(Index_beg);
                Index_end=RawData.IndexMax(j);%index of the end of contraction
                time_end=RawData.Time(Index_end);
                RawData.ContracTime(j)=time_end-time_beg;%time of contraction

                RawData.ContractSpeedStrain(j)=RawData.AmplitudeStrain(j)/RawData.ContracTime(j); %speed of contraction=amplitude/time
            end

            %About relaxation
            for j=1:length_index-1
                rel90{1,i}(j)=(RawData.ContractMaxStrain(j)-RawData.ContractMinStrain(j+1))*0.1;%we take the first min and the second max, amplitude of relaxation at 90%
                for k=RawData.IndexMax(j):RawData.IndexMin(j+1)-1%we look for the index of rel90(i,j) between the min and max of the ECT length
                    if RawData.FilteredStrain(k)<(rel90{1,i}(j)+RawData.ContractMinStrain(j+1)+epsilon90) && (RawData.FilteredStrain(k)>rel90{1,i}(j)+RawData.ContractMinStrain(j+1)-epsilon90)%we look for the index of the value of rel90(i,j)
                        Index90{1,i}(j)=k;
                    end
                end
                if Index90{1,i}(j)~=0
                    RawData.RelaxTime(j)=RawData.Time(Index90{1,i}(j))-RawData.Time(RawData.IndexMax(j));%time of relaxation=time(rel90(i,j))-time(max_contraction)
                    RawData.RelaxSpeedStrain(j)=(RawData.AmplitudeStrain(j)-rel90{1,i}(j))/RawData.RelaxTime(j);
                    RawData.IndexRel90(j)=Index90{1,i}(j);
                    if Index90{1,i}(j)<RawData.IndexMin(j)
                        RawData.RestTime(j)=RawData.Time(RawData.IndexMin(j))-RawData.Time(Index90{1,i}(j));
                    else
                        RawData.RestTime(j)=RawData.Time(RawData.IndexMin(j+1))-RawData.Time(Index90{1,i}(j));
                    end
                end
            end
        end


        %Mean parameters
        ExportVal=struct();
        ExportVal.Name=RawData.Name;
        ExportVal.Duration=RawData.Duration;
        ExportVal.FrameRate=RawData.FrameRate;
        ExportVal.NumFrames=RawData.NumFrames;
        ExportVal.Height=RawData.Height;
        ExportVal.Width=RawData.Width;
        ExportVal.ContractionFreqFromFFT=RawData.ContractionFreqFromFFT;

%        ExportVal.MeanContractMaxStrain=RawData.MeanContractMaxStrain;

        ExportVal.MeanAmpliStrain=mean(RawData.AmplitudeStrain(:));
        ExportVal.MeanRiseTimeStrain=RawData.MeanRiseTimeStrain; 
        ExportVal.MeanFallTimeStrain=RawData.MeanFallTimeStrain;     

        ExportVal.MeanContracTime=mean(RawData.ContracTime(:));
        ExportVal.MeanContractSpeedStrain=mean(RawData.ContractSpeedStrain(:));  
        ExportVal.MeanRelaxTime=mean(RawData.RelaxTime(:));
        ExportVal.MeanRelaxSpeedStrain=mean(RawData.RelaxSpeedStrain(:));
        ExportVal.RestTime=mean(RawData.RestTime(:));
        ExportVal.RestTimeStd=std(RawData.RestTime(:));
        


    %plot points max, min and rel90 as a function of time
    figure,

        plot (RawData.Time(:),RawData.FilteredStrain,'LineWidth',2)
        hold on
        plot (RawData.Time(RawData.IndexMin),RawData.ContractMinStrain,'*','MarkerSize',8,'LineWidth',1.5) %points min
        hold on
        plot (RawData.Time(RawData.IndexMax),RawData.ContractMaxStrain,'*','MarkerSize',8,'LineWidth',1.5) %points max
        hold on
        for j=1:length(rel90{1,i})
            if Index90{1,i}(j)~=0 && rel90{1,i}(j)~=0
                plot (RawData.Time(RawData.IndexRel90(j)),RawData.FilteredStrain(Index90{1,i}(j)),'k*','MarkerSize',8,'LineWidth',1.5) %points rel90
                hold on
            end
        end


    title ('Contraction of the tissue, min, max and relaxation 90%','FontSize', 18)
    xlabel('Time (s)','FontSize', 16)
    ylabel('Ratio of contracted area (A.U.)','FontSize',16)
                %% Another way to calculate the speeds: max and min of the derivative

    %filter the local speed
    figure,  


        Nber_points=length(RawData.FilteredStrain);
        l=RawData.FilteredStrain; %resampled contraction length (delta) of the tissue
        time=RawData.Time;
        l_dot=zeros(Nber_points-1,1);

        for t=1:Nber_points-1
            l_dot(t)=(l(t+1)-l(t))/(time(t+1)-time(t));
        end

        RawData.LocSpeedLength=l_dot; %local speed of the contraction of the tissue

        %Butterworth filter
        delta=mean(diff(RawData.Time));
        Fs=1/delta;%data sampling
    %     Filter lowpass 20Hz
        fc=20;
        [b,a] = butter(2, fc/Fs,'low');
        RawData.FilteredLocSpeed= filtfilt(b,a,RawData.LocSpeedLength);





        MH=max(RawData.FilteredLocSpeed(10:end-1))/3;
        MH2=max(-RawData.FilteredLocSpeed(10:end-1))/3;%relaxation slower
        [pkmax,locmax] = findpeaks(RawData.FilteredLocSpeed,'MinPeakDistance',MP,'MinPeakHeight',MH);
        RawData.MaxContracSpeedStrain=pkmax;
        RawData.IndexMaxContracSpeedStrain=locmax;
        ExportVal.MaxContracSpeedStrain=mean(pkmax);
        [pkmin,locmin] = findpeaks(-RawData.FilteredLocSpeed,'MinPeakDistance',MP,'MinPeakHeight',MH2);
        RawData.MaxRelaxSpeedStrain=pkmin;
        RawData.IndexMaxRelaxSpeedStrain=locmin;
        ExportVal.MaxRelaxSpeedStrainStrain=mean(pkmin);

        plot(RawData.Time(1:end-1),RawData.LocSpeedLength,'k','LineWidth',2)
        hold on,    
        plot(RawData.Time(1:end-1),RawData.FilteredLocSpeed,'k','LineWidth',2)
        hold on,
        plot(RawData.Time(locmax),pkmax,'vr','MarkerSize',8,'LineWidth',1.5)
        hold on,
        plot(RawData.Time(locmin),-pkmin,'vb','MarkerSize',8,'LineWidth',1.5)
        axis tight
        title ('Contraction speed','FontSize', 18)
        xlabel('Time (s)','FontSize', 16)
        ylabel('Contraction speed (A.U./sec)','FontSize', 16)

        %% Save Data


        %Arrhythmia: standard deviation of the time between two maxima
        ExportVal.StdTimeBtwMax=std(diff(RawData.Time(Index_maxima)));

        cd(SavePath);
        SaveName=FinalExcelName;
        TableNorm2=struct2table(ExportVal);
        writetable(TableNorm2,char(SaveName),'WriteMode','append');


        %Saves the data in a matlab variable
        cd(strcat(rootPath,'/motion_analysis')); 
        SaveName2=strcat('Motion Analysis_',fileName,'.mat');
        save(char(SaveName2),'RawData');


end

