%% Config
clear, clc;
addpath('utils')

%% Create Connection

configs.TBV_IP = '192.168.1.200';
configs.TBV_PORT = 55555;

tbvNetInt = TBVNetworkInterface( TBVclient( configs.TBV_IP, configs.TBV_PORT ) );

tbvNetInt.createConnection();

%% ROI Coordinates

n_rois = tbvNetInt.tGetNrOfROIs();

coordsOfVoxelsOfROI = cell([1 n_rois]);

for i=1:n_rois
    coordsOfVoxelsOfROI{i} = tbvNetInt.tGetAllCoordsOfVoxelsOfROI(i);
end

%% Configs
timePoint = 100; %To iterate after
[xDim, yDim, zDim] = tbvNetInt.tGetDimsOfFunctionalData();

%% Method 1 - TBV
tic

ROImeansM1 = zeros(1,n_rois);

for i=1:n_rois
    ROImeansM1(i) = tbvNetInt.tGetMeanOfROIAtTimePoint( i, timePoint );
end

timeM1 = toc;

%% Method 2 - All Voxels
tic

voxelsValsM2 = tbvNetInt.tGetValueOfAllVoxelsAtTime( timePoint );
% voxelsMatrix = reshape(voxelsVals, [zDim yDim xDim]);

ROImeansM2 = zeros(1,n_rois);

for i=1:n_rois
    
   ROImeansM2(i) = mean( voxelsValsM2(coordsOfVoxelsOfROI{i}(:,3)*xDim*yDim ...
       + coordsOfVoxelsOfROI{i}(:,2)*xDim + coordsOfVoxelsOfROI{i}(:,1) + 1) );
    
end

timeM2 = toc;

%% Method 3 - All Raw Voxels
tic

voxelsValsM3 = tbvNetInt.tGetRawValueOfAllVoxelsAtTime( timePoint );

ROImeansM3 = zeros(1,n_rois);

for i=1:n_rois
    
   ROImeansM3(i) = mean( voxelsValsM3(coordsOfVoxelsOfROI{i}(:,3)*xDim*yDim ...
       + coordsOfVoxelsOfROI{i}(:,2)*xDim + coordsOfVoxelsOfROI{i}(:,1) + 1) );
    
end

timeM3 = toc;

%% Method 4 - Voxels 1by1
tic

ROImeansM4 = zeros(1,n_rois);

for i=1:n_rois
    
    temp = zeros(1,size(coordsOfVoxelsOfROI{i},1));
    for j=1:size(coordsOfVoxelsOfROI{i},1)
        
        temp(j) = tbvNetInt.tGetValueOfVoxelAtTime(coordsOfVoxelsOfROI{i}(j,1) ...
            ,coordsOfVoxelsOfROI{i}(j,2),coordsOfVoxelsOfROI{i}(j,3),timePoint);
    end
    
    ROImeansM4(i) = mean( temp );
    
end

timeM4 = toc;

%% Time iteration

time = 1;
counter = 1;
currentTime = tbvNetInt.tGetCurrentTimePoint;
expectedTime = tbvNetInt.tGetExpectedNrOfTimePoints;

ROImeansM4 = zeros(expectedTime,n_rois);

while time <= currentTime
    
    if time == currentTime
        pause(0.5)
        counter = counter + 1;
        if counter == 5
            break;
        end
        
    else
        time = time + 1;
        
        %Methods

        for i=1:n_rois
            
            temp = zeros(1,size(coordsOfVoxelsOfROI{i},1));
            for j=1:size(coordsOfVoxelsOfROI{i},1)
                
                temp(j) = tbvNetInt.tGetValueOfVoxelAtTime(coordsOfVoxelsOfROI{i}(j,1) ...
                    ,coordsOfVoxelsOfROI{i}(j,2),coordsOfVoxelsOfROI{i}(j,3),time);
            end
            
            ROImeansM4(time,i) = mean( temp );
            
        end
        
        
    end
    
    
    
end



% 
% 
% %% Correlations from TBV
% 
% n_rois = tbvNetInt.tGetNrOfROIs();
% n_con = (n_rois*(n_rois-1))/2; % #Connections
% 
% windowSize = 30;
% timePoint = 100;
% 
% [wSize, pearsonC] = tbvNetInt.tGetPearsonCorrelation(windowSize);
% fprintf('---> Pearson Correlation with window size %d = %f %f %f \n',wSize,pearsonC)
% 
% [wSize , tPoint , pearsonC_t] = tbvNetInt.tGetPearsonCorrelationAtTimePoint(windowSize,timePoint);
% fprintf('---> Pearson Correlation at time point %d, window size %d = %f %f %f \n',wSize,tPoint,pearsonC_t)
% 
% [wSize, partialC] = tbvNetInt.tGetPartialCorrelation(windowSize);
% fprintf('---> Partial Correlation with window size %d = %f %f %f \n',wSize,partialC)
% 
% [wSize , tPoint , partialC_t] = tbvNetInt.tGetPartialCorrelationAtTimePoint(windowSize,timePoint);
% fprintf('---> Partial Correlation at time point %d, window size %d = %f %f %f \n',wSize,tPoint,partialC_t)
% 
% %% Fig 1
% figure()
%     plot(pearsonC_t,'o'); hold on
%     plot(partialC_t,'o'); hold on
%     xlim([0 4]); ylim([-1 1]); grid on;
%     legend('Pearson','Partial','Location','Southwest')
%     
% %% Correlations from ROI Means
% clear roi1 roi2 ROIcorr
% 
% ROIcorr = zeros(n_con,1);
% 
% C = combnk(1:n_con,2);
% 
% for c=1:n_con
%     roi1 = zeros(windowSize,1);
%     roi2 = zeros(windowSize,1);
%     idx = 1;
%     for i=timePoint-windowSize-1:timePoint
%         roi1(idx) = tbvNetInt.tGetMeanOfROIAtTimePoint(C(c,1),i);
%         roi2(idx) = tbvNetInt.tGetMeanOfROIAtTimePoint(C(c,2),i);
%         idx = idx+1;
%     end
%     ROIcorr(c) = corr(roi1,roi2);
%     
% end
% 
% %% Real time
% i = 1;
% currentTime = tbvNetInt.tGetCurrentTimePoint;
% 
% while i < currentTime
%     
%     if i ~= currentTime
%         [~ , ~ , RTpearson(i,:)] = tbvNetInt.tGetPearsonCorrelationAtTimePoint(windowSize,i);
%         i = i+1;
%     else
%         pause(0.5)
%     end
% 
% end
% 
% %% Fig 2
% figure()
%     plot(RTpearson(:,1),'--o'); hold on
%     plot(RTpearson(:,2),'--o'); hold on
%     plot(RTpearson(:,3),'--o'); hold on
%     grid on;
%     legend('ROI 12','ROI 13','ROI 23','Location','Southwest')
%     