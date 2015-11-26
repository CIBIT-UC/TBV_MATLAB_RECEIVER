%% Config
clear, clc;
addpath('utils')
addpath('functions')

%% Create Connection

configs.TBV_IP = '192.168.1.200';
configs.TBV_PORT = 55555;

tbvNetInt = TBVNetworkInterface( TBVclient( configs.TBV_IP, configs.TBV_PORT ) );

tbvNetInt.createConnection();

%% ROI Coordinates

n_rois = tbvNetInt.tGetNrOfROIs();

coordsOfVoxelsOfROI = cell([1 n_rois]);

for i=1:n_rois
    coordsOfVoxelsOfROI{i} = tbvNetInt.tGetAllCoordsOfVoxelsOfROI(i - 1);
end

%% Configs
[xDim, yDim, zDim] = tbvNetInt.tGetDimsOfFunctionalData();

%% Time iteration

time = 1;
counter = 0;

currentTime = tbvNetInt.tGetCurrentTimePoint;
expectedTime = tbvNetInt.tGetExpectedNrOfTimePoints;

ROImeansM1 = zeros(expectedTime,n_rois);
ROImeansM2 = zeros(expectedTime,n_rois);
ROImeansM3 = zeros(expectedTime,n_rois);
ROImeansM4 = zeros(expectedTime,n_rois);

timeM1 = zeros(expectedTime,1);
timeM2 = zeros(expectedTime,1);
timeM3 = zeros(expectedTime,1);
timeM4 = zeros(expectedTime,1);

while time <= currentTime
    
    if time == currentTime && currentTime ~= expectedTime
        disp('Waiting...')
        pause(5)
        counter = counter + 1;
        if counter == 5
            break;
        end
        
    else
        %---Method 1 --> Get Mean ROI
        [ROImeansM1(time,:),timeM1(time)] = method1(n_rois,tbvNetInt,time-1);
        
        %---Method 2 --> All Voxels
        [ROImeansM2(time,:),timeM2(time)] = method2(n_rois,tbvNetInt,coordsOfVoxelsOfROI,xDim,yDim,time-1);
        
        %---Method 3 --> All Voxels Raw
        [ROImeansM3(time,:),timeM3(time)] = method3(n_rois,tbvNetInt,coordsOfVoxelsOfROI,xDim,yDim,time-1);
        
        %---Method 4 --> Voxels 1by1
        [ROImeansM4(time,:),timeM4(time)] = method4(n_rois,tbvNetInt,coordsOfVoxelsOfROI,time-1);
        
        fprintf('Time %d \n',time);
        
        time = time + 1;
        counter = 0;
    end
    
    currentTime = tbvNetInt.tGetCurrentTimePoint;
    
end

%% Plots

timeaxis = 1:expectedTime;
figure()
    plot(timeaxis,timeM1,'--o',timeaxis,timeM2,'--o', ... 
        timeaxis,timeM3,'--o',timeaxis,timeM4,'--o')
    grid on; xlim([1 expectedTime]);
    legend('Method 1','Method 2','Method 3','Method 4');

