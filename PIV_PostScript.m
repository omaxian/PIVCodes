nBins = 20;
WindowSize = 15; % Has to be ODD
Name = "NMYNikonEd";
EmbryoWindows;
AvgPl=1;
IndsPl=0;

%for iT=6:-1:2
% define last frames for late maintenance
t2 = EndMaint - 60;%30*(iT-1);
t1 = EndMaint - 120;%30*iT; % Last minute of maintenance phase
nFrames=nFrames-1;

% compute x spacing on per-embryo basis
w = x2-x1;
d = w/nBins;
nTrial=length(x1);
n1 = zeros(nTrial,nBins+1);
OneBrt1 = zeros(nTrial,nBins+1);
BDPos = nan*zeros(nTrial,201);
MaxFlowSpeed = nan*zeros(nTrial,201);
ZeroWindowIndex = 101;
MeanIntensities = nan*zeros(nTrial,2001);
ZeroIndex=1001;

% compute x spacing on per-embryo basis
for i = 1:nTrial
    filename = strcat(Name,'_r',num2str(i),'.mat');
    filenameX = strcat(Name,'_x',num2str(i),'.mat');
    load(filename);
    load(filenameX);
    dxHave = xi(2)-xi(1);
    xHaveStart=xi(1);

    % average across y
    try
    mid = ceil(size(ru,3)/2);
    w = floor(mid/2);
    m = mean(ru(:,:,(mid-w):(mid+w)),3,'omitnan');
    catch
    mid = ceil(size(r,3)/2);
    w = floor(mid/2);
    m = mean(r(:,:,(mid-w):(mid+w)),3,'omitnan');
    end
    % average over time
    if (t1(i) > 0)
        m1 = mean(m(t1(i):t2(i),:),'omitnan');
    else
        if (t2(i) < 0)
            m1 = nan*m(1,:);
        else
            m1 = mean(m(1:t2(i),:),'omitnan');
        end
    end
    
    x = x1(i):d(i):x2(i);
    f = 1 + floor((x-xHaveStart)/dxHave); 
    r = rem((x-xHaveStart),dxHave)/dxHave;
    n1(i,:) = (1-r).*m1(f) + r.*m1(1+f);
    
    % Plot boundary position and max flows over time
%     AllFl = (1-r).*m(:,f) + r.*m(:,1+f);
%     AllFl(isnan(AllFl))=0;
%     % Average over 15 s window, centered at 0
%     MiddleInds = [fliplr(EndMaint(i)-WindowSize:-WindowSize:1) ...
%         EndMaint(i):WindowSize:nFrames(i)];
%     LeftInt = floor(WindowSize/2);
%     StartAvgInds = max(MiddleInds-LeftInt,1);
%     EndAvgInds = min(MiddleInds+LeftInt,nFrames(i));
%     nWind = length(StartAvgInds);
%     FlEveryWind = zeros(nWind,nBins+1);
%     for iW=1:nWind
%         FlEveryWind(iW,:)=mean(AllFl(StartAvgInds(iW):EndAvgInds(iW),:),'omitnan');
%     end
%     % Find window corresponding to zero
%     ZeroWindow = find(MiddleInds==EndMaint(i));
%     offset = ZeroWindowIndex-ZeroWindow;
%     MaxFlowSpeed(i,offset+1:offset+nWind)=max(abs(FlEveryWind'));
end
   
% Plots
% Individual flows and mean
% figure(1)
Fac = MicronsPerPixel*60;
n1 = n1*Fac;
if (IndsPl)
plot((0:nBins)/nBins,n1,'Color',LightColor,'LineWidth',1.0)
hold on
end
if (AvgPl)
errorbar((0:nBins)/nBins,mean(n1,'omitnan'),std(n1,'omitnan')/sqrt(nTrial),'-',...
    'LineWidth',2.0,'Color',DarkColor);
hold on
end
% Max flow over time
% figure(2)
% tvals = -WindowSize*100:WindowSize:WindowSize*100;
% nnz = nTrial-sum(isnan(MaxFlowSpeed));
% MaxFlowSpeed=MaxFlowSpeed*Fac;
% if (IndsPl)
% plot(MaxFlowSpeed,tvals,'Color',LightColor,'LineWidth',1.0)
% hold on
% end
% if (AvgPl)
% errorbar(mean(MaxFlowSpeed,'omitnan'),tvals,[],[],...
%     std(MaxFlowSpeed,'omitnan')./sqrt(nnz),std(MaxFlowSpeed,'omitnan')./sqrt(nnz),...
%     'LineWidth',2.0,'Color',DarkColor)
% hold on
% end
% set(gca,'YDir','Reverse')
% ylim([-400 0])
% xlabel('Max flow speed $\mu$m/min')
% ylabel('Time (s)')
%end