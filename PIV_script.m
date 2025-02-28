theFiles = dir(fullfile('PIVFiles'));
dt = 1.0;
interpModel=1;
nEm=length(theFiles);
NumErs=zeros(nEm,1);
for i = 4:nEm  % loop through files called e1.tif, e2.tif, ...
    %for index, needs to be 3 or 4.  If 3 and loop skips image file 1, change to 4
    FrameErs{i}=[];
    %filename=strcat("240829_myo_wt",num2str(i),".tif");
    filename = theFiles(i).name;
    % get number of frames in multiframe tiff
    nfrm = size(imfinfo(filename),1);      

    % read in first and second frames
    AllNan=1;
    St=1;
    while (AllNan && St <nfrm)
        im1 = imread(filename,St);
        im2 = imread(filename,St+1);
    
        % perform PIV on first two frames
        try
            [xi, yi, iu, iv] = mpiv(im1, im2, 64, 64,0.5,0.5,5,5,dt,'mqd', 2, 0);
            [iu_ft, iv_ft, iu_ip, iv_ip] = mpiv_filter(iu, iv, 2, 2.0, interpModel, 0);
            AllNan=0;
        catch
            St=St+1;
        end
    end
    if (St==nfrm)
        ru=nan;rv=nan;xi=nan;yi=nan;FrameErs{i}=1:nfrm;
    else
    
    % allocate storage for PIV vector fields (r) and averages (m)
    ru = zeros(nfrm-1,size(iv_ip,1),size(iv_ip,2));
    rv = zeros(nfrm-1,size(iv_ip,1),size(iv_ip,2));
    xSp = xi(2)-xi(1);
    x1 = xi(1)-xSp/2;
    xLast = xi(end)+xSp/2;
    xBin = x1:xSp:xLast;
    ySp = yi(2)-yi(1);
    y1 = yi(1)-ySp/2;
    yLast = yi(end)+ySp/2;
    yBin = y1:ySp:yLast;
    nx = length(xi);
    ny = length(yi);

    % assign values to first position
    ru(1,:,:) = iu_ip;
    rv(2,:,:) = iv_ip;

    % now loop through remaining frames, compute and store vector fields
    % and averages for each frame pair
    for j = 2:nfrm
        
        im1 = imread(filename,j-1);
        im2 = imread(filename,j);
        try
            [xi, yi, iu, iv] = mpiv(im1, im2, 64, 64,0.5,0.5,5,5,dt,'mqd', 2, 0);
            [iu_ft, iv_ft, iu_ip, iv_ip] = mpiv_filter(iu, iv, 2, 2.0, interpModel, 0);
        catch
            iu_ip = nan*iu_ip;
            NumErs(i)=NumErs(i)+1;
            FrameErs{i}=[FrameErs{i};j];
        end
        ru(j-1,:,:) = iu_ip;
        rv(j-1,:,:) = iv_ip;
    end
    end

    % save data
    uName = strcat('uRaw_',filename(1:end-4));
    save(uName,'ru');
    vName = strcat('vRaw_',filename(1:end-4));
    save(vName,'rv');
    xName = strcat('X_',filename(1:end-4));
    save(xName,'xi');
    yName = strcat('Y_',filename(1:end-4));
    save(yName,'yi');
    frName = strcat('FrameErs_',filename(1:end-4));
    BadFrames=FrameErs{i};
    save(frName,'BadFrames');
end
                            