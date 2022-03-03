function [out]= EgoCentricRateMapMultiFrame(root, varargin)
% out = EgoCentricRateMap_fastDev(root, varargin)
% 
% Creates egocentric ratemaps for any number of reference frames. 
% 
% arguments:
%   videoSamp: Sample every nth video frame (1)
%   degSamp: For edges, how many degrees in each "spoke" (3)
%   heading: 0 to use HD, 1 to use heading (1)
%   distanceBins: Spatial binning in cm (0:2.5:50)
%   smooth: 2D Gaussian smoothing [std_x std_y, width] ([5 5 5])
%   boundaryMode: How to determine the structure of the environment:
%                 0: Use the extremities of x and y trajectory (open field
%                 only)
%                 1: Interactive assignment (see below)
%                 Structure: Use a predefined structure (eg from out.envStruct)
% 
% 
% For interactive definition of the environment, follow these steps:
%   1) Trajectory plot will pop up. Controls:
%      Left click: to set a control point (eg: corner of an
%      environment/hole/object/point)
%      Space Bar: To end the current object and continue to the next one
%      Escape key: To finish assigning points
% 
%   2) Command window will now ask for a radius to assign to any "point"
%   objects (eg: those with only one control point). 
% 
%   3) Command window will now ask for an assignment of types to each
%   structure. The options are:
% 
%       A) Environment: The outer boundary of the enclosure
% 
%       B) Exclusion: Some boundary inside the enclosure which the animal
%       can not enter, and should be treated as a "wall" just as the outer
%       boundaries are (eg: a cut out, or walled off internal section of
%       environment)
% 
%       C) Object: Includes any internal objects in the environment which
%       should each generate their own reference frame. This includes 
%       multi-walled objects and point objects. (The primary difference
%       being how the closest point on the object to the animal is
%       calculated).
% 
%   Note that an exclusion is treated as part of the environment -- they
%   are used together to create a single ratemap. An additional ratemap is
%   generated for each object.
   
    %%
    warning('off','stats:glmfit:IterationLimit');
    warning('off','curvefit:fit:noStartPoint');
    warning('off','curvefit:fittype:sethandles:WeibullxMustBePositive');
    warning('off','stats:glmfit:BadScaling');
    warning('off','MATLAB:nargchk:deprecated');
    
    %% Setup & Parse
    p = inputParser;
    p.addParameter('videoSamp', 1);                    % calculate every X frames of video
    p.addParameter('degSamp', 3);                      % Degree bins
    p.addParameter('heading', 1);                       % use heading (0--> Head direction)
    p.addParameter('distanceBins', 0:2.5:62.5);           % How far to look (cm)
    p.addParameter('boundaryMode', 1);                  % 0-> autolines, 1-> click, mat->useit
    p.addParameter('smooth', [5 5 5])
    p.addParameter('envStruct', []);                    % to pass in instead of GUI
    p.addParameter('spikeMatch', []); 
    p.parse(varargin{:});
   
    fn = fieldnames(p.Results);
    for i = 1:length(fn)
        eval([fn{i} '= (p.Results.' fn{i} ');']);
    end  
    
    %% Get behavior information
    my = max(root.b_y);
    %root.b_y = my-(root.b_y); #TODO:Why was this causing multiple repeats
    %of dis, and offsetting theta?
    
    rx = CMBHOME.Utils.ContinuizeEpochs(root.x);
    ry =  CMBHOME.Utils.ContinuizeEpochs(root.y);
    brx = root.b_x; bry = root.b_y;
    
    maxDist = max([pdist2([max(rx),min(ry)],[min(rx),max(ry)]),pdist2([min(rx),min(ry)],[max(rx),max(ry)])])*root.spatial_scale;
    
    if heading==0
        hd = wrapTo2Pi(deg2rad(CMBHOME.Utils.ContinuizeEpochs(root.headdir)));
        bhd = wrapTo2Pi(deg2rad(root.b_headdir));
    else
        hd = [0; atan2(diff(ry),diff(rx))];   % use this for heading
        hd = wrapTo2Pi(hd);
        bhd = wrapTo2Pi([0; atan2(diff(root.b_y), diff(root.b_x))]);
    end
    
    %% Get structure of environment
    if isempty(envStruct)
    
            if numel(boundaryMode)==1 && ~isstruct(boundaryMode)
                if boundaryMode == 0 
                    % Headless / auto. Automatically detect the edges. Only works for a //
                    % rectangular box with no insertions.

                    p = [-1000 -1000];
                    d = (brx-p(1)).^2 + (bry-p(2)).^2;
                    [~,ind] = min(d);
                    ll = [brx(ind) bry(ind)];

                    p = [1000 -1000];
                    d = (brx-p(1)).^2 + (bry-p(2)).^2;
                    [~,ind] = min(d);
                    lr = [brx(ind) bry(ind)];

                    p = [1000 1000];
                    d = (brx-p(1)).^2 + (bry-p(2)).^2;
                    [~,ind] = min(d);
                    ur = [brx(ind) bry(ind)];

                    p = [-1000 1000];
                    d = (brx-p(1)).^2 + (bry-p(2)).^2;
                    [~,ind] = min(d);
                    ul = [brx(ind) bry(ind)];


                    QP = [ll;lr;ur;ul];
                    edg = splitter(QP);
                    envStruct.edg = edg;
                    envStruct.rad = NaN;
                    envStruct.type = 1;

                elseif boundaryMode == 1 && ~isstruct(boundaryMode)
                    envStruct = findEdges(root);
                end
            else
                envStruct = boundaryMode;
            end
    else
        
    end
    
    %% Calculate distances
    [dis, ang] = subfunc(brx,bry,bhd, envStruct, degSamp);
    dis{1} = dis{1}*root.spatial_scale;
    dis{2} = cellfun(@(x) x*root.spatial_scale, dis{2},'UniformOutput',0);
    dis{3} = cellfun(@(x) x*root.spatial_scale, dis{3},'UniformOutput',0);   
    [dis{4},angMin] = min(dis{1},[],2); % distance and angle to closest wall
    
    %% Calculate raw maps:
    thetaBins = -pi:deg2rad(degSamp):pi;
    
    if numel(dis{1}) > 0
        distanceBins = 0:2.5:62.5;
        [occ{1}, nspk{1}, rm_ns{1}] = CalcMaps(dis{1}, distanceBins, thetaBins, root, spikeMatch);
    else
        occ{1} = []; nspk{1} = []; rm_ns{1} = [];
    end
    
    for i = 1:length(dis{2})
        distanceBins = 0:2.5:maxDist;
        [occ{end+1}, nspk{end+1}, rm_ns{end+1}] = CalcMaps_Point(dis{2}{i}, ang{2}{i}, distanceBins, thetaBins, root, spikeMatch);
    end
    
    for i = 1:length(dis{3})
        distanceBins = 0:2.5:maxDist;
        [occ{end+1}, nspk{end+1}, rm_ns{end+1}] = CalcMaps(dis{3}{i}, distanceBins, thetaBins, root, spikeMatch);
    end
    
    if numel(dis{4}) > 0
%         distanceBins = 0:2.5:62.5;
        [occ{end+1}, nspk{end+1}, rm_ns{end+1}] = CalcMaps_Point(dis{4}, thetaBins(angMin)', distanceBins, thetaBins, root, spikeMatch);
    else
        occ{end+1} = []; nspk{end+1}=[]; rm_ns{end+1}=[];
    end    
    
    %% Smoothing
    occ_sm = cellfun(@(x) smoothing(x, smooth(1:2), smooth(3)), occ,'UniformOutput',0);
    nspk_sm = cellfun(@(x) smoothing(x, smooth(1:2), smooth(3)), nspk,'UniformOutput',0);
    rm = cellfun(@(x,y) x./y*root.fs_video, nspk, occ,'UniformOutput',0);
    rm_sm = cellfun(@(x) smoothing(x, smooth(1:2), smooth(3)), rm,'UniformOutput',0);
    
    %% EBC test stuff
    
    fs = root.fs_video;
    [MRL, PO, sDistMin, sDistMax]= cellfun(@(x,y) ebcTests(x, y, fs), nspk_sm, occ_sm, 'UniformOutput',0);
    
    %% Package the output
    out.distanceBins = distanceBins;
    out.envStruct = envStruct;

    out.EgoCentricRateMapParams = varargin;
    out.rm = rm;
    out.occ = occ;
    out.nspk = nspk;
    out.occ_sm = occ_sm;
    out.nspk_sm = nspk_sm;
    out.rm_sm = rm_sm;
    
    out.thetaBins = thetaBins;
    out.distanceBins = distanceBins;
    
    out.ang = ang;
    out.dis = dis;
    
    out.MRL = MRL;
    out.PO = PO;
    out.sMin = sDistMin;
    out.sMax = sDistMax;
    
    
end

function envStruct = findEdges(root)
    ifEscape = 0;
    h=figure();

    while ~ifEscape  
        disp('- Left Click to set boundary of wall')
        disp('- Press space to begin new object')
        disp('- For circles center on loc, press space after one click')
        disp('- Escape to finish')
        disp('** Do NOT "close" environments')
        figure(h); 
        clf
        
        [occupancy, xdim, ydim]=root.Occupancy([],[],1,2);
        set(gca,'YDir','Normal'); %colormap(jet);
        clim=get(gca,'clim');set(gca,'clim',clim/50);
        hold on
        plot(CMBHOME.Utils.ContinuizeEpochs(root.b_x),CMBHOME.Utils.ContinuizeEpochs(root.b_y),'k');
        QP = [];
        
        set(h,'Name','Select Corners of Walls. \n Esc--> done. \n **Do not complete!**')

        button = 1;

        while button~=27
            [x,y,button] = ginput(1);

            clf
            
            %imagesc(xdim,ydim,occupancy'); 
            set(gca,'YDir','Normal'); %colormap(jet);
            clim=get(gca,'clim');set(gca,'clim',clim/50);
            hold on
            plot(CMBHOME.Utils.ContinuizeEpochs(root.b_x),CMBHOME.Utils.ContinuizeEpochs(root.b_y),'k');
            
            if ~isempty(QP)
                plot(QP(:,1),QP(:,2),'r')
                plot(QP(:,1),QP(:,2),'ro','MarkerFaceColor','r')
            end

            if button == 32 %space bar
                QP = [QP; NaN NaN];
            elseif button~=27
                QP = [QP; x y];
            end

            plot(QP(:,1),QP(:,2),'r')
            plot(QP(:,1),QP(:,2),'ro','MarkerFaceColor','r')

        end

        
        %% Ask for verification
        edg = splitter(QP);
        points = (cellfun(@(x) size(x,1)==1,edg));
        for i = 1:length(points)
            if points(i)
                rad(i) = input(['Enter radius for point ' num2str(i) ': ']);
            else
                rad(i) = NaN;
            end
        end
        
        %%
        colors = {'r','g','b',[.5 .5 .5], 'y','m','c'};
        clf;
        set(h,'Name','Verify. 0--> Try again; 1--> Confirm')
        hold on
        
        for m = 1:numel(edg)
            cm = mod(m, length(colors)) + 1;   
            if size(edg{m},1) > 1
                for n = 1:size(edg{m},1)
                    sp = squeeze(edg{m}(n,:,1));
                    ep = squeeze(edg{m}(n,:,2));
                    plot([sp(1) ep(1)],[sp(2) ep(2)],'o','MarkerFaceColor',colors{cm},'MarkerEdgeColor',colors{cm})
                    plot([sp(1) ep(1)],[sp(2) ep(2)],'MarkerFaceColor',colors{cm})
                end
                h2 = fill(edg{m}(:,1), edg{m}(:,2),colors{cm},'EdgeColor',colors{cm});
                set(h2,'facealpha',0.5);
                
            else
                X = [edg{m}(1)-rad(m) edg{m}(2)-rad(m) 2*rad(m) 2*rad(m)];
                h2 = rectangle('Position',X,'Curvature',[1 1],'FaceColor',colors{cm},'EdgeColor',colors{cm});
                %set(h,'facealpha',0.5);
            end
        end
        plot(CMBHOME.Utils.ContinuizeEpochs(root.b_x),CMBHOME.Utils.ContinuizeEpochs(root.b_y),'k');

        % set or repeat
        while button ~=48 && button~=49
            [~,~,button]=ginput(1);
        end
        ifEscape = button==49;

    end

    %% 
    for m = 1:numel(edg)
        type(m) = input(['Region ' num2str(m) ' Environment(1), Exclusion (2), or Object (3)? ']);
    end
    
    %% 
    for m = 1:numel(edg)
        envStruct(m).edg = edg{m};
        envStruct(m).rad = rad(m);
        envStruct(m).type = type(m);
    end
    
    close(h);
    drawnow();
    
end

%%
function edg = splitter(QP)
    
    inds = find(isnan(QP(:,1)));
    xs=CMBHOME.Utils.SplitVec(QP(:,1), @(x) isnan(x));
    ys=CMBHOME.Utils.SplitVec(QP(:,2), @(x) isnan(x));
    
    % split corners
    for m = 1:size(xs,1)
        QP2{m} = [xs{m} ys{m}];
        QP2{m}(find(isnan(QP2{m}(:,1))),:) = [];
    end
    
    for m = 1:numel(QP2)
        for n = 1:size(QP2{m},1)
            sp = n;ep=n+1;
            if ep>size(QP2{m},1), ep=1;end
            edg{m}(n,:,1) = [QP2{m}(sp,1) QP2{m}(sp,2)];
            edg{m}(n,:,2) = [QP2{m}(ep,1) QP2{m}(ep,2)];
        end
    end

end

%%
function [dis, ang] = subfunc(rx,ry,hd, envStruct, degSamp)
    
    enc = arrayfun(@(x) x.type == 1 | x.type==2, envStruct);
    
    if any(enc > 0)
        dists = calcDisAng_Env(rx, ry, hd, envStruct(enc), degSamp);
    else
        dists = [];
    end
    
    points = find(arrayfun(@(x) x.type==3 & ~isnan(x.rad), envStruct));
    for i = 1:length(points)
        [dist2{i}, ang2{i}] = calcDisAng_Points(rx,ry,hd, envStruct(points(i)));
    end
    if numel(points) == 0, dist2={}; ang2={};end
    
    
    objs = find(arrayfun(@(x) x.type==3 & isnan(x.rad), envStruct));
    for i = 1:length(objs)
        [dist3{i}] = calcDisAng_Obj(rx,ry,hd, envStruct(objs(i)), degSamp);
    end
    if numel(objs) == 0, dist3={}; end
    
    dis = {dists dist2 dist3};
    ang = {[], ang2, []};
    
end

%%
function [dis] = calcDisAng_Env(rx,ry,hd, envStruct, degSamp)
    % Calculate distance to a boundary, using search spokes. Necessary when
    % any part on an edge could be responsible for a cell's firing.

    if numel(envStruct)==1
        edg = envStruct.edg;
    else
        % concat to edgepoint (TODO: for more than 1 exclusion region)
        edg = cat(1,envStruct(1).edg, envStruct(2).edg);
    end
    
    mxd = sqrt((max(rx)-min(rx))^2 + (max(ry)-min(ry))^2);
    degs = deg2rad(-180:degSamp:180);
        
    dis = NaN(numel(rx),size(edg,1), numel(degs));
    dir = dis;
    
    for i = 1:size(edg,1)
        x1=edg(i,1,1);x2=edg(i,1,2);
        y1=edg(i,2,1);y2=edg(i,2,2);
        
        for h = 1:numel(degs)
            hdof=degs(h);
            y3=ry;x3=rx;
            y4=ry+mxd*sin(hd+hdof);
            x4=rx+mxd*cos(hd+hdof);
            
            %https://urldefense.proofpoint.com/v2/url?u=https-3A__en.wikipedia.org_wiki_Line-25E2-2580-2593line-5Fintersection-23Intersection-5Fof-5Ftwo-5Flines&d=DwIGAg&c=-35OiAkTchMrZOngvJPOeA&r=X8VhOQfcglc0O635i6NluA&m=iclUM1400u7J2kCALU0EpG45mkV-2p-N32VFkpC77ro_pIFR5Yc56Ss8Lw-CTq7m&s=IRFqbN3oi9QPMu7r8T-h9584h5gg39Xay9UYErVSMGY&e= 
            px1 = (x1.*y2-y1.*x2).*(x3-x4) - (x1-x2).*(x3.*y4-y3.*x4);
            px2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
            px  = px1./px2;
            
            py1 = (x1.*y2-y1.*x2).*(y3-y4) - (y1-y2).*(x3.*y4-y3.*x4);
            py2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
            py = py1./py2;

            d = sqrt((ry-py).^2 + (rx-px).^2);
            dis(:,i,h) = d;
            
            % need to filter down to the right direction ...
            dir(:,i,h) = wrapToPi(atan2(py-ry,px-rx)-(hd+hdof));
            
            % oh ... we were allowing forever.... filter by bounding box
            bb = [min([x1 x2]) max([x1 x2]); min([y1 y2]) max([y1 y2])];  
            % |xmin, xmax|
            % |ymin, ymax|
            indexes = ~(px>=bb(1,1) & px<=bb(1,2) & py>=bb(2,1) & py<=bb(2,2));
            dis(indexes,i,h) = NaN;
        end
        
    end
    
    
    dis(dis>mxd) = NaN;
    dis(abs(dir)>pi/4) = NaN;
    
    %% output
    dis=squeeze(nanmin(dis,[],2));
    dd=repmat(degs,size(rx,1),1) + repmat(hd,1,numel(degs));
    dx=dis.*cos(dd); dy=dis.*sin(dd);
    ey=dy+repmat(ry,1,numel(degs));
    ex=dx+repmat(rx,1,numel(degs));
    
end

%%
function [dis] = calcDisAng_Obj(rx,ry,hd, envStruct, degSamp)

   % Calculate distance to a boundary, using search spokes. Necessary when
    % any part on an edge could be responsible for a cell's firing.

    if numel(envStruct)==1
        edg = envStruct.edg;
    else
        % concat to edgepoint (TODO: for more than 1 exclusion region)
        edg = cat(1,envStruct(1).edg, envStruct(2).edg);
    end
    
    mxd = sqrt((max(rx)-min(rx))^2 + (max(ry)-min(ry))^2);
    degs = deg2rad(-180:degSamp:180);
        
    dis = NaN(numel(rx),size(edg,1), numel(degs));
    dir = dis;
    
    for i = 1:size(edg,1)
        x1=edg(i,1,1);x2=edg(i,1,2);
        y1=edg(i,2,1);y2=edg(i,2,2);
        
        for h = 1:numel(degs)
            hdof=degs(h);
            y3=ry;x3=rx;
            y4=ry+mxd*sin(hd+hdof);
            x4=rx+mxd*cos(hd+hdof);
            
            %https://urldefense.proofpoint.com/v2/url?u=https-3A__en.wikipedia.org_wiki_Line-25E2-2580-2593line-5Fintersection-23Intersection-5Fof-5Ftwo-5Flines&d=DwIGAg&c=-35OiAkTchMrZOngvJPOeA&r=X8VhOQfcglc0O635i6NluA&m=iclUM1400u7J2kCALU0EpG45mkV-2p-N32VFkpC77ro_pIFR5Yc56Ss8Lw-CTq7m&s=IRFqbN3oi9QPMu7r8T-h9584h5gg39Xay9UYErVSMGY&e= 
            px1 = (x1.*y2-y1.*x2).*(x3-x4) - (x1-x2).*(x3.*y4-y3.*x4);
            px2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
            px  = px1./px2;
            
            py1 = (x1.*y2-y1.*x2).*(y3-y4) - (y1-y2).*(x3.*y4-y3.*x4);
            py2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
            py = py1./py2;

            d = sqrt((ry-py).^2 + (rx-px).^2);
            dis(:,i,h) = d;
            
            % need to filter down to the right direction ...
            dir(:,i,h) = wrapToPi(atan2(py-ry,px-rx)-(hd+hdof));
            
            % oh ... we were allowing forever.... filter by bounding box
            bb = [min([x1 x2]) max([x1 x2]); min([y1 y2]) max([y1 y2])];  
            % |xmin, xmax|
            % |ymin, ymax|
            indexes = ~(px>=bb(1,1) & px<=bb(1,2) & py>=bb(2,1) & py<=bb(2,2));
            dis(indexes,i,h) = NaN;
        end
        
    end

    dis(dis>mxd) = NaN;
    dis(abs(dir)>pi/4) = NaN;
    
    %% output
    dis=squeeze(nanmin(dis,[],2));
    dd=repmat(degs,size(rx,1),1) + repmat(hd,1,numel(degs));
    dx=dis.*cos(dd); dy=dis.*sin(dd);
    
    %% Ang and Dis to NaN if inside of the object
    in = inpolygon(rx,ry,edg(:,1),edg(:,2));
    dis(in,:) = NaN;
end

%%
function [dis, ang] = calcDisAng_Points(rx,ry,hd, envStruct)
    % Calculate distance and angle to nearest part of a circular object. If
    % animal is inside of the object, distance is NaN and the sample is
    % omitted from any future calculations
    envStruct = [envStruct.edg(1) envStruct.edg(2) envStruct.rad];
    
    % get closest point on circle to animal
    dx = (rx - envStruct(1));
    dy = (ry - envStruct(2));
    dis = sqrt(dx.^2 + dy.^2);
    ex = envStruct(3) .* [dx dy] ./ dis + [envStruct(1) envStruct(2)]; % intersection point
    ey = ex(:,2); 
    ex = ex(:,1);
    
    % Get distance and angle to surface point
    dx = rx-ex; 
    dy = ry-ey;
    ang = wrapToPi(hd - atan2(dy, dx)); % relative angle
    dis = sqrt(dx.^2+dy.^2);
    
    % NaN points inside
    inside = sqrt((rx - envStruct(1)).^2 + (ry- envStruct(2)).^2) <= envStruct(3);
    dis(inside) = NaN;
    ang(inside) = NaN;
    
end
%%

function [occ, nspk, rm_ns] = CalcMaps(dis, distanceBins, thetaBins, root, spikeMatch)
    occ = NaN(length(thetaBins), length(distanceBins));
    nspk = occ;
    distanceBins(end+1) = Inf;
    if iscell(root.ind)
        in = root.ind{1}(1);
    else
        in = root.ind(1); 
    end
    ci = CMBHOME.Utils.ContinuizeEpochs(root.cel_i) - in + 1;
    
    isIn = zeros(size(dis,1),1); isIn(CMBHOME.Utils.ContinuizeEpochs(root.ind)) = 1;
    isIn = logical(isIn);
    
    
    if ~isempty(spikeMatch)
        ci = sort(ci(randperm(length(ci),spikeMatch)));
    end
    
    for i = 1:length(thetaBins)
        t = dis(:,i);
        for k = 1:length(distanceBins)-1
            inds = t>=distanceBins(k) & t<distanceBins(k+1) & isIn;
            occ(i,k) = sum (inds);
            inds = find(inds);
            nspk(i,k) = length(intersect(inds,ci));
        end
    end
    distanceBins = distanceBins(1:end-1);
        
    % bring back to original dims
    occ = occ(:,1:end-1); occ=occ';
    nspk = nspk(:,1:end-1); nspk=nspk';
    rm_ns = (nspk./occ) * root.fs_video; % non-smoothed ratemap

end

%%
function [occ, nspk, rm_ns] = CalcMaps_Point(dis, ang, distanceBins, thetaBins, root, spikeMatch)

    if iscell(root.ind)
        in = root.ind{1}(1);
    else
        in = root.ind(1); 
    end
    ci = CMBHOME.Utils.ContinuizeEpochs(root.cel_i) - in + 1;
    
    isIn = zeros(size(dis,1),1); isIn(CMBHOME.Utils.ContinuizeEpochs(root.ind)) = 1;
    isIn = logical(isIn);
    
    if ~isempty(spikeMatch)
        ci = sort(ci(randperm(length(ci),spikeMatch)));
    end
    [occ] = CMBHOME.Utils.histcn([dis(isIn,:) ang(isIn,:)], distanceBins, thetaBins);
    [nspk] = CMBHOME.Utils.histcn([dis(ci) ang(ci)], distanceBins, thetaBins);
    
%     distCutOff = min(find(sum(occ,2)==0));
%     if ~isempty(distCutOff)
%         nspk = nspk(1:distCutOff,:);
%         occ = occ(1:distCutOff,:);
%     end
    
    rm_ns = (nspk./occ) * root.fs_video;
end

%%
function mp = smoothing(mp, k1, k2)
    nd = size(mp,2);
    mp = [mp mp mp];
    mp = CMBHOME.Utils.SmoothMat(mp, k1, k2);
    mp = mp(:, nd+1:2*nd);
end

% occ_sm = cellfun(@(x) smoothing(x, smooth(1:2), smooth(3)), occ,'UniformOutput',0);

function [MRL, PO, spikeDistMin, spikeDistMax]= ebcTests(nspk, occ,fs)

    import CMBHOME.Utils.circ.*
    
    % Restrict MRL calculation to range of distance bins in which a spike occurred.
    spikeDistMin = min(find(nansum(nspk,2)>0)); 
    spikeDistMax = max(find(nansum(nspk,2)>0));
        
    NFR = (nansum(nspk(spikeDistMin:spikeDistMax,:)) ./ nansum(occ(spikeDistMin:spikeDistMax,:))) * fs;
    thetaBins = deg2rad(linspace(-180,180,size(NFR,2)));

    MRL = circ_r(-pi:2*pi/(length(thetaBins)-1):pi, NFR, 2*pi/(length(thetaBins)-1), 2);

    PO = circ_mean(-pi:2*pi/(length(thetaBins)-1):pi, NFR, 2);

end





