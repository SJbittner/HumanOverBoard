clear variables;
close all;

%% setup

% Load Video
vidObj = VideoReader('MAH01462.MP4');
dt = 1 / vidObj.FrameRate; % time between frames
vidObj.CurrentTime = 0.4;

% load camera parameters
load cameraParams;

% output video
videoWriter = VideoWriter('out', 'Motion JPEG 2000');
videoWriter.FrameRate = vidObj.FrameRate;
open(videoWriter);

% canny edge detection
sigma1 = 0.0019;
sigma2 = 0.0016;

% structuring element for dilation
SE = [0 1 0; 1 1 1; 0 1 0];

% detect Buoy
initialDetectedBuoy = false;
buoyLost =false;
maxMoveBuoy = 0;

%for selecting the buoy manually
selectedBuoy = true;

% distances list
distances = [0];
% how many distances to average
dist_avg_K = 10;
% counting how many times the buoy is seen
buoy_seen_i = 0;

%% Read video frame by frame
while hasFrame(vidObj)

    % obtain video frame
    raw_frame = undistortImage(readFrame(vidObj), cameraParams);
    frame = im2gray(raw_frame);

    %% horizon detection

    % edge detection
    canny_ut = ut_edge(frame, 'c', 's', 3, 'h', [sigma1, sigma2]);
    canny_ut_gray = mat2gray(canny_ut);
    IM = canny_ut_gray;

    % dilate the horizontal edges
    for i=1:5
        IM = imdilate(IM, SE);
    end

    % Apply hough
    [H,T,R] = hough(IM,'RhoResolution',1,'Theta',-90:0.5:89);
    P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
    lines = houghlines(IM,T,R,P,'FillGap',5,'MinLength',7);

    % find longest straight line
    max_len = 0;
    for k = 1:length(lines)

        xy = [lines(k).point1; lines(k).point2];

        % Determine the endpoints of the longest line segment
        len = norm(lines(k).point1 - lines(k).point2);
        if (len > max_len)
            max_len = len;
            xy_long = xy;
        end
    end

    % Stabilize the frame by placing the horizon in the image center

    % interpolate linear function for the horizon segment we found
    dx = xy_long(2, 1) - xy_long(1, 1);
    dy = xy_long(2, 2) - xy_long(1, 2);
    x0 = xy_long(1, 1);
    y0 = xy_long(1, 2);
    c = -(dy / dx) * x0 + y0;

    % determine translation and rotation to stabilize the horizon
    [im_height, im_width] = size(IM);
    vp_y = im_height / 2;

    translationY = vp_y - ((im_width/2) * (dy/dx) + y0);
    theta = acosd(dy / dx) - 90;

    % apply first translation and then rotation to the image
    IM_translated = imtranslate(raw_frame, [0, translationY]);
    IM_stable = imrotate(IM_translated, -theta, 'bilinear', 'crop');

    %% buoy detection

    % smaller ROI cropping
    cropX = 668 - 50;
    cropY = 584 - 25;
    cropWidth = 200;
    cropHeight = 30;
    IM_cropped = imcrop(IM_stable, [cropX, cropY, cropWidth, cropHeight]);

    % thresholding, so only light pixels are tracked
    IM_cropped_mask = IM_cropped;
    light_threshold = 160;
    mask = im2gray(IM_cropped_mask) < light_threshold;
    IM_cropped_mask(mask) = 0;

    corners = detectHarrisFeatures(rgb2gray(IM_cropped_mask));

    % Select Buoy
    if selectedBuoy == false
        imshow(insertMarker(IM_cropped_mask, corners.Location, 'o'));
        try
            [xBuoy,yBuoy] = getpts;
            save buoyX xBuoy;
            save buoyY yBuoy;
            selectedBuoy = true;
        end
        continue; % sacrifice this frame; used for initialization
    end

    % Loading X and Y position of buoy in the begining
    if ~initialDetectedBuoy
        load buoyX;
        load buoyY;
    else
        xBuoy = buoyVector(1);
        yBuoy = buoyVector(2);
    end

    %% buoy tracking

    potentialCorners = corners.Location;
    differenceToBuoy = [];

    %calculating euclidean distances
    for i = 1:length(corners)
        XCorn = potentialCorners(i,1);
        YCorn = potentialCorners(i,2);
        differenceToBuoy(end+1) = [sqrt((xBuoy-XCorn)^2+(yBuoy-YCorn)^2)];
    end

    % thresholding the detection distance from last buoy coords
    buoyDistThresh = 11;

    % Getting closest distance point to buoy
    if length(corners) >= 1

        MinDist = min(differenceToBuoy);

        if MinDist > buoyDistThresh
            buoyLost = true;
            disp('Boy is drowning, oh noo!');
        else

            if buoyLost
                disp('Found the boii!');
            end

            % counting how many times it's seen
            buoy_seen_i = buoy_seen_i + 1;
            buoyLost = false;
            indexBuoy = find(differenceToBuoy==MinDist);
            buoyVector = potentialCorners(indexBuoy,:);
            initialDetectedBuoy = true;
            IM_cropped = insertMarker(IM_cropped, buoyVector , 'o');

            if MinDist > maxMoveBuoy
                maxMoveBuoy = MinDist;
            end
        end
    end

    % cropping for the buoy
    buoy_roi_dim = 40;
    IM_buoy_roi = imcrop(IM_cropped, [buoyVector-buoy_roi_dim/2, buoy_roi_dim, buoy_roi_dim]);

    %% estimate distance
    load cameraParams;

    % our camera focal lengths in pixels in x and y directions
    f_x = cameraParams.Intrinsics.FocalLength(1);
    f_y = cameraParams.Intrinsics.FocalLength(2);

    [im_height, im_width, dummy] = size(IM_stable);
    buoy_pos = buoyVector + [cropX, cropY]; % position of the buoy in the image
    cam_pos = [im_width / 2, im_height]; % position of the camera's vantage point in the image (center bottom)

    r = 6371; % globally-average radius of the earth in kilometers
    circumference_earth = 2 * r * pi; % spherical earth model
    h = 2.5 / 1000; % height of the camera above sea level in kilometers

    % determine length of straight line from camera to horizon
    l = sqrt((h + r)^2 - r^2);
    % determine angle between earth normal vector and the line to horizon
    theta = atan(r / l);
    beta = (pi / 2) - theta; % the other acute angle of the right triangle

    % distance from the camera to the horizon
    D_horizon = beta * (circumference_earth / (2 * pi));

    % linear function
    dx = cam_pos(2) - buoy_pos(2); % flip coordinate system to find horizon point
    dy = cam_pos(1) - buoy_pos(1);
    x0 = cam_pos(2);
    y0 = cam_pos(1);
    b = -(dy / dx) * x0 + y0;

    horizon_point = [(dy / dx) * (im_height / 2) + b, im_height / 2]; % unflipped coordinates

    % determine cartesian distance between the buoy and the horizon in pixels
    D_prime_buoy_horizon = norm(buoy_pos - horizon_point);

    % determine rate of increase of line through origin and our points
    a = dx / dy; % re-flip our flipped coordinates

    % determine focal length in direction of buoy (assume ellipse-shaped lens)
    a_prime = (1 / f_x^2) + (a^2 / f_y^2);
    f_xy_x = sqrt(4 * a_prime) / (2 * a_prime);
    f_xy_y = sqrt((1 - (f_xy_x^2 / f_x^2)) * f_y^2);
    f_xy = norm([f_xy_x, f_xy_y]);

    phi = atan(D_prime_buoy_horizon / f_xy);

    % linear function describing straight line from camera to buoy in real world
    a = -tan((pi / 2) - (theta - phi));
    c = r + h;

    % find intersection of that line and our spherical earth
    candidate_x_1 = (-(2*a*c) + sqrt((2*a*c)^2 - (4 * (a^2 + 1) * (c^2 - r^2)))) / (2 * (a^2 + 1));
    candidate_x_2 = (-(2*a*c) - sqrt((2*a*c)^2 - (4 * (a^2 + 1) * (c^2 - r^2)))) / (2 * (a^2 + 1));
    x0 = min(candidate_x_1, candidate_x_2);

    mu = asin(x0 / r);

    % We have our distance to the buoy in kilometers
    D_buoy = mu * (circumference_earth / (2 * pi));

    % output image
    IM_out = IM_stable;

    if ~buoyLost
        % store distance
        distances(end+1) = D_buoy;
        % take average of last K
        % this is for the first couple of frames, where there aren't
        % dist_avg_K distances in the list
        last_K = min(buoy_seen_i, dist_avg_K);
        D_buoy = mean(nonzeros(distances(end-last_K:end)));

        % draw line from here to buoy, labeled with distance
        IM_out = insertShape(IM_out, 'Line', [cam_pos(1) cam_pos(2), buoy_pos(1) buoy_pos(2)], 'LineWidth', 4, 'Color', 'yellow');
        text_pos = cam_pos + ([buoy_pos(1) - cam_pos(1), buoy_pos(2) - cam_pos(2)] / 2) + [20, -20];
        IM_out = insertText(IM_out, text_pos, round(D_buoy * 1000) + "m", 'FontSize', 32, 'TextColor', 'yellow', 'BoxOpacity', 0, 'AnchorPoint', 'LeftBottom');
    else
        % notify user that the buoy is lost
        IM_out = insertText(IM_out, cam_pos - [0, 20], "BUOY LOST", 'FontSize', 32, 'TextColor', 'red', 'BoxColor', 'yellow', 'AnchorPoint', 'CenterBottom');
    end

    % displaying final image ...
    imshow(IM_out);
    % ... with a red rectangle around the buoy
    rectangle('Position', ...
        [cropX+buoyVector(1)-buoy_roi_dim/2, ...
        cropY+buoyVector(2)-buoy_roi_dim/2, ...
        buoy_roi_dim buoy_roi_dim], ...
        'EdgeColor', 'r');
    hold on;
    % displaying the buoy ROI on top of original image
    image([100, 300], [100, 280], IM_buoy_roi);
    hold off;

    writeVideo(videoWriter, IM_out);

    pause(dt);
end

%% clean up
close(videoWriter);
