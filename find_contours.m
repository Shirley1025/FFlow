function [filtered_contours,x,y,x1,x2,y1,y2] = find_contours(Image1,threshold)
gray_img = rgb2gray(Image1);
binary_img = imbinarize(gray_img, threshold/255);  % Normalize the threshold to between [0,1]
filled_img = imfill(binary_img, 'holes');  % Filling holes

contours = bwboundaries(filled_img); % Find all eligible profiles
all_coords = [];

% Eliminate oversized or undersized contours
filtered_contours = {};
for i = 1:length(contours)
    area = size(contours{i}, 1);
    if area > 100 && area < 3000000
        filtered_contours{end+1} = contours{i};
    end
end

% Determine the edge pixel position of the contour
for k = 1:length(filtered_contours)
    inner_contour = filtered_contours{k};
    all_coords = [all_coords;inner_contour];  % Merge the coordinates of the current contour together    
end

% x = inner_contour(:, 2);
% y = inner_contour(:, 1);

x = all_coords(:, 2);
y = all_coords(:, 1);

% Get the edge position of the entire outline
x1=min(x, [], 'all')-10;
x2=max(x, [], 'all')+10;
y1=min(y, [], 'all')-10;
y2=max(y, [], 'all')+10;

