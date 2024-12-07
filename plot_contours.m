figure;
imshow(Image1);
hold on;
for k = 1:length(filtered_contours)
    contour = filtered_contours{k};
    plot(contour(:, 2), contour(:, 1), 'g', 'LineWidth', 2);
end

% Calculate the width and height of a rectangle
width = x2 - x1;
height = y2 - y1;
% Drawing rectangles on an image
rectangle('Position', [x1, y1, width, height], 'EdgeColor', 'r', 'LineWidth', 2);
title('All Contours in Image');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
set(get(gca, 'Title'), 'FontName', 'Times New Roman', 'FontSize', 14);
hold off;