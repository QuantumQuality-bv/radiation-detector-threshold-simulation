% plot_detector_threshold_figures.m
% Recreate report-ready detector-threshold figures from generated CSV outputs.
% Run from the repository root or from the matlab/ directory.

clear; clc;

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
outputDir = fullfile(repoRoot, 'outputs');
figureDir = fullfile(repoRoot, 'figures');
imageDir = fullfile(repoRoot, 'images');
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end
if ~exist(imageDir, 'dir')
    mkdir(imageDir);
end

background = readtable(fullfile(outputDir, 'background_only_counts.csv'));
signal = readtable(fullfile(outputDir, 'signal_plus_background_counts.csv'));
summary = readtable(fullfile(outputDir, 'false_positive_summary.csv'));

backgroundMean = mean(background.observed_counts);
threshold3 = summary.threshold_counts(summary.threshold_sigma == 3);
threshold5 = summary.threshold_counts(summary.threshold_sigma == 5);

% Figure 1: background-only counts.
fig = newReportFigure();
plot(background.time_min, background.observed_counts, 'o-', ...
    'Color', [0.12 0.47 0.71], 'LineWidth', 1.8, 'MarkerSize', 4, ...
    'DisplayName', 'Observed counts');
hold on;
yline(backgroundMean, '--', 'Estimated background mean', ...
    'Color', [0.25 0.25 0.25], 'LineWidth', 1.8, 'DisplayName', 'Estimated background mean');
yline(threshold3, ':', '3-sigma threshold', ...
    'Color', [0.84 0.15 0.16], 'LineWidth', 2.0, 'DisplayName', '3-sigma threshold');
yline(threshold5, '-.', '5-sigma threshold', ...
    'Color', [0.58 0.40 0.74], 'LineWidth', 2.0, 'DisplayName', '5-sigma threshold');
title('Background-Only Synthetic Detector Counts', 'Color', 'k', 'FontWeight', 'bold');
xlabel('Time (minutes)');
ylabel('Observed counts per 60 s');
styleAxes(gca);
legend('Location', 'best', 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.7 0.7 0.7]);
saveFigure(fig, figureDir, imageDir, 'background_only_counts_matlab.png');

% Figure 2: signal-plus-background counts.
fig = newReportFigure();
plot(signal.time_min, signal.observed_counts, 'o-', ...
    'Color', [0.12 0.47 0.71], 'LineWidth', 1.8, 'MarkerSize', 4, ...
    'DisplayName', 'Observed counts');
hold on;
shadeEventIntervals(signal);
yline(threshold3, ':', '3-sigma threshold', ...
    'Color', [0.84 0.15 0.16], 'LineWidth', 2.0, 'DisplayName', '3-sigma threshold');
yline(threshold5, '-.', '5-sigma threshold', ...
    'Color', [0.58 0.40 0.74], 'LineWidth', 2.0, 'DisplayName', '5-sigma threshold');
title('Signal-Plus-Background Synthetic Detector Counts', 'Color', 'k', 'FontWeight', 'bold');
xlabel('Time (minutes)');
ylabel('Observed counts per 60 s');
styleAxes(gca);
legend('Location', 'best', 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.7 0.7 0.7]);
saveFigure(fig, figureDir, imageDir, 'signal_plus_background_counts_matlab.png');

% Figure 3: zoom around first elevated-count interval.
zoomMask = signal.time_min >= 38 & signal.time_min <= 62;
zoomData = signal(zoomMask, :);
fig = newReportFigure();
plot(zoomData.time_min, zoomData.observed_counts, 'o-', ...
    'Color', [0.12 0.47 0.71], 'LineWidth', 1.8, 'MarkerSize', 4, ...
    'DisplayName', 'Observed counts');
hold on;
shadeEventIntervals(zoomData);
yline(threshold3, ':', '3-sigma threshold', ...
    'Color', [0.84 0.15 0.16], 'LineWidth', 2.0, 'DisplayName', '3-sigma threshold');
yline(threshold5, '-.', '5-sigma threshold', ...
    'Color', [0.58 0.40 0.74], 'LineWidth', 2.0, 'DisplayName', '5-sigma threshold');
markFirstCrossing(zoomData, threshold3, [0.84 0.15 0.16], 'First 3-sigma crossing');
markFirstCrossing(zoomData, threshold5, [0.58 0.40 0.74], 'First 5-sigma crossing');
title('Zoom Around First Synthetic Elevated-Count Interval', 'Color', 'k', 'FontWeight', 'bold');
xlabel('Time (minutes)');
ylabel('Observed counts per 60 s');
styleAxes(gca);
legend('Location', 'best', 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.7 0.7 0.7]);
saveFigure(fig, figureDir, imageDir, 'threshold_crossing_zoom_matlab.png');

% Figure 4: empirical and analytical background threshold-crossing probabilities.
fig = newReportFigure();
x = summary.threshold_sigma;
empirical = summary.false_positive_fraction;
analyticalTail = summary.poisson_tail_probability;
atLeastOne = summary.probability_at_least_one_background_crossing_per_120;
positiveValues = [empirical(empirical > 0); analyticalTail(analyticalTail > 0); atLeastOne(atLeastOne > 0)];
visualFloor = max(min(positiveValues) / 5, 1e-8);
empiricalForPlot = empirical;
empiricalForPlot(empiricalForPlot <= 0) = visualFloor;
semilogy(x, empiricalForPlot, 'o-', ...
    'Color', [0.12 0.47 0.71], 'LineWidth', 1.8, 'MarkerSize', 6, ...
    'DisplayName', 'Empirical fraction (0 shown at floor)');
hold on;
semilogy(x, analyticalTail, 's-', ...
    'Color', [0.84 0.15 0.16], 'LineWidth', 1.8, 'MarkerSize', 6, ...
    'DisplayName', 'Analytical Poisson P(X > T_n)');
semilogy(x, atLeastOne, '^-', ...
    'Color', [0.17 0.63 0.17], 'LineWidth', 1.8, 'MarkerSize', 6, ...
    'DisplayName', 'P(at least one crossing in 120 bins)');
zeroRows = empirical == 0;
for k = find(zeroRows)'
    text(x(k), empiricalForPlot(k) * 1.7, '0 observed', ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.25 0.25 0.25]);
end
title({'Background Threshold-Crossing Probability', 'vs Threshold Multiplier'}, 'Color', 'k', 'FontWeight', 'bold');
xlabel('Threshold multiplier n');
ylabel('Background crossing probability / fraction');
xticks([2 3 4 5]);
ylim([visualFloor / 2, max(atLeastOne) * 2]);
styleAxes(gca);
legend('Location', 'best', 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.7 0.7 0.7]);
saveFigure(fig, figureDir, imageDir, 'false_positive_rate_vs_threshold_matlab.png');

disp('MATLAB figure export complete. Files written to figures/ and images/.');

function fig = newReportFigure()
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 900 540]);
    set(fig, 'InvertHardcopy', 'off');
end

function styleAxes(ax)
    set(ax, 'Color', 'w');
    set(ax, 'XColor', 'k', 'YColor', 'k');
    set(ax, 'FontSize', 11, 'LineWidth', 1.0);
    grid(ax, 'on');
    box(ax, 'on');
    ax.GridColor = [0.85 0.85 0.85];
    ax.GridAlpha = 1.0;
end

function shadeEventIntervals(tbl)
    eventMask = tbl.signal_expected_counts > 0;
    eventTimes = tbl.time_min(eventMask);
    if isempty(eventTimes)
        return;
    end
    yLimits = ylim;
    starts = eventTimes([true; diff(eventTimes) > 1]);
    ends = eventTimes([diff(eventTimes) > 1; true]) + 1;
    for idx = 1:numel(starts)
        patch([starts(idx) ends(idx) ends(idx) starts(idx)], ...
              [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
              [0.85 0.85 0.85], 'FaceAlpha', 0.45, 'EdgeColor', 'none', ...
              'DisplayName', 'Synthetic elevated interval');
    end
    lines = findobj(gca, 'Type', 'line');
    uistack(lines, 'top');
end

function markFirstCrossing(tbl, threshold, markerColor, labelText)
    crossingRows = find(tbl.observed_counts > threshold, 1, 'first');
    if isempty(crossingRows)
        return;
    end
    scatter(tbl.time_min(crossingRows), tbl.observed_counts(crossingRows), ...
        70, markerColor, 'filled', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.7, 'DisplayName', labelText);
end

function saveFigure(fig, figureDir, imageDir, fileName)
    figurePath = fullfile(figureDir, fileName);
    imagePath = fullfile(imageDir, fileName);
    exportgraphics(fig, figurePath, 'Resolution', 300, 'BackgroundColor', 'white');
    copyfile(figurePath, imagePath);
    close(fig);
end
