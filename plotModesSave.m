%to be used with modalAnalysisDynamicStiffString.m
% modesSave(modesSave<=0) = nan;

% end
% figure('Position', [489 511 560 346])
% figure('Position', [489 620 560 237])
figure ('Position', [440 394 516 234]); % for paper

plot(real(expectedF(1:n, :)) * 0.001, '--', 'color', [1, 0, 0, 1], 'Linewidth', 1)
hold on;

h = plot(real(modesSave) * 0.001, 'k', 'Linewidth', 1);

% if ~plotMulti
%     colours = [];
%     for colLoop = 1:floor(length(h))
%         if mod(colLoop,2) == 0
%             colours = [colours; 0,0,1];
%         else
%             colours = [colours; 1,0,0];
%         end
%     end
%     %     figure
%     set(h, {'color', 'Linewidth'}, [num2cell(colours, 2), num2cell(2 * ones(floor(length(h)), 1))])
% else
%     if interpolation == "quadratic" && lowPassConnection
%         set(h, {'color', 'Linewidth', 'Linestyle'}, [{[0.5, 0.5, 0.5]}, {1}, {'--'}]);
%     elseif interpolation == "quadratic"
%         set(h, {'color', 'Linewidth', 'Linestyle'}, [{[0.1, 0.1, 0.1]}, {1}, {'-'}]);
%     elseif interpolation == "sinc"
%         if fullSinc == 0
%             set(h, {'color', 'Linewidth', 'Linestyle'}, [{[0.8, 0.8, 0.8]}, {1}, {'-'}]);
%         elseif fullSinc == 1
%             set(h, {'color', 'Linewidth', 'Linestyle'}, [{[0.8, 0.8, 0.8]}, {1}, {'-'}]);
%         elseif fullSinc == 2
%             set(h, {'color', 'Linewidth', 'Linestyle'}, [{[0.6, 0.6, 0.6]}, {1}, {'-'}]);
%         elseif fullSinc == 3
%             set(h, {'color', 'Linewidth', 'Linestyle'}, [{[0.4, 0.4, 0.4]}, {1}, {'-'}]);
%         end
%     end
% end
title ("Modal Analysis $N = " + Nstart + " \rightarrow" + Nend + "$", 'interpreter', 'latex');
xlabelsave = num2cell(Nstart:sign(Nend-Nstart):Nend);
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', nChangeSave, ...
    'xticklabel', xlabelsave, 'TickLabelInterpreter', 'latex', ...
    'Position', [0.0950 0.1197 0.8876 0.8643]);
%     'Position', [0.0821 0.0809 0.8947 0.9135]);
labelX = xlabel("$\mathcal{N}^n$", 'interpreter', 'latex');
labelY = ylabel("Frequency (kHz)", 'interpreter', 'latex');
labelX.Position = [578.1680   -1.0568    1.0000];
labelY.Position = [-52.8035 11.0250 1.0000];

ylabel("Frequency (kHz)", 'interpreter', 'latex')
ylim([0, fs / 2000])
xlim([0, size(modesSave, 1)+1])
grid on