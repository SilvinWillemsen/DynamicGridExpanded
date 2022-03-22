%to be used with modalAnalysis.m
modesSaveRange = 1:size (modesSave, 1);

figure('Position', [489 511 560 346])
plot(real(expectedF(1:n, :)) * 0.001, '--', 'color', [1, 0, 0, 1], 'Linewidth', 1)
hold on;
colorMap = zeros(length(modesSaveRange), 3);
normSigmaSave = abs(sigmaSave) / max(max(abs(sigmaSave)));
for scatLoop = 1:size(modesSave, 2)
    damp = exp(sigmaSave(:, scatLoop));
%     colorMap = repmat(1-damp, 1, 3);
    colorMap = [(1-damp) * 0.33,(1-damp) * 0.33, 1-damp];
    sz = damp*2 + 0.5;
%     sz = 5;
    for plotSection = 1:size(modesSave,1)-1
        if ~isnan(modesSave(plotSection+1, scatLoop))
            plot([plotSection, plotSection+1], real(modesSave(plotSection:plotSection+1, scatLoop))*0.001, ...
                'Linewidth', sz(plotSection), 'color', abs(colorMap(plotSection, :)));
        end
    end
%     scatter(modesSaveRange, real(modesSave(modesSaveRange, scatLoop))/1000, sz, colorMap, 'filled');
    hold on;
end
hold off
% title ("Modal Analysis $\mathcal{N}^n = " + N + " \rightarrow" + loopNend + "$", 'interpreter', 'latex');
xlabelsave = num2cell(Nstart:sign(Nend-Nstart):Nend);
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', nChangeSave, ... %loopStart(loopStartRange1)
    'xticklabel', xlabelsave, 'TickLabelInterpreter', 'latex', ...
    'Position', [0.0821 0.0809 0.8947 0.9135]);
labelX = xlabel("$\mathcal{N}^n$", 'interpreter', 'latex');
labelY = ylabel("Frequency (kHz)", 'interpreter', 'latex');
labelX.Position = [577.0325 -0.6538 -1.0000];
labelY.Position = [-39.0239 11.0250 -1.0000];
xlim([0, length(modesSaveRange)+1])
ylim([0, fs / 2000])
grid on