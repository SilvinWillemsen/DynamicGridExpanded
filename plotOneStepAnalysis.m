%to be used with modalAnalysis.m
modesSaveRange = 1:size (modesSave, 1);

colorMap = zeros(length(modesSaveRange), 3);
imagesc(exp(sigmaSave))
normSigmaSave = abs(sigmaSave) / max(max(abs(sigmaSave)));
for scatLoop = 1:size(modesSave, 2)
    damp = exp(sigmaSave(:, scatLoop));
%     colorMap = repmat(1-damp, 1, 3);
    colorMap = [(1-damp) * 0.33,(1-damp) * 0.33, 1-damp];
    sz = 7.5 * damp + 2;
%     sz = 5;
    scatter(modesSaveRange, real(modesSave(modesSaveRange, scatLoop))/1000, sz, colorMap, 'filled');
    hold on;
end
hold off
% title ("Modal Analysis $\mathcal{N}^n = " + N + " \rightarrow" + loopNend + "$", 'interpreter', 'latex');
xlabelsave = num2cell(NinitSave:sign(Nend-NinitSave):Nend);
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', loopStart(loopStartRange1(1:end-1)), ... %loopStart(loopStartRange1)
    'xticklabel', xlabelsave, 'TickLabelInterpreter', 'latex', ...
    'Position', [0.0969, 0.1100, 0.8799, 0.7869]);
labelX = xlabel("$\mathcal{N}^n$", 'interpreter', 'latex');
ylabel("Frequency (kHz)", 'interpreter', 'latex');
labelX.Position(2) = -1;
labelX.Position(1) = 0.575 * size(modesSave, 1);
xlim([0, length(loopAmountRange)])
ylim([0, fs / 2000])
grid on