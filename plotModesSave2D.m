%to be used with modalAnalysis.m
% modesSave(modesSave<=0) = nan;
% if (Nend-Nstart) < 0
%     modesSaveRange = 2:size(modesSave, 1);
%     loopStartRange1 = 2:length(loopStart)-1;
% else
%     modesSaveRange = 1:size(modesSave,1);
%     loopStartRange1 = 1:length(loopStart)-1;
% 
% end
% if interpolation == "quadratic" && plotMulti && ~lowPassConnection 
%     hold on;
% else
% end
% if interpolation == "linear" && plotMulti && fullSinc ~= 0
%     hold on;
% end
figure('Position', [489 511 560 346])

% plot(real(expectedF(1:n, :)) * 0.001, '--', 'color', [1, 0, 0, 1], 'Linewidth', 1)
hold on;

for i = 1:length(nChangeSave)-1
    plot(nChangeSave(i):nChangeSave(i+1)-1, modesSave(nChangeSave(i):nChangeSave(i+1)-1, :) * 0.001, 'k', 'Linewidth', 1)
    hold on;
    plot(nChangeSave(i):nChangeSave(i+1)-1, expectedF(nChangeSave(i):nChangeSave(i+1)-1, :) * 0.001, '--r' ,'Linewidth', 1)
    hold on;

end
% h = plot(real(modesSave) * 0.001, 'k', 'Linewidth', 1);

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
if NxStart - NxEnd ~= 0
    Nstart = NxStart;
    Nend = NxEnd;
else
    Nstart = NyStart;
    Nend = NyEnd;
end

title ("Modal Analysis $N = " + Nstart + " \rightarrow" + Nend + "$", 'interpreter', 'latex');
xlabelsave = num2cell(Nstart:sign(Nend-Nstart):Nend);
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', nChangeSave, ...
    'xticklabel', xlabelsave, 'TickLabelInterpreter', 'latex', ...  
    'Position', [0.0821 0.0809 0.8947 0.9135]);
labelX = xlabel("$\mathcal{N}_x^n$", 'interpreter', 'latex');
labelY = ylabel("Frequency (kHz)", 'interpreter', 'latex');
labelX.Position = [537.0724 -0.4445 -1.0000];
labelY.Position = [-39.0239 11.0250 -1.0000];

ylabel("Frequency (kHz)", 'interpreter', 'latex')
ylim([0, fs / 2000])
xlim([0, size(modesSave, 1)+1])
grid on