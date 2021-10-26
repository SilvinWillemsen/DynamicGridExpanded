if ~exist('firstPlot','var')
    firstPlot = true;
end
if ~exist('plotDampingAgainstFrequency','var')
    plotDampingAgainstFrequency = false;
end
if ~exist('noDamping','var')
    noDamping = false;
end
if ~exist('QisGiven','var')
    QisGiven = false;
end
if ~noDamping
    if firstPlot
        figure('Position', [173 578 827 220]) 
    end
    plot1Pos = [0.0632 0.2091 0.4115 0.6909];
else
    if firstPlot
        figure('Position', [173 578 435 220]);
    end

    plot1Pos = [0.1172 0.2091 0.8621 0.6909];
end


color = 'k';
% create Q matrix (one-step form)
if ~QisGiven
    Q = [Amat \ B, Amat \ C;
         speye(size(B)), sparse(zeros(size(B)))];
end
% obtain complex frequencies
s = 1/k * log(eig(full(Q)));

% obtain positive frequencies and sort them
s = s(imag(s) >= 0);
[~, order] = sort(imag(s));
s = s(order);

if ~noDamping
    subplot(121)
end
if firstPlot
    plot(0, 0)
end
hold on
scatter(1:length(s), imag(s)/(2*pi), color, 'Linewidth', 2)
grid on
xlabel("Mode number $p$", 'interpreter', 'latex')
ylabel("$f_p$ [Hz]", 'interpreter', 'latex')

xlim([1, length(s)])
title("Modal frequency")
set(gca, 'Fontsize', 16, 'Linewidth', 2, 'FontName', 'times', ...
        'Position', plot1Pos)
    
if ~noDamping
    subplot(122)
    if firstPlot
        color = 'k';
    else
        color = 'r';
    end
    if firstPlot
        plot(0, 0)
    end
    hold on
    if plotDampingAgainstFrequency
        scatter(imag(s)/(2*pi), real(s), color, 'Linewidth', 2)
        xlabel("Modal frequency [Hz]", 'interpreter', 'latex')

    else
        scatter(1:length(s), real(s), color, 'Linewidth', 2)
        xlabel("Mode number $p$", 'interpreter', 'latex')
        xlim([1, length(s)])

    end
    grid on
    % xlim([1, length(s)])
    title("Damping per mode")
    %     xlim([1, length(s)])
%
    ylabel("$\sigma_p$ [s$^{-1}$]", 'interpreter', 'latex')
    set(gca, 'Fontsize', 16, 'Linewidth', 2, 'FontName', 'times', ...
        'Position', [0.5732 0.2091 0.4115 0.6909])
end