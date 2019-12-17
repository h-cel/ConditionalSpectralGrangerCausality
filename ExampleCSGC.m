%%Example of how condtional spectral Granger causality (CSGC) is determined
%%for 1 location

% load data
load('DATA.mat')
Variables = fields(DATA);
N = size(DATA.(Variables{1}), 1);
R = size(DATA.(Variables{1}), 2);

for i1 = 1:length(Variables)
    X{i1} = DATA.(Variables{i1});
end

% calculate CSGC
[CSGC, Fr] =  CalculateCSGC(X, 'param', 6, 'order', 1);
Scale = 1./(Fr*12);

% Plot spectrum
Colours = {'b', 'r', 'g'};
figure()
for i1 = 1:(length(Variables)-1)
    plot(Scale, (1-exp(-CSGC{1, i1+1}))*100, [Colours{i1}, '-'], 'linewidth', 2)
    if i1 == 1
        hold on
    end
end
set(gca, 'xscale', 'log', 'xtick', [0.1, 1, 10],...
    'xticklabel', {'0.1', '1', '10'},...
    'ylim', [0, 100], 'ytick', 0:20:100)
xlabel('Scale [Years]')
ylabel('Percentage of variance in LAI explained [%]')
legend('Precipitation', 'Air temperature', 'Net radiation')

% Plot fraction of variance explained for scales of interest
SOI = {'Monthly', 'Seasonal', 'Interannual'};
Scales_of_interest.Monthly = [0, 0.32];
Scales_of_interest.Seasonal = [0.32, 1.54];
Scales_of_interest.Interannual = [1.54, 9];

for i1 = 1:length(Variables)
    for i2 = 1:length(Variables)
        if (i1~=i2) && ((i1 == 1) || (i2 == 1))
            for i3 = 1:length(SOI)
                Filter = (Scale > Scales_of_interest.(SOI{i3})(1)) &...
                    (Scale <= Scales_of_interest.(SOI{i3})(2));
                
                Cause.(Variables{i1}).To.(Variables{i2}).(SOI{i3}) =...
                    nanmean((1-exp(-CSGC{i2, i1}(Filter)))*100);
            end
        end
    end
end

figure()
for i3 = 1:3
    subplot(1, 3, i3)
    for i2 = 1:(length(Variables)-1)
        fill(i2 - [-0.3, -0.3, 0.3, 0.3],...
            [0, ones(1, 2)*Cause.(Variables{i2+1}).To.LAI.(SOI{i3}), 0],...
            Colours{i2}, 'edgecolor', 'none')
        if i2 == 1
            hold on
        end
    end
    set(gca, 'xlim', [0.5, 3.5], 'xtick', 1:3, 'xticklabels',...
        Variables(2:end), 'ylim', [0, 40], 'ytick', 0:10:40)
    ylabel('Percentage of variance in LAI explained [%]')
    title(SOI{i3})
end
        
    






