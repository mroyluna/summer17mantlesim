function combined = trackStream(previousStreams, newStreams, n, lower, upper)
numStreams = length(newStreams);
xEndPoints = zeros(1,numStreams);
for ind = 1:numStreams
    xData = get(newStreams(ind), 'XData');
    xEndPoints(ind) = xData(end);
end

figure(10); clf; hold on;
subplot(1,2,1); hold on; title('Bin Dist at Current step');
hist(xEndPoints, n); xlim([lower,upper]);

combined = [previousStreams; histc(xEndPoints, linspace(lower,upper,n))]
subplot(1,2,2); hold on; title('Cumulative distribution');
area(combined'); xlim([1,n]); %ylim([0,100]);
end

%   function combTracers=trackStream(newStreams, n, lower, upper)
%         numStreams = length(newStreams);
%         xEndPoints = zeros(1,numStreams);
%         for index = 1:numStreams
%             xData = get(newStreams(index), 'XData');
%             xEndPoints(index) = xData(end);
%         end
%         figure(10); clf; hold on;
%         subplot(2,1,1); hold on; title('Bin Dist at Current step');
%         hist(xEndPoints, n); xlim([lower,upper]);
%         
%         combTracers = [combTracers; histc(xEndPoints, linspace(lower,upper,n))];
%         subplot(2,1,2); hold on; title('Cumulative distribution');
%         area(combTracers'); xlim([1,n]);
%     end