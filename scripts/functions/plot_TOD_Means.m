function [] = plot_TOD_Means(t,TOD_Means)
%Mean by time of day figures:

%BandPower:
figure
subplot(7,1,1);
plot(t,TOD_Means{1,1});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,1});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,1});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,1});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,1});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,1});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,1});
datetick('x',15)
axis tight
suptitle('BandPower Means')

%Freq_Entropy

figure
subplot(7,1,1);
plot(t,TOD_Means{1,2});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,2});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,2});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,2});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,2});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,2});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,2});
datetick('x',15)
axis tight
suptitle('FreqEntropy Means')

%Freq_Pairs1

figure
subplot(7,1,1);
plot(t,TOD_Means{1,3});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,3});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,3});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,3});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,3});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,3});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,3});
datetick('x',15)
axis tight
suptitle('FreqPairs1 Means')

%Freq_Pairs2

figure
subplot(7,1,1);
plot(t,TOD_Means{1,4});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,4});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,4});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,4});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,4});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,4});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,4});
datetick('x',15)
axis tight
suptitle('FreqPairs2 Means')

%Med_Freq

figure
subplot(7,1,1);
plot(t,TOD_Means{1,5});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,5});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,5});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,5});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,5});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,5});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,5});
datetick('x',15)
axis tight
suptitle('MedFreq Means')

%SMA

figure
subplot(7,1,1);
plot(t,TOD_Means{1,6});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,6});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,6});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,6});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,6});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,6});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,6});
datetick('x',15)
axis tight
suptitle('SMA Means')

%Wavelets

figure
subplot(7,1,1);
plot(t,TOD_Means{1,7});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,7});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,7});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,7});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,7});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,7});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,7});
datetick('x',15)
axis tight
suptitle('Wavelet Means')
end

