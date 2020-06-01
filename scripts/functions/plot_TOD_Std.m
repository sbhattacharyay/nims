function [] = plot_TOD_Std(t,TOD_Std)
% STD by time of day figures:

%BandPower:
figure
subplot(7,1,1);
plot(t,TOD_Std{1,1});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,1});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,1});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,1});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,1});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,1});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,1});
datetick('x',15)
axis tight
suptitle('BandPower STD')

%Freq_Entropy

figure
subplot(7,1,1);
plot(t,TOD_Std{1,2});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,2});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,2});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,2});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,2});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,2});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,2});
datetick('x',15)
axis tight
suptitle('FreqEntropy STD')

%Freq_Pairs1

figure
subplot(7,1,1);
plot(t,TOD_Std{1,3});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,3});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,3});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,3});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,3});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,3});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,3});
datetick('x',15)
axis tight
suptitle('FreqPairs1 STD')

%Freq_Pairs2

figure
subplot(7,1,1);
plot(t,TOD_Std{1,4});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,4});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,4});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,4});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,4});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,4});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,4});
datetick('x',15)
axis tight
suptitle('FreqPairs2 STD')

%Med_Freq

figure
subplot(7,1,1);
plot(t,TOD_Std{1,5});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,5});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,5});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,5});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,5});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,5});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,5});
datetick('x',15)
axis tight
suptitle('MedFreq STD')

%SMA

figure
subplot(7,1,1);
plot(t,TOD_Std{1,6});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,6});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,6});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,6});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,6});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,6});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,6});
datetick('x',15)
axis tight
suptitle('SMA STD')

%Wavelets

figure
subplot(7,1,1);
plot(t,TOD_Std{1,7});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,7});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,7});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,7});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,7});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,7});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,7});
datetick('x',15)
axis tight
suptitle('Wavelet STD')
end
