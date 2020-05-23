function [] = plot_TOD_Diff(t,TOD_Diff)
% Diff by time of day figures:
t = t(1:end-1);

%BandPower:
figure
subplot(7,1,1);
plot(t,TOD_Diff{1,1});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,1});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,1});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,1});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,1});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,1});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,1});
datetick('x',15)
axis tight
suptitle('BandPower Diff')

%Freq_Entropy

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,2});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,2});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,2});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,2});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,2});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,2});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,2});
datetick('x',15)
axis tight
suptitle('FreqEntropy Diff')

%Freq_Pairs1

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,3});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,3});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,3});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,3});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,3});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,3});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,3});
datetick('x',15)
axis tight
suptitle('FreqPairs1 Diff')

%Freq_Pairs2

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,4});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,4});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,4});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,4});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,4});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,4});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,4});
datetick('x',15)
axis tight
suptitle('FreqPairs2 Diff')

%Med_Freq

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,5});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,5});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,5});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,5});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,5});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,5});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,5});
datetick('x',15)
axis tight
suptitle('MedFreq Diff')

%SMA

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,6});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,6});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,6});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,6});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,6});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,6});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,6});
datetick('x',15)
axis tight
suptitle('SMA Diff')

%Wavelets

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,7});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,7});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,7});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,7});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,7});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,7});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,7});
datetick('x',15)
axis tight
suptitle('Wavelet Diff')
end
