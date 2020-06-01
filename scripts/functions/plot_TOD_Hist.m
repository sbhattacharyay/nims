function [] = plot_TOD_Hist(TOD_Means)

%BandPower:
figure
subplot(7,1,1);
histogram(TOD_Means{1,1});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,1});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,1});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,1});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,1});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,1});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,1});

axis tight
suptitle('BandPower Means Hist')

%Freq_Entropy

figure
subplot(7,1,1);
histogram(TOD_Means{1,2});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,2});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,2});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,2});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,2});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,2});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,2});

axis tight
suptitle('FreqEntropy Means Hist')

%Freq_Pairs1

figure
subplot(7,1,1);
histogram(TOD_Means{1,3});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,3});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,3});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,3});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,3});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,3});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,3});

axis tight
suptitle('FreqPairs1 Means Hist')

%Freq_Pairs2

figure
subplot(7,1,1);
histogram(TOD_Means{1,4});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,4});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,4});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,4});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,4});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,4});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,4});

axis tight
suptitle('FreqPairs2 Means Hist')

%Med_Freq

figure
subplot(7,1,1);
histogram(TOD_Means{1,5});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,5});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,5});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,5});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,5});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,5});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,5});

axis tight
suptitle('MedFreq Means Hist')

%SMA

figure
subplot(7,1,1);
histogram(TOD_Means{1,6});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,6});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,6});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,6});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,6});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,6});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,6});

axis tight
suptitle('SMA Means Hist')

%Wavelets

figure
subplot(7,1,1);
histogram(TOD_Means{1,7});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,7});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,7});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,7});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,7});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,7});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,7});

axis tight
suptitle('Wavelet Means Hist')

end
