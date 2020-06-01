function waveOutput = get_wavelets(x,y,z)
[c1,l1] = wavedec(x,6,'db5');
[c2,l2] = wavedec(y,6,'db5');
[c3,l3] = wavedec(z,6,'db5');
xNorms = 0;
yNorms = 0;
zNorms = 0;

for p = 2:6
    xNorms=xNorms+norm(detcoef(c1,l1,p))^2;
    yNorms=yNorms+norm(detcoef(c2,l2,p))^2;
    zNorms=zNorms+norm(detcoef(c3,l3,p))^2;
end
waveOutput=rssq([xNorms,yNorms,zNorms]);
end
