[width,height]=ERP_visualangle([340 340].*2,120,[1920 1080],[63.5 35])

[ visang ] = pix2visang(390,63.5,1920,120)./2
[ visang ] = pix2visang(390,63.5,1920,120)./2


% 360 = 11.2°
% 340 = 10.4° (as PNAS)

% pixels:
% 12 = 0.38°
% 8 = 0.25° (as PNAS, 120 pr RDK)


[width,height]=ERP_visualangle([360 360].*2,120,[1920 1080],[63.5 35])
[width,height]=ERP_visualangle([12 12].*2,120,[1920 1080],[63.5 35])
[width,height]=ERP_visualangle([1 1].*2,120,[1920 1080],[63.5 35])
