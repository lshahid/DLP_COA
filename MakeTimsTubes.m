% Make Tim's Tubes

%% Latex 1
L=130;
E=1.272;
h=0.3;
D=15.9;
A=pi*D^2/4;

segment_latex1=segment;

segment_latex1.E=10*E*10^6;
segment_latex1.h=h/10;
segment_latex1.data(2,3)=L;
segment_latex1.data(:,14)=[A; A];

%% Latex 2
L=130;
E=1.672;
h=0.3;
D=12.7;
A=pi*D^2/4;

segment_latex2=segment;

segment_latex2.E=10*E*10^6;
segment_latex2.h=h/10;
segment_latex2.data(2,3)=L;
segment_latex2.data(:,14)=[A, A];

%% Tygon
L=130;
E=2.065;
h=3.2;
D=19.1;
A=pi*D^2/4;

segment_tygon=segment;

segment_tygon.E=10*E*10^6;
segment_tygon.h=h/10;
segment_tygon.data(2,3)=L;
segment_tygon.data(:,14)=[A, A];

%% Silicone
L=130;
E=3.345;
h=3.2;
D=19.1;
A=pi*D^2/4;

segment_silicone=segment;

segment_silicone.E=10*E*10^6;
segment_silicone.h=h/10;
segment_silicone.data(2,3)=L;
segment_silicone.data(:,14)=[A, A];


