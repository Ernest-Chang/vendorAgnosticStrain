
function  [amp,phase,ny] = speckleFFT (varargin)


%step 0

nF=varargin{1};
my1=varargin{2};
K=fft(my1);


   
%step 1
%find fundamental frequency and filter
MK=max(K(3:end-2));
MKi=find(MK==K);


if length(varargin)>2
   fo=varargin{3};
   MKii=fo;
else
    
    MKii=MKi(1);
end
MKi2=nF-MKii+2;









%step 2
%find angle and amplitude
amp= abs(K(MKii))/nF*4; %peak to peak amplitude
phase=angle(K(MKii)); %radian


%step 3
%plot filtered waveform
nY=zeros(nF,1);
nY(MKii-1:MKii+1)=K(MKii-1:MKii+1);
nY(MKi2-1:MKi2+1)=K(MKi2-1:MKi2+1);
ny=ifft(nY);