function  [varargout] = lpfFilter (varargin)

imy1=varargin{1};
fc = varargin{2};

% %1
% try zero-padding next time
% fc=6;
% fimy1=fft(imy1,length(imy1));
% fimy1(fc:end-fc+1)=0;
% my1=abs(ifft(fimy1));

%2
my1=sgolayfilt(imy1,4,9);

%3
% my1=lowpass(imy1,17,80,'ImpulseResponse','iir','Steepness',0.80);



%4
%   Wn=1000; %hz
% [b,a] = cheby2(15,20,Wn); % Example of ChebyII Filter
% my1 = filtfilt(b,a,imy1); % needed to create 0 phase offset
% 
%    

%5 polyfit
% x=1:length(imy1);
% p=polyfit(x,imy1,7);
% my1=polyval(p,x);

varargout{1}=my1;
