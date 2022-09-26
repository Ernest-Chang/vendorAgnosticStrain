
close all;


close all;
clear all;
[fileName,pathName]=uigetfile('*.*');
info =VideoReader([pathName fileName]);
nF=info.NumFrames;
Y=read(info,1);
vidFrame=Y;
% vidFrame=squeeze(Y(:,:,1)/3+Y(:,:,2)/3+Y(:,:,3)/3);


v=info;

currAxes = axes;
vidFrame1=readFrame(v);
imagesc(vidFrame1);

[y,x]=ginput(2);
y=ceil(y);
x=ceil(x);
z=improfile (vidFrame1,x,y);

zs=z(:,:,1);
zsp=[];
while hasFrame(v)
    vidFrame = readFrame(v);
    image(vidFrame, 'Parent', currAxes);
    z=improfile (vidFrame,y,x);

%     zs=z(:,:,1);
    zs=uint8(z);
    currAxes.Visible = 'off';
    pause(1/v.FrameRate);
    zsp=[zsp,zs];
end

figure;imagesc(zsp)