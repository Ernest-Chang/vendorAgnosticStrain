function  [nnvidFrame2,nnvidFrame3,nF,mStrain,cycleL]= speckle_strain (init,fin,nF,data,pathName,fileName)


% close all;
mymap=[
             0         0         0
    0.1250         0         0
    0.2500         0         0
    0.3750         0         0
    0.5000         0         0
    0.6250         0         0
    0.7500         0         0
    0.8750         0         0
    1.0000         0         0
    1.0000    0.0625         0
    1.0000    0.1250         0
    1.0000    0.1875         0
    1.0000    0.2500         0
    1.0000    0.3125         0
    1.0000    0.3750         0
    1.0000    0.4375         0
    1.0000    0.5000         0
    1.0000    0.5625         0
    1.0000    0.6250         0
    1.0000    0.6875         0
    1.0000    0.7500         0
    1.0000    0.8125         0
    1.0000    0.8750         0
    1.0000    0.9375         0
    1.0000    1.0000         0
    0.9375    1.0000    0.0625
    0.8750    1.0000    0.1250
    0.8125    1.0000    0.1875
    0.7500    1.0000    0.2500
    0.6875    1.0000    0.3125
    0.6250    1.0000    0.3750
    0.5625    1.0000    0.4375
    0.5000    1.0000    0.5000
    0.4375    1.0000    0.5625
    0.3750    1.0000    0.6250
    0.3125    1.0000    0.6875
    0.2500    1.0000    0.7500
    0.1875    1.0000    0.8125
    0.1250    1.0000    0.8750
    0.0625    1.0000    0.9375
         0    1.0000    1.0000
         0    0.9375    1.0000
         0    0.8750    1.0000
         0    0.8125    1.0000
         0    0.7500    1.0000
         0    0.6875    1.0000
         0    0.6250    1.0000
         0    0.5625    1.0000
         0    0.5000    1.0000
         0    0.4375    1.0000
         0    0.3750    1.0000
         0    0.3125    1.0000
         0    0.2500    1.0000
         0    0.1875    1.0000
         0    0.1250    1.0000
         0    0.0625    1.0000
         0         0    1.0000
         0         0    0.9375
         0         0    0.8750
         0         0    0.8125
         0         0    0.7500
         0         0    0.6875
         0         0    0.6250
         0         0    0.5625];

% [fileName,pathName]=uigetfile('*.*');
% info =VideoReader([pathName fileName]);
% nF=nF

% init=1;
% init=140;clc
vidFrame=data(:,:,init);

% for EchoNet data
% vidFrame=imresize(vidFrame,6);

% --------- get diastolic frame-------------
% for i=15:50
%     figure;
%     Y=read(info,i);
%     imagesc(Y)
% end
% adf

% vidFrame1=medfilt2(squeeze(vidFrame(:,:,1)),[5 5]);

vidfile = VideoWriter([pathName fileName '1'],'MPEG-4');
open(vidfile);


f=figure;
imagesc(vidFrame),colormap(gray);
h=drawfreehand;
xs=h.Position(:,1);
ys=h.Position(:,2);
% close(f);

% [xs,ys]=snake_read(vidFrame);
xs=round(xs);
xs(end)=[];
ys=round(ys);
ys(end)=[];


n=12;

% even spacing and interpolation
xf=linspace(1,length(xs),n);

p=polyfit(1:length(xs),xs',5);
xf2=round(p(1)*xf.^5+p(2)*xf.^4+p(3)*xf.^3+p(4)*xf.^2+p(5)*xf.^1+p(6));
% xf2(end)=[];

yf=linspace(1,length(ys),n);
p=polyfit(1:length(ys),ys',5);
yf2=round(p(1)*yf.^5+p(2)*yf.^4+p(3)*yf.^3+p(4)*yf.^2+p(5)*yf.^1+p(6));

bs=max(xs)-min(xs);
bs2=max(ys)-min(ys);
if bs>bs2
    msz=bs;
else
    msz=bs2;
end

msz=ceil(msz/n)*4;


% fprintf('msz is %02d \n',msz);

% xf2=linspace(min(xs),max(xs),n);
% p=polyfit(xs,ys,);
% yf2=polyval(p,xf2);
% 
% figure;plot(xs,ys);
% figure;plot(xf2,yf2);
% afafasdf

% yf2(end)=[];

% xf3=resample(xs(1:end-1),20,length(xs)-1);
% yf3=resample(ys(1:end-1),20,length(ys)-1);
% 

[s1,s2]=size(xf2);
xf3=zeros(s1,s2);
yf3=xf3;
xf3(1)=xf2(1);


% figure; plot(xs,ys,xf2,yf2,'ro');

[a,b,f]=size(vidFrame);
ytt=zeros(f,nF);
xtt=ytt;
yt=[];
xt=[];

    vidFrame2=squeeze(vidFrame(:,:,1));
%     imagesc(vidFrame2);
    ix=0;
    iy=0;
% vidFrame=medfilt2(vidFrame, [4 4]);
% vidFrame=vidFrame-mean(mean(vidFrame));

% iBox=28;
% iBox=15;
% iFrame=40;

% 

if msz>100
    iFrame=80;
    iBox=50;
else
    iFrame=ceil(msz);
    iBox=ceil(msz*0.6);
end
% iBox=45;
% iFrame=70;


% mouse studies
% iBox=70
% iFrame=200;

%------Gaussian PDF
n=(iBox+iFrame)*2+1;
x=1:n;
y=exp(-(x-floor(n/2)).^2/(iFrame)^2); %human
% y=exp(-(x-floor(n/2)).^2/(iFrame/1.8)^2); %experimental

y=exp(-(x-floor(n/2)).^2/(iFrame/2)^2); %mice
% y=exp(-(x-floor(n/2)).^2/(iFrame/4)^2); %mice

z= ones(n,n);
y1=repmat(y,n,1);
y2=repmat(y',1,n);
zg=(z.*y1.*y2);
% *.9+0.1;
%-------------------%

crBox=iBox+iFrame;
vidFrame=medfilt2(vidFrame,[2 2]);


% for i=1


%this can go into speckleTrack function
    box1 = vidFrame(yf2(1)-iBox:yf2(1)+iBox,xf2(1)-iBox:xf2(1)+iBox);
    box2 = vidFrame(yf2(2)-iBox:yf2(2)+iBox,xf2(2)-iBox:xf2(2)+iBox);
    box3 = vidFrame(yf2(3)-iBox:yf2(3)+iBox,xf2(3)-iBox:xf2(3)+iBox);
    box4 = vidFrame(yf2(4)-iBox:yf2(4)+iBox,xf2(4)-iBox:xf2(4)+iBox);
    box5 = vidFrame(yf2(5)-iBox:yf2(5)+iBox,xf2(5)-iBox:xf2(5)+iBox);
    box6 = vidFrame(yf2(6)-iBox:yf2(6)+iBox,xf2(6)-iBox:xf2(6)+iBox);
    box7 = vidFrame(yf2(7)-iBox:yf2(7)+iBox,xf2(7)-iBox:xf2(7)+iBox);
    
    box8 = vidFrame(yf2(8)-iBox:yf2(8)+iBox,xf2(8)-iBox:xf2(8)+iBox);
    box9 = vidFrame(yf2(9)-iBox:yf2(9)+iBox,xf2(9)-iBox:xf2(9)+iBox);
    box10 = vidFrame(yf2(10)-iBox:yf2(10)+iBox,xf2(10)-iBox:xf2(10)+iBox);
    box11 = vidFrame(yf2(11)-iBox:yf2(11)+iBox,xf2(11)-iBox:xf2(11)+iBox);
    box12 = vidFrame(yf2(12)-iBox:yf2(12)+iBox,xf2(12)-iBox:xf2(12)+iBox);

    

%     box1=box1-mean(mean(box1));
%     box1= vidFrame1(ys(i)-10:ys(i)+10,xs(i)-10:xs(i)+10);
yt=[];
xt=[];

iy=0; ix=0;

% v= VideoReader('A4C.mov');
% figure
% imagesc(vidFrame)
% hold on
% plot(xf2(1)+iy,yf2(1)+ix,'rs')
% hold off
% i=1;



%initialization
ny1=yf2(1);
nx1=xf2(1);
ny2=yf2(2);
nx2=xf2(2);
ny3=yf2(3);
nx3=xf2(3);
ny4=yf2(4);
nx4=xf2(4);
ny5=yf2(5);
nx5=xf2(5);
ny6=yf2(6);
nx6=xf2(6);
ny7=yf2(7);
nx7=xf2(7);
ny8=yf2(8);
nx8=xf2(8);
ny9=yf2(9);
nx9=xf2(9);
ny10=yf2(10);
nx10=xf2(10);
ny11=yf2(11);
nx11=xf2(11);
ny12=yf2(12);
nx12=xf2(12);



mat1x=ones(nF,1)*nx1;
mat2x=ones(nF,1)*nx2;
mat3x=ones(nF,1)*nx3;
mat4x=ones(nF,1)*nx4;
mat5x=ones(nF,1)*nx5;
mat6x=ones(nF,1)*nx6;
mat7x=ones(nF,1)*nx7;
mat8x=ones(nF,1)*nx8;
mat9x=ones(nF,1)*nx9;
mat10x=ones(nF,1)*nx10;
mat11x=ones(nF,1)*nx11;
mat12x=ones(nF,1)*nx12;

Mmatx=[mat1x';mat2x';mat3x';mat4x';mat5x';mat6x';mat7x';mat8x';mat9x';mat10x';mat11x';mat12x'];

mat1y=ones(nF,1)*ny1;
mat2y=ones(nF,1)*ny2;
mat3y=ones(nF,1)*ny3;
mat4y=ones(nF,1)*ny4;
mat5y=ones(nF,1)*ny5;
mat6y=ones(nF,1)*ny6;
mat7y=ones(nF,1)*ny7;
mat8y=ones(nF,1)*ny8;
mat9y=ones(nF,1)*ny9;
mat10y=ones(nF,1)*ny10;
mat11y=ones(nF,1)*ny11;
mat12y=ones(nF,1)*ny12;
Mmaty=[mat1y';mat2y';mat3y';mat4y';mat5y';mat6y';mat7y';mat8y';mat9y';mat10y';mat11y';mat12y'];


ix1=0;ix2=ix1;ix3=ix1;ix4=ix1;ix5=ix1;ix6=ix1;ix7=ix1;ix8=ix1;ix9=ix1;ix10=ix1;ix11=ix1;ix12=ix1;
iy1=0;iy2=iy1;iy3=iy1;iy4=iy1;iy5=iy1;iy6=iy1;iy7=iy1;iy8=iy1;iy9=iy1;iy10=iy1;iy11=iy1;iy12=iy1;

if fin==nF
    fin=nF-1;
end


g=figure;

for v= init:fin
%     fprintf('%2d iteration',v);


    vidFrame=data(:,:,v+1);


    % normal strain
    vidFrameMinus=data(:,:,v);
    %--------------------------

    %  EchoNet Data--------------
%     vidFrame=imresize(vidFrame,6);
%     vidFrameMinus=data(:,:,v);
%     vidFrameMinus=imresize(vidFrameMinus,6);
    %-------------------------------------


    
%     vidFrame=squeeze(frame(:,:,1)/3+frame(:,:,2)/3+frame(:,:,3)/3);
    vidFrame=medfilt2(vidFrame,[2 2]);
    vidFrameMinus=medfilt2(vidFrameMinus,[2 2]);

    
    %clause to update iBox information after each iteration
    
%     %trial 1
% if v==init
%     continue
% else

    [ix1, iy1]= speckleTrack (vidFrame,yf2(1)+ix1,xf2(1)+iy1,iFrame,box1,crBox,iBox,zg);
    [ix2, iy2]= speckleTrack (vidFrame,yf2(2)+ix2,xf2(2)+iy2,iFrame,box2,crBox,iBox,zg);
    [ix3, iy3]= speckleTrack (vidFrame,yf2(3)+ix3,xf2(3)+iy3,iFrame,box3,crBox,iBox,zg);
    [ix4, iy4]= speckleTrack (vidFrame,yf2(4)+ix4,xf2(4)+iy4,iFrame,box4,crBox,iBox,zg);
    [ix5, iy5]= speckleTrack (vidFrame,yf2(5)+ix5,xf2(5)+iy5,iFrame,box5,crBox,iBox,zg);
    [ix6, iy6]= speckleTrack (vidFrame,yf2(6)+ix6,xf2(6)+iy6,iFrame,box6,crBox,iBox,zg);
    [ix7, iy7]= speckleTrack (vidFrame,yf2(7)+ix7,xf2(7)+iy7,iFrame,box7,crBox,iBox,zg);
    [ix8, iy8]= speckleTrack (vidFrame,yf2(8)+ix8,xf2(8)+iy8,iFrame,box8,crBox,iBox,zg);
    [ix9, iy9]= speckleTrack (vidFrame,yf2(9)+ix9,xf2(9)+iy9,iFrame,box9,crBox,iBox,zg);
    [ix10, iy10]= speckleTrack (vidFrame,yf2(10)+ix10,xf2(10)+iy10,iFrame,box10,crBox,iBox,zg);
    [ix11, iy11]= speckleTrack (vidFrame,yf2(11)+ix11,xf2(11)+iy11,iFrame,box11,crBox,iBox,zg);
    [ix12, iy12]= speckleTrack (vidFrame,yf2(12)+ix12,xf2(12)+iy12,iFrame,box12,crBox,iBox,zg);

    

% end

% trial 2

    ix=[ix1,ix2,ix3,ix4,ix5,ix6,ix7,ix8,ix9,ix10,ix11,ix12];
    iy=[iy1,iy2,iy3,iy4,iy5,iy6,iy7,iy8,iy9,iy10,iy11,iy12];
    

    box1 = Boxing(vidFrameMinus, yf2(1),xf2(1), iBox);
    box2 = Boxing(vidFrameMinus, yf2(2),xf2(2), iBox);
    box3 = Boxing(vidFrameMinus, yf2(3),xf2(3), iBox);
    box4 = Boxing(vidFrameMinus, yf2(4),xf2(4), iBox);
    box5 = Boxing(vidFrameMinus, yf2(5),xf2(5), iBox);
    box6 = Boxing(vidFrameMinus, yf2(6),xf2(6), iBox);
    box7 = Boxing(vidFrameMinus, yf2(7),xf2(7), iBox);
    box8 = Boxing(vidFrameMinus, yf2(8),xf2(8), iBox);
    box9 = Boxing(vidFrameMinus, yf2(9), xf2(9), iBox);
    box10= Boxing(vidFrameMinus, yf2(10),xf2(10), iBox);
    box11= Boxing(vidFrameMinus, yf2(11),xf2(11), iBox);
%     box12= Boxing(vidFrameMinus, yf2(12),xf2(12), iBox);
    
         [xf2,yf2]=UpdateLine(xf2,yf2,ix,iy);


    [ix1, iy1]= speckleTrack (vidFrame,yf2(1),xf2(1),iFrame,box1,crBox,iBox,zg);
    [ix2, iy2]= speckleTrack (vidFrame,yf2(2),xf2(2),iFrame,box2,crBox,iBox,zg);
    [ix3, iy3]= speckleTrack (vidFrame,yf2(3),xf2(3),iFrame,box3,crBox,iBox,zg);
    [ix4, iy4]= speckleTrack (vidFrame,yf2(4),xf2(4),iFrame,box4,crBox,iBox,zg);
    [ix5, iy5]= speckleTrack (vidFrame,yf2(5),xf2(5),iFrame,box5,crBox,iBox,zg);
    [ix6, iy6]= speckleTrack (vidFrame,yf2(6),xf2(6),iFrame,box6,crBox,iBox,zg);
    [ix7, iy7]= speckleTrack (vidFrame,yf2(7),xf2(7),iFrame,box7,crBox,iBox,zg);
    [ix8, iy8]= speckleTrack (vidFrame,yf2(8),xf2(8),iFrame,box8,crBox,iBox,zg);
    [ix9, iy9]= speckleTrack (vidFrame,yf2(9),xf2(9),iFrame,box9,crBox,iBox,zg);
    [ix10, iy10]= speckleTrack (vidFrame,yf2(10),xf2(10),iFrame,box10,crBox,iBox,zg);
    [ix11, iy11]= speckleTrack (vidFrame,yf2(11),xf2(11),iFrame,box11,crBox,iBox,zg);
    [ix12, iy12]= speckleTrack (vidFrame,yf2(12),xf2(12),iFrame,box12,crBox,iBox,zg);
 

%              [xf2,yf2]=UpdateLine(xf2,yf2,ix,iy);

%     
%     


    % Adaptive reference point box shifting for next iteration
%     ny1=ny1+ix1;
%     nx1=nx1+iy1;
%     ny2=ny2+ix2;
%     nx2=nx2+iy2;
%     ny3=ny3+ix3;
%     nx3=nx3+iy3;
%     ny4=ny4+ix4;
%     nx4=nx4+iy4;
%     ny5=ny5+ix5;
%     nx5=nx5+iy5;
%     ny6=ny6+ix6;
%     nx6=nx6+iy6;
%     ny7=ny7+ix7;
%     nx7=nx7+iy7;
% 
%     
%     box1 = vidFrame(ny1-iBox:ny1+iBox,nx1-iBox:nx1+iBox);
%     box2 = vidFrame(ny2-iBox:ny2+iBox,nx2-iBox:nx2+iBox);
%     box3 = vidFrame(ny3-iBox:ny3+iBox,nx3-iBox:nx3+iBox);
%     box4 = vidFrame(ny4-iBox:ny4+iBox,nx4-iBox:nx4+iBox);
%     box5 = vidFrame(ny5-iBox:ny5+iBox,nx5-iBox:nx5+iBox);
%     box6 = vidFrame(ny6-iBox:ny6+iBox,nx6-iBox:nx6+iBox);
%     box7 = vidFrame(ny7-iBox:ny7+iBox,nx7-iBox:nx7+iBox);
%     
%     xx=[nx1, nx2, nx3, nx4, nx5, nx6, nx7];
%     yy=[ny1, ny2, ny3, ny4, ny5, ny6, ny7];
    
    %-------------------------------------------------------%

    
    
% ------- sourcecode for speckleTrack --------------------------%
%     frame2=vidFrame(yf2(i)-iFrame:yf2(i)+iFrame,xf2(i)-iFrame:xf2(i)+iFrame);
% 
%     crcoeff=(normxcorr2(box1,frame2));
%     mcrcoeff=max(crcoeff(:));
%     
%     [ssr,snd] = max(crcoeff);
%     [ij,ji] = ind2sub(size(crcoeff),snd);
% 
%     [ix,iy]=find(mcrcoeff==crcoeff);
%     ix=ix-crBox;
%     iy=iy-crBox;
% 
%     yt=[yt,iy(1)];
%     xt=[xt,ix(1)];
%---------------------------------------------------------------%


xx=[xf2(1)+iy1,xf2(2)+iy2,xf2(3)+iy3,xf2(4)+iy4,xf2(5)+iy5,xf2(6)+iy6,xf2(7)+iy7,xf2(8)+iy8, xf2(9)+iy9, xf2(10)+iy10, xf2(11)+iy11, xf2(12)+iy12];
yy=[yf2(1)+ix1,yf2(2)+ix2,yf2(3)+ix3,yf2(4)+ix4,yf2(5)+ix5,yf2(6)+ix6,yf2(7)+ix7,yf2(8)+ix8, yf2(9)+ix9, yf2(10)+ix10, yf2(11)+ix11, yf2(12)+ix12];  
    
% xx=xf;
% yy=yf;



%update xf2 and yf2




%   fprintf('size is %02d',size(Mmatx(:,6)));
%     fprintf('length is %02d',xx');

    Mmatx(:,v)=xx';
    Mmaty(:,v)=yy';




% ind=v-1;



% this step removes a sudden jump in tracking
% 
% difx=abs(Mmatx(:,v)-Mmatx(:,ind));
% dify=abs(Mmaty(:,v)-Mmaty(:,ind));
% 
% 
% ncrBox=round(crBox/2);
% 
% subi=find(difx>ncrBox);
% if isempty(subi)==0
%     for l=1:length(subi)
%         Mmatx(subi(l),v)=Mmatx(subi(l),ind);
%     end
% end
% 
% 
% subiy=find(dify>ncrBox);
% if isempty(subiy)==0
%     for l=1:length(subiy)
%         Mmaty(subiy(l),v)=Mmaty(subiy(l),ind);
%     end
% end
% xx=Mmatx(:,v)';
% yy=Mmaty(:,v)';

% for i=1:n
%     
%     Mmatx(n,v)=xx(1) ;
% 
%     mx=Mmatx(i,:);
%     if (mx(v)-mx(ind))>crBox
%         mx(v)=mx(ind);
%         MMatx(i,v)=mx(v);
%     end
% 
%         my=Mmaty(i,:);
%     if (my(v)-my(ind))>crBox
%         my(v)=my(ind);
%         Mmaty(i,v)=my(v);
%     end
%     
% end

mat1x(v)=xx(1);
mat2x(v)=xx(2);
mat3x(v)=xx(3);
mat4x(v)=xx(4);
mat5x(v)=xx(5);
mat6x(v)=xx(6);
mat7x(v)=xx(7);
mat8x(v)=xx(8);
mat9x(v)=xx(9);
mat10x(v)=xx(10);
mat11x(v)=xx(11);
mat12x(v)=xx(12);


mat1y(v)=yy(1);
mat2y(v)=yy(2);
mat3y(v)=yy(3);
mat4y(v)=yy(4);
mat5y(v)=yy(5);
mat6y(v)=yy(6);
mat7y(v)=yy(7);
mat8y(v)=yy(8);
mat9y(v)=yy(9);
mat10y(v)=yy(10);
mat11y(v)=yy(11);
mat12y(v)=yy(12);



imagesc(vidFrame), colormap('gray');
hold on
% plot(xf2(1)+iy1,yf2(1)+ix1,'gs-',xf2(2)+iy2,yf2(2)+ix2,'gs-',xf2(3)+iy3,yf2(3)+ix3,'gs-',xf2(4)+iy4,yf2(4)+ix4,'gs-',xf2(5)+iy5,yf2(5)+ix5,'gs-',xf2(6)+iy6,yf2(6)+ix6,'gs-',xf2(7)+iy7,yf2(7)+ix7,'gs-', ...
%     xf2(8)+iy8,yf2(8)+ix8,'gs-', xf2(9)+iy9,yf2(9)+ix9,'gs-',xf2(10)+iy10,yf2(10)+ix10,'gs-', xf2(11)+iy11,yf2(11)+ix11,'gs-',xf2(12)+iy12,yf2(12)+ix12,'gs-')

plot(xx(1),yy(1),'gs-',xx(2),yy(2),'gs-',xx(3),yy(3),'gs-',xx(4),yy(4),'gs-',xx(5),yy(5),'gs-',xx(6),yy(6),'gs-',xx(7),yy(7),'gs-', ...
    xx(8),yy(8),'gs-', xx(8),yy(8),'gs-',xx(9),yy(9),'gs-',xx(10),yy(10),'gs-',xx(11),yy(11),'gs-',xx(12),yy(12),'gs-')

plot(xx,yy,'r','LineWidth',2); %human

xxx=[xx xx(1)];
yyy=[yy yy(1)];
% fprintf('%2d \n',polyarea(xxx,yyy));
% plot(xx,yy,'r','LineWidth',1); %mice


hold off

F = getframe(gca);
writeVideo(vidfile, F);
    

% end

% figure;imagesc(box1)
% figure;imagesc(frame2)
% figure;imagesc(crcoeff);
% xtt(i,:)=xt-mean(xt);
% ytt(i,:)=yt-mean(yt);

end

close(vidfile)


% Calculations: filter and strain

%LPF of longitudinal displacement
% n=ceil(nF);
fc= 6;
imy1=mat1y(init:end);
my1=lpfFilter(imy1,fc);
% [~,~,imy1]=speckleFFT(nF, my1);


imy2=mat2y(init:end);
my2=lpfFilter(imy2,fc);

% [~,~,imy2]=speckleFFT(nF, my2);


imy3=mat3y(init:end);
my3=lpfFilter(imy3,fc);

% [~,~,imy3]=speckleFFT(nF, my3);


imy4=mat4y(init:end);
my4=lpfFilter(imy4,fc);

% [~,~,imy4]=speckleFFT(nF, my4);


imy5=mat5y(init:end);
my5=lpfFilter(imy5,fc);

% [~,~,imy5]=speckleFFT(nF, my5);


imy6=mat6y(init:end);
my6=lpfFilter(imy6,fc);

% [~,~,imy6]=speckleFFT(nF, my6);


imy7=mat7y(init:end);
my7=lpfFilter(imy7,fc);

% [~,~,imy7]=speckleFFT(nF, my7);


imy8=mat8y(init:end);
my8=lpfFilter(imy8,fc);

% [~,~,imy8]=speckleFFT(nF, my8);


imy9=mat9y(init:end);
my9=lpfFilter(imy9,fc);

% [~,~,imy9]=speckleFFT(nF, my9);


imy10=mat10y(init:end);
my10=lpfFilter(imy10,fc);

% [~,~,imy10]=speckleFFT(nF, my10);


imy11=mat11y(init:end);
my11=lpfFilter(imy11,fc);

% [~,~,imy11]=speckleFFT(nF, my11);


imy12=mat12y(init:end);
my12=lpfFilter(imy12,fc);

% [~,~,imy12]=speckleFFT(nF, my12);




%LPF of horizontal displacement
imx1=mat1x(init:end);
mx1=lpfFilter(imx1,fc);

% [~,~,imx1]=speckleFFT(nF, mx1);


imx2=mat2x(init:end);
mx2=lpfFilter(imx2,fc);

% [~,~,imx2]=speckleFFT(nF, mx2);

imx3=mat3x(init:end);
mx3=lpfFilter(imx3,fc);

% [~,~,imx3]=speckleFFT(nF, mx3);

imx4=mat4x(init:end);
mx4=lpfFilter(imx4,fc);

% [~,~,imx4]=speckleFFT(nF, mx4);

imx5=mat5x(init:end);
mx5=lpfFilter(imx5,fc);

% [~,~,imx5]=speckleFFT(nF, mx5);

imx6=mat6x(init:end);
mx6=lpfFilter(imx6,fc);

% [~,~,imx6]=speckleFFT(nF, mx6);

imx7=mat7x(init:end);
mx7=lpfFilter(imx7,fc);

% [~,~,imx7]=speckleFFT(nF, mx7);

imx8=mat8x(init:end);
mx8=lpfFilter(imx8,fc);

% [~,~,imx8]=speckleFFT(nF, mx8);

imx9=mat9x(init:end);
mx9=lpfFilter(imx9,fc);

% [~,~,imx9]=speckleFFT(nF, mx9);

imx10=mat10x(init:end);
mx10=lpfFilter(imx10,fc);

% [~,~,imx10]=speckleFFT(nF, mx10);


imx11=mat11x(init:end);
mx11=lpfFilter(imx11,fc);

% [~,~,imx11]=speckleFFT(nF, mx11);


imx12=mat12x(init:end);
mx12=lpfFilter(imx12,fc);

% [~,~,imx12]=speckleFFT(nF, mx12);

% ydel1=imy2-imy1;
% ydel2=imy3-imy2;
% ydel3=imy4-imy3;
% ydel4=imy5-imy4;
% ydel5=imy6-imy5;
% ydel6=imy7-imy6;
% ydel7=imy8-imy7;
% ydel8=imy9-imy8;
% ydel9=imy10-imy9;
% ydel10=imy11-imy10;
% ydel11=imy12-imy11;
% ydel12=imy1-imy12;

ydel1=my2-my1;
ydel2=my3-my2;
ydel3=my4-my3;
ydel4=my5-my4;
ydel5=my6-my5;
ydel6=my7-my6;
ydel7=my8-my7;
ydel8=my9-my8;
ydel9=my10-my9;
ydel10=my11-my10;
ydel11=my12-my11;
ydel12=my1-my12;


% xdel1=imx2-imx1;
% xdel2=imx3-imx2;
% xdel3=imx4-imx3;
% xdel4=imx5-imx4
% xdel5=imx6-imx5;
% xdel6=imx7-imx6;
% xdel7=imx8-imx7;
% xdel8=imx9-imx8;
% xdel9=imx10-imx9;
% xdel10=imx11-imx10;
% xdel11=imx12-imx11;
% xdel12=imy1-imy12;

xdel1=mx2-mx1;
xdel2=mx3-mx2;
xdel3=mx4-mx3;
xdel4=mx5-mx4;
xdel5=mx6-mx5;
xdel6=mx7-mx6;
xdel7=mx8-mx7;
xdel8=mx9-mx8;
xdel9=mx10-mx9;
xdel10=mx11-mx10;
xdel11=mx12-mx11;
xdel12=mx1-12;



%TIssue Doppler code
% figure;plot(my1),title('annulus motion');
% diss=sqrt(my1.^2+mx1.^2);
% figure;plot(-diff(my1)),title('annulus vel');
% figure;plot(diss); title('xy annulus motion');
% figure;plot(-diff(diss)),title('xy annnulus vel');
% figure;plot(-diff(my12)),title('annulus vel2');


%--full 2-3 cycles in the clip
dis1=sqrt(ydel1.^2+xdel1.^2);
bdis1=dis1-dis1(1);

normdis1=dis1./max(dis1);
% normdis1=bdis1/dis1(1);

% normdis1=dis1;
% normdis1(round(nF/2):end)=[];
str1=max(normdis1)-min(normdis1);
[str1,ndis1]=normdis(dis1);

dis2=sqrt(ydel2.^2+xdel2.^2);
bdis2=dis2-dis2(1);
normdis2=dis2./max(dis2);

% normdis2=dis2;
str2=max(normdis2)-min(normdis2);
[str2,ndis2]=normdis(dis2);


dis3=sqrt(ydel3.^2+xdel3.^2);
bdis3=dis3-dis3(1);
normdis3=dis3./max(dis3);

str3=max(normdis3)-min(normdis3);
[str3,ndis3]=normdis(dis3);



dis4=sqrt(ydel4.^2+xdel4.^2);
bdis4=dis4-dis4(1);
normdis4=dis4./max(dis4);
str4=max(normdis4)-min(normdis4);
[str4,ndis4]=normdis(dis4);

dis5=sqrt(ydel5.^2+xdel5.^2);
bdis5=dis5-dis5(1);
normdis5=dis5./max(dis5);
str5=max(normdis5)-min(normdis5);
[str5,ndis5]=normdis(dis5);


dis6=sqrt(ydel6.^2+xdel6.^2);
bdis6=dis6-dis6(1);
normdis6=dis6./max(dis6);
% normdis6(round(nF/2):end)=[];
str6=max(normdis6)-min(normdis6);
[str6,ndis6]=normdis(dis6);


dis7=sqrt(ydel7.^2+xdel7.^2);
bdis7=dis7-dis7(1);
normdis7=dis7./max(dis7);
str7=max(normdis7)-min(normdis7);
[str7,ndis7]=normdis(dis7);


dis8=sqrt(ydel8.^2+xdel8.^2);
bdis8=dis8-dis8(1);
normdis8=dis8./max(dis8);
str8=max(normdis8)-min(normdis8);
[str8,ndis8]=normdis(dis8);


dis9=sqrt(ydel9.^2+xdel9.^2);
bdis9=dis9-dis9(1);
normdis9=dis9./max(dis9);
str9=max(normdis9)-min(normdis9);
[str9,ndis9]=normdis(dis9);



dis10=sqrt(ydel10.^2+xdel10.^2);
bdis10=dis10-dis10(1);
normdis10=dis10./max(dis10);
str10=max(normdis10)-min(normdis10);
[str10,ndis10]=normdis(dis10);


dis11=sqrt(ydel11.^2+xdel11.^2);
bdis11=dis11-dis11(1);
normdis11=dis11./max(dis11);
str11=max(normdis11)-min(normdis11);
[str11,ndis11]=normdis(dis11);



dis12=sqrt(ydel12.^2+xdel12.^2);
bdis12=dis12-dis12(1);
normdis12=dis12./max(dis12);
str12=max(normdis12)-min(normdis12);
[str12,ndis12]=normdis(dis12);


mstr=(str1+str2+str3+str4+str5+str6+str7+str8+str9+str10+str11)/11;



% mat1y=mat1y-(max(mat1y));
% mat2y=mat2y-(max(mat2y));
% mat3y=mat3y-(max(mat3y));
% mat4y=mat4y-(max(mat4y));
% mat5y=mat5y-(max(mat5y));
% mat6y=mat6y-(max(mat6y));
% mat7y=mat7y-(max(mat7y));
% 
% mat1x=mat1x-(max(mat1x));
% mat2x=mat2x-(max(mat2x));
% mat3x=mat3x-(max(mat3x));
% mat4x=mat4x-(max(mat4x));
% mat5x=mat5x-(max(mat5x));
% mat6x=mat6x-(max(mat6x));
% mat7x=mat7x-(max(mat7x));



% need to input f0

% f0=5;
% NF=length(init:nF);
% [amp1,phase1,nds1] = speckleFFT (NF, normdis1,f0);
% [amp2,phase2,nds2] = speckleFFT (NF, normdis2,f0);
% [amp3,phase3,nds3] = speckleFFT (NF, normdis3,f0);
% [amp4,phase4,nds4] = speckleFFT (NF, normdis4,f0);
% [amp5,phase5,nds5] = speckleFFT (NF, normdis5,f0);
% [amp6,phase6,nds6] = speckleFFT (NF, normdis6,f0);
% [amp7,phase7,nds7] = speckleFFT (NF, normdis7,f0);
% [amp8,phase8,nds8] = speckleFFT (NF, normdis8,f0);
% [amp9,phase9,nds9] = speckleFFT (NF, normdis9,f0);
% [amp10,phase10,nds10] = speckleFFT (NF, normdis10,f0);
% [amp11,phase11,nds11] = speckleFFT (NF, normdis11,f0);
% [amp12,phase12,nds12] = speckleFFT (NF, normdis12,f0);
% 
% 
% amps=[amp1,amp2, amp3, amp4, amp5, amp6, amp7, amp8, amp9, amp10, amp11, amp12];
% phases=[phase1,phase2,phase3,phase4,phase5,phase6,phase7,phase8,phase9,phase10,phase11, phase12];
% 







x=init:nF;


h=figure;
xx=1:length(bdis1);
plot(xx,(bdis1+bdis2)/2,xx,(bdis3+bdis4)/2,xx,(bdis5+bdis6)/2,xx,(bdis7+bdis8)/2,xx,(bdis9+bdis10+bdis11)/3);
legend('seg1','seg2','seg3','seg4','seg5');


% figure;
% plot(xx,ndis1,xx,ndis2,xx,ndis3,xx,ndis4,xx,ndis5,xx,ndis6);

save([pathName fileName '2.mat'],'bdis1','bdis2','bdis4','bdis6','bdis10')

% ,xx,normdis6,xx,normdis7),...
%     xx,normdis8,xx,normdis9,xx,normdis10);

normdismat=[normdis1,normdis2,normdis3,normdis4,normdis5,normdis6,normdis7,normdis8,normdis9,normdis10,normdis11,normdis12];

[ax,~]=ginput(2);
% fprintf('%02d',ax);
normdismatN=normdismat(round(ax(1)):round(ax(2)),:);
cycleL=round(ax(2))-round(ax(1));
temp1 = diff(smoothdata(normdismatN,1));
[min_der, min_i]=min(smoothdata(temp1,1));


figure;
yyy=(bdis1+bdis2+bdis3+bdis4+bdis5+bdis6+bdis7+bdis8+bdis9+bdis10+bdis11)/11;
plot(xx,(bdis1+bdis2+bdis3+bdis4+bdis5+bdis6+bdis7+bdis8+bdis9+bdis10+bdis11)/11)
figure;plot(diff(yyy));
strainData=max(normdismatN)-min(normdismatN);

iii=strainData>0.6;
strainData(iii)=0;

indi=strainData~=0;
mStrain=mean(strainData(indi));
sprintf('mean strain is -%02f percent',mStrain*100)


% 
% % Find First Derivative and time of maxium
% temp2 = diff(temp,1,3); % first derivative
% [max_der max_i] = max(temp2,[],3); % find location of max derivative
% 
% 
% 
% figure;
% plot(x,nds1,x,nds2,x,nds3,x,nds4,x,nds5,x,nds6,x,nds7,...
% x, nds8, x, nds9, x, nds10, x, nds11,x,nds12);
% legend({'1','2','3','4','5','6','7','8','9','10','11','12'},'Location','southwest')
% 




% 
nvidFrame=uint8(vidFrame);
% nnvidFrame=zeros(size(nvidFrame));
% nnvidFrame(min(yf2(1),yf2(2)):max(yf2(1),yf2(2)),min(xf2(1),xf2(2)):max(xf2(1),xf2(2)))=phases(1);
% nnvidFrame(min(yf2(2),yf2(3)):max(yf2(2),yf2(3)),min(xf2(2),xf2(3)):max(xf2(2),xf2(3)))=phases(2);
% nnvidFrame(min(yf2(3),yf2(4)):max(yf2(3),yf2(4)),min(xf2(3),xf2(4)):max(xf2(3),xf2(4)))=phases(3);
% nnvidFrame(min(yf2(4),yf2(5)):max(yf2(4),yf2(5)),min(xf2(4),xf2(5)):max(xf2(4),xf2(5)))=phases(4);
% nnvidFrame(min(yf2(5),yf2(6)):max(yf2(5),yf2(6)),min(xf2(5),xf2(6)):max(xf2(5),xf2(6)))=phases(5);
% nnvidFrame(min(yf2(6),yf2(7)):max(yf2(6),yf2(7)),min(xf2(6),xf2(7)):max(xf2(6),xf2(7)))=phases(6);
% nnvidFrame(min(yf2(7),yf2(8)):max(yf2(7),yf2(8)),min(xf2(7),xf2(8)):max(xf2(7),xf2(8)))=phases(7);
% nnvidFrame(min(yf2(8),yf2(9)):max(yf2(8),yf2(9)),min(xf2(8),xf2(9)):max(xf2(8),xf2(9)))=phases(8);
% nnvidFrame(min(yf2(9),yf2(10)):max(yf2(9),yf2(10)),min(xf2(9),xf2(10)):max(xf2(9),xf2(10)))=phases(9);
% nnvidFrame(min(yf2(10),yf2(11)):max(yf2(10),yf2(11)),min(xf2(10),xf2(11)):max(xf2(10),xf2(11)))=phases(10);
% nnvidFrame(min(yf2(11),yf2(12)):max(yf2(11),yf2(12)),min(xf2(11),xf2(12)):max(xf2(11),xf2(12)))=phases(11);
% nnvidFrame(min(yf2(12),yf2(1)):max(yf2(12),yf2(1)),min(xf2(12),xf2(1)):max(xf2(12),xf2(1)))=phases(12);
% figure;imagesc(nnvidFrame),colormap(hsv),caxis([-pi pi]);



nnvidFrame2=zeros(size(nvidFrame));
nnvidFrame2(min(yf2(1),yf2(2)):max(yf2(1),yf2(2)),min(xf2(1),xf2(2)):max(xf2(1),xf2(2)))=strainData(1);
nnvidFrame2(min(yf2(2),yf2(3)):max(yf2(2),yf2(3)),min(xf2(2),xf2(3)):max(xf2(2),xf2(3)))=strainData(2);
nnvidFrame2(min(yf2(3),yf2(4)):max(yf2(3),yf2(4)),min(xf2(3),xf2(4)):max(xf2(3),xf2(4)))=strainData(3);
nnvidFrame2(min(yf2(4),yf2(5)):max(yf2(4),yf2(5)),min(xf2(4),xf2(5)):max(xf2(4),xf2(5)))=strainData(4);
nnvidFrame2(min(yf2(5),yf2(6)):max(yf2(5),yf2(6)),min(xf2(5),xf2(6)):max(xf2(5),xf2(6)))=strainData(5);
nnvidFrame2(min(yf2(6),yf2(7)):max(yf2(6),yf2(7)),min(xf2(6),xf2(7)):max(xf2(6),xf2(7)))=strainData(6);
nnvidFrame2(min(yf2(7),yf2(8)):max(yf2(7),yf2(8)),min(xf2(7),xf2(8)):max(xf2(7),xf2(8)))=strainData(7);
nnvidFrame2(min(yf2(8),yf2(9)):max(yf2(8),yf2(9)),min(xf2(8),xf2(9)):max(xf2(8),xf2(9)))=strainData(8);
nnvidFrame2(min(yf2(9),yf2(10)):max(yf2(9),yf2(10)),min(xf2(9),xf2(10)):max(xf2(9),xf2(10)))=strainData(9);
nnvidFrame2(min(yf2(10),yf2(11)):max(yf2(10),yf2(11)),min(xf2(10),xf2(11)):max(xf2(10),xf2(11)))=strainData(10);
nnvidFrame2(min(yf2(11),yf2(12)):max(yf2(11),yf2(12)),min(xf2(11),xf2(12)):max(xf2(11),xf2(12)))=strainData(11);
% nnvidFrame2(min(yf2(12),yf2(1)):max(yf2(12),yf2(1)),min(xf2(12),xf2(1)):max(xf2(12),xf2(1)))=amps(12);
% figure;imagesc(nnvidFrame2),colormap(jet), caxis([0 0.40])


nnvidFrame3=zeros(size(nvidFrame));
nnvidFrame3(min(yf2(1),yf2(2)):max(yf2(1),yf2(2)),min(xf2(1),xf2(2)):max(xf2(1),xf2(2)))=min_i(1);
nnvidFrame3(min(yf2(2),yf2(3)):max(yf2(2),yf2(3)),min(xf2(2),xf2(3)):max(xf2(2),xf2(3)))=min_i(2);
nnvidFrame3(min(yf2(3),yf2(4)):max(yf2(3),yf2(4)),min(xf2(3),xf2(4)):max(xf2(3),xf2(4)))=min_i(3);
nnvidFrame3(min(yf2(4),yf2(5)):max(yf2(4),yf2(5)),min(xf2(4),xf2(5)):max(xf2(4),xf2(5)))=min_i(4);
nnvidFrame3(min(yf2(5),yf2(6)):max(yf2(5),yf2(6)),min(xf2(5),xf2(6)):max(xf2(5),xf2(6)))=min_i(5);
nnvidFrame3(min(yf2(6),yf2(7)):max(yf2(6),yf2(7)),min(xf2(6),xf2(7)):max(xf2(6),xf2(7)))=min_i(6);
nnvidFrame3(min(yf2(7),yf2(8)):max(yf2(7),yf2(8)),min(xf2(7),xf2(8)):max(xf2(7),xf2(8)))=min_i(7);
nnvidFrame3(min(yf2(8),yf2(9)):max(yf2(8),yf2(9)),min(xf2(8),xf2(9)):max(xf2(8),xf2(9)))=min_i(8);
nnvidFrame3(min(yf2(9),yf2(10)):max(yf2(9),yf2(10)),min(xf2(9),xf2(10)):max(xf2(9),xf2(10)))=min_i(9);
nnvidFrame3(min(yf2(10),yf2(11)):max(yf2(10),yf2(11)),min(xf2(10),xf2(11)):max(xf2(10),xf2(11)))=min_i(10);
nnvidFrame3(min(yf2(11),yf2(12)):max(yf2(11),yf2(12)),min(xf2(11),xf2(12)):max(xf2(11),xf2(12)))=min_i(11);
% nnvidFrame3(min(yf2(12),yf2(1)):max(yf2(12),yf2(1)),min(xf2(12),xf2(1)):max(xf2(12),xf2(1)))=min_i(12);
ind4=nnvidFrame3~=0;
mTime=mean(nnvidFrame3(ind4));
sprintf('mean activation time is %02f',mTime)




axl=round(ax(end)-ax(1));
% figure;imagesc(nnvidFrame3),colormap(mymap), caxis([0 axl]);
save([pathName fileName '1.mat'],'nnvidFrame3')

% 
% nnvidFrame(yf2(2),xf2(2))=phases(2);
% nnvidFrame(yf2(3),xf2(3))=phases(3);
% nnvidFrame(yf2(4),xf2(4))=phases(4);
% nnvidFrame(yf2(5),xf2(5))=phases(5);
% nnvidFrame(yf2(6),xf2(6))=phases(6);
% nnvidFrame(yf2(7),xf2(7))=phases(7);
% nnvidFrame(yf2(8),xf2(8))=phases(8);
% nnvidFrame(yf2(9),xf2(9))=phases(9);
% nnvidFrame(yf2(10),xf2(10))=phases(10);
% nnvidFrame(yf2(11),xf2(11))=phases(11);
% mean(phases);


% figure
% nvidFrame=uint8(vidFrame);
% imagesc(nvidFrame)
% % , colormap('gray');
% hold on
% plot(xf2(1)+iy1,yf2(1)+ix1,'gs-',xf2(2)+iy2,yf2(2)+ix2,'gs-',xf2(3)+iy3,yf2(3)+ix3,'gs-',xf2(4)+iy4,yf2(4)+ix4,'gs-',xf2(5)+iy5,yf2(5)+ix5,'gs-',xf2(6)+iy6,yf2(6)+ix6,'gs-',xf2(7)+iy7,yf2(7)+ix7,'gs-')
% plot(xx,yy,'r');



% figure;
% plot(x,imy1,x,imy2,x,imy3,x,imy4,x,imy5,x,imy6,x,imy7),title('y_filtered')
% figure;
% plot(x,imy1,x,imy2,x,imy3,x,imy4,x,imy5,x,imy6,x,imy7),title('y_original')


% plot(x2,str1pc,x2,str2pc,x2,str3pc,x2,str4pc,x2,str5pc,x2,str6pc,x2,str7pc);

% dely1=smoothdata(mat1y(2:end))-smoothdata(mat2y(2:end));
% dely2=smoothdata(mat2y(2:end))-smoothdata(mat3y(2:end));
% dely3=smoothdata(mat3y(2:end))-smoothdata(mat4y(2:end));
% dely4=smoothdata(mat4y(2:end))-smoothdata(mat5y(2:end));
% dely5=smoothdata(mat5y(2:end))-smoothdata(mat6y(2:end));
% dely6=smoothdata(mat6y(2:end))-smoothdata(mat7y(2:end));


% odely1=abs(xf2(1)-xf2(2));
% odely2=abs(xf2(2)-xf2(3));
% odely3=abs(xf2(3)-xf2(4));
% odely4=abs(xf2(4)-xf2(5));
% odely5=abs(xf2(5)-xf2(6));
% odely6=abs(xf2(6)-xf2(7));
% 
% 
% figure;
% plot(x,dely1,x,dely2,x,dely3,x,dely4,x,dely5,x,dely6);
% 
% 
% 
% figure;
% plot(x,smoothdata(mat1x(2:end)),x,smoothdata(mat2x(2:end)),x,smoothdata(mat3x(2:end)),x,smoothdata(mat4x(2:end)),x,smoothdata(mat5x(2:end)),x,smoothdata(mat6x(2:end)),x,smoothdata(mat7x(2:end)));
% title('x');
% hold on
% plot(medfilt1(xtt(2,:)-xtt(1,:),3))
% plot(medfilt1(xtt(3,:)-xtt(2,:)))
% plot(medfilt1(xtt(4,:)-xtt(3,:)))
% plot(medfilt1(xtt(5,:)-xtt(3,:)))
% plot(medfilt1(xtt(6,:)-xtt(5,:)))
% plot(medfilt1(xtt(7,:)-xtt(6,:)))
% plot(medfilt1(xtt(8,:)-xtt(7,:)))
% plot(medfilt1(xtt(9,:)-xtt(8,:)))
% plot(medfilt1(xtt(10,:)-xtt(9,:),3))
% 
% for w=1:20
%     plot(xtt(:,w),ytt(:,w));
% end
close(f);
close(g);
close(h);