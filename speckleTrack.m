    
function [ix, iy]= speckleTrack (vidFrame, yIndex,xIndex,iFrame,box1,crBox,iBox,zg)


frame2=vidFrame(yIndex-iFrame:yIndex+iFrame,xIndex-iFrame:xIndex+iFrame);

%     frame2=vidFrame2(ys(i)-50:ys(i)+50,xs(i)-50:xs(i)+0);

%    adaptive box 1
%     box1 = vidFrame(yIndex-iBox:yIndex+iBox,xIndex-iBox:xIndex+iBox);


% modify
% xx=normxcorr2(box1,frame2);
% save(['temp2.mat'],'xx')
% 
% adfadf
    

% iterative method _ optimization
%     a=normxcorr2(box1,frame2).*zg;
%     crcoeff=iterate(a,iBox,iFrame);
%     

% average method: medfilt2
    crcoeff=medfilt2(normxcorr2(box1,frame2),[5 5]).*zg;
    
    
%    
%     [ssr,snd] = max(crcoeff);
% %     [ij,ji] = ind2sub(size(crcoeff),snd);

    mcrcoeff=max(crcoeff(:));
    [ix,iy]=find(mcrcoeff==crcoeff);
    
%     fprintf('ix is %02d',ix);
    
%     ix=ix-crBox;
%     iy=iy-crBox;
    ix=ix(1)-crBox;
    iy=iy(1)-crBox;
    


%     yt=[yt,iy(1)];
%     xt=[xt,ix(1)];