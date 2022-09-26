
% c is corrcoeff;

function [out]= iterate (c,iBox,iFrame)


threshold=0.7;

for i=1:4
d=mean(c,2); %ydim
e=mean(c,1); %xdim


% step1
maxd=max(d);
hmaxd=maxd*threshold;
% hhmaxd=maxd*0.7;
ind=find(d>hmaxd);
mindy=round(mean(ind));

%step2
maxe=max(e);
hmaxe=maxe*threshold;
% hhmaxe=maxe*0.7;
inde=find(e>hmaxe);
mindx=round(mean(inde));


% step3: make PDF at mindx and mindy

% iBox=35;
% iFrame=60;

n=(iBox+iFrame)*2+1;
x=1:n;


yx=exp(-(x-mindx).^2/(iFrame)); %mice
yy=exp(-(x-mindy).^2/(iFrame)); 

z= ones(n,n)*maxe;
y1=repmat(yx,n,1);
y2=repmat(yy',1,n);
zg=z.*y1.*y2;
c=c+zg;
% figure;imagesc(c);

end

out=c;
% 
% 
% figure;imagesc(zg);
% figure;imagesc(a.xx);
% figure;imagesc(a.xx+zg);


