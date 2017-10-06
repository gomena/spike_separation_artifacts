function PlotArtifact(Art,cond,center,patterns,nh,nv,numfig,sizes,axiss,factor,colors)
%Plots artifact (condition cond, centered at center, with stimulating
%electrodes paterns (on a different scale, on a nhxnv portion of the array,
%on figure numfig
%Gonzalo Mena May/2016

hy0=.1;
hx0=0;
%% Make ducks in the pond
%defaults nh=18,nv=14

hf=figure(numfig);

set(hf,'units','normalized');
set(hf,'position',sizes)
set(hf,'color',[1 1 1])
[xc,yc] = getElectrodeCoords512();
[a b]=sort(abs(unique(xc)-xc(center)),'ascend');
[c d]=sort(abs(unique(yc)-yc(center)),'ascend');
pat1=find(abs(xc-xc(center))<=a(nh));
pat2=find(abs(yc-yc(center))<=c(nv));
pat=intersect(pat1,pat2);
xc=xc(pat);
yc=yc(pat);
xc=(xc-min(xc));
xc=xc/max(xc)*(1-2*hx0)+hx0;
yc=(yc-min(yc));
yc=yc/max(yc)*(1-2*hy0)+hy0;
hx=0.075;
hy=0.8;

h=[hx hy];

axis 'off'
%set electrodes for plotting

for i=1:length(pat)


if(~isempty(find(pat(i)==patterns)))
    
pos2(1:2)=[xc(i) yc(i)];
pos2(2)=pos2(2)-0.20;
pos2(3:4)=[h(1) h(2)*factor(1)];
   

h3=axes('position',pos2);
axis 'off'
hold on

    plot(squeeze(Art(cond,pat(i),1:20))','color',colors{1},'linewidth',2)
axis(axiss{1})
%axis([1 20 -100 200])
axis 'off'
hold on
else
    
pos2(1:2)=[xc(i) yc(i)];
pos2(2)=pos2(2)-0.25;
pos2(3:4)=h;
h3=axes('position',pos2);
axis 'off'
hold on
plot(squeeze(Art(cond,pat(i),1:20))','color',colors{2},'linewidth',2)
axis(axiss{2})
end
axis off
hold on
end

