%test tanh shape
MeshHeight = 1.0;
MeshWidth  = 1.0;
LABHeight  = 0.9;

x = [0:0.01:1];

disc = (0.2 * MeshHeight).^2 - (x - 0.5 * MeshWidth).^2;
%hump = sqrt(disc);

shape = LABHeight - hump;

figure(1);clf
%plot(x,shape,'o'); set(gca,'xlim',[0 1],'ylim',[0 1]); hold on

height = 0.1;
width = 0.2;
scale = 0.05;
ridge1 = height*(1-tanh((x - (0.5 + width)* MeshWidth)/scale));
ridge2 = height*(1-tanh((x - (0.5 - width)* MeshWidth)/scale));
ridge = ridge1 - ridge2;

shape = LABHeight - ridge1 + ridge2;

plot(x,shape,'o'); set(gca,'xlim',[0 1],'ylim',[0 1]);

figure
plot(x,ridge1,'ro'); set(gca,'xlim',[0 1],'ylim',[0 1]);hold on
plot(x,ridge2,'o'); set(gca,'xlim',[0 1],'ylim',[0 1]);

figure
xnew     = x - MeshWidth*.5;
ridgenew = height*( - tanh((x-0.5-width)/scale) + tanh((x-0.5+width)/scale));
depth = 0.1 + ridge;
depthnew = 0.1 + ridgenew;
plot(xnew,depth,'o'); set(gca,'xlim',[-0.5 0.5],'ylim',[0 1]);hold on
plot(xnew,depthnew,'k--'); set(gca,'xlim',[-0.5 0.5],'ylim',[0 1]);


