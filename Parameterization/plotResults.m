function M=plotResults(numPoints,dataForPlotting,m,r,A5_0)


theta4plot=[0:0.1:2*pi];

figure
colors = get(gca,'colororder');
view([-65 10]);
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'xlim',[0 3],'ylim',[-1.5 1.5],'zlim',[0 3])
axis equal


hold on
for i=1:numPoints
    lineH(i)=line([0,0],[0,0],[0,0],'color',colors(i,:),'marker','.');
    circleH(i)=line([0,0],[0,0],[0,0],'color',colors(i,:));
    centerH(i)=line(0,0,0,'color',colors(i,:),'marker','.');
    gammaRMat(:,i)=m.( ['gammaR' num2str(i)] );
end
for soln=1:size(dataForPlotting,1)
    centerPts = equationSet(dataForPlotting(soln,:));
    for ptn=1:numPoints
        
        %Update the circle x/y/z data
        %circleParam is not vectorized, so have to loop
        for j=1:length(theta4plot)
            cData(:,j)=circleParam(theta4plot(j),r,centerPts(1,ptn),centerPts(2,ptn),centerPts(3,ptn),m.phiA(soln),m.psiA(soln));      
        end      
        set(circleH(ptn),'xdata',cData(1,:),'ydata',cData(2,:),'zdata',cData(3,:))
        
       %Update the line x/y/z data
       l=lineParam(A5_0{ptn},gammaRMat(soln,ptn));
       pLaser=A5_0{ptn};
       pLaser=pLaser(1:3,4);
       set(lineH(ptn),'xdata',[pLaser(1) l(1)],'ydata',[pLaser(2) l(2)],'zdata',[pLaser(3) l(3)]);
       
       %Update Center Position
       set(centerH(ptn),'xdata',[centerPts(1,ptn)],'ydata',[centerPts(2,ptn)],'zdata',[centerPts(3,ptn)]);
        
    end
    drawnow
    M(soln)=getframe;
    
end
