function plotClusterInfo(clusterInfo)
    figure();
    clusterSizeInfo = clusterInfo(:,1);
    clusterConductInfo = clusterInfo(:,2);
    subplot(3,1,1);
    hist(clusterSizeInfo);
    title('Histogram of Cluster Sizes');
    subplot(3,1,2);
    hist(clusterConductInfo);
    title('Histogram of Cluster Conductances');
    subplot(3,1,3);
    
    %//amap = [.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;];
    %alpha(.5);
    plot(clusterSizeInfo, clusterConductInfo, 'x');%, 'FaceAlpha', .);%, 'MarkerFaceColor', 'blue');%, 'Alphamap', alphamap);
    
    %//alphamap(amap);
    title('Conductance vs Cluster Size Scatter Plot');
    hold off;