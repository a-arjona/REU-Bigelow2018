%% Plot All Cluster (Annual)
clear all
close all
clc

lat = 60.125  : .25 : 80.875;
lon = -43.875 : .25 : 17.875;

[LON, LAT] = meshgrid(lon,lat);
years = [1998:2017];

for year = 1:20
    
    subplot(4,5,year)
    
    c_year = num2str(years(year));
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'metadata.txt');
    % open metadata
    T = readtable(filename);
    
    maxo        = table2array(T(:,1));
    xo          = table2array(T(:,2));
    yo          = table2array(T(:,3));
    good_points = table2array(T(:,4));
    
    % open clusters
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'clusters.txt');
    T = readtable(filename, 'Delimiter','space', 'ReadVariableNames',false);
    
    cluster     = table2array(T);
    
    % plot map
    map_ccc = nan(size(LON));
    
    for ii = 1:length(xo)
        map_ccc(yo(ii),xo(ii)) = cluster(ii,1);
    
    end
    set(0,'DefaultFigureRenderer','zbuffer')
    m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]);
    col = cbrewer('qual','Set1',max(cluster(:,1)));
     
    colormap(col)
    [X,Y] = m_ll2xy(LON,LAT);
    h = pcolor(X,Y,map_ccc);
    set(h,'edgecolor','none')
    
    tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
    set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])
    
    m_grid()
    m_gshhs_l('patch',[0.8 0.8 0.8],'edgecolor','k');
    title(c_year)
    
    
            
end

%pic_name = strcat('C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/m_figure/test_all_clusters.png');
%saveas(gcf,pic_name)
%% Time Series
for year = 1:20
    
    c_year = num2str(years(year));
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'metadata.txt');
    % open metadata
    T = readtable(filename);
    
    maxo        = table2array(T(:,1));
    xo          = table2array(T(:,2));
    yo          = table2array(T(:,3));
    good_points = table2array(T(:,4));
    
    % open clusters
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'clusters.txt');
    T = readtable(filename, 'Delimiter','space', 'ReadVariableNames',false);
    
    cluster     = table2array(T);
    
    % open timeseries
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'timeseries.txt');
    T = readtable(filename);
    
    timeseries  = table2array(T);
    
    % plot time series
    figure(year + 1)
    
    M = [3:10];
    D = ones(1,length(M));
    Y = repmat(2000,1,length(M));
    datem = datenum('01-01-2000')+ [3:8:365];
    xtick = datenum(Y,M,D);
    
    cpt = 0;
    sub = [1:6];
    for ccc = [1 2 3 4 5 6]
        f = find(cluster(:,1) == ccc);
        cpt = cpt + 1;
        subplot(3,3,sub(cpt))
        plot(datem(9:35), mean(timeseries(f,:))*mean(maxo(f)),'-','color',col(ccc,:),'LineWidth',2)
        hold on
        plot(datem(9:35), mean(timeseries(f,:))*mean(maxo(f))+std(timeseries(f,:))*mean(maxo(f)),'k-')
        plot(datem(9:35), mean(timeseries(f,:))*mean(maxo(f))-std(timeseries(f,:))*mean(maxo(f)),'k-')
        set(gca,'Ylim',[0 3],'Xtick',xtick,'Xticklabel',datestr(xtick,'mmm'),'Xlim',[xtick(1) xtick(end)+10],...
            'Xgrid','on','Ygrid','on','fontsize',9)
        ylabel('[Chl-a] (mg m^-^3)','fontsize',9)
        %title(c_year)
    end
   pic_name = strcat('C:/Files/Work/Bigelow/Data/figures/Annual/',c_year,'ts.png');
   saveas(gcf,pic_name)
end
%%  Silhouette Analysis
clear all
close all
clc
cpt = 0;

for year = 1:20 
    
    years = [1998:2017];
    c_year = num2str(years(year));
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'metadata.txt');
    % open metadata
    T = readtable(filename);
    
    maxo        = table2array(T(:,1));
    xo          = table2array(T(:,2));
    yo          = table2array(T(:,3));
    good_points = table2array(T(:,4));
    
    % open clusters
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'clusters.txt');
    T = readtable(filename, 'Delimiter','space', 'ReadVariableNames',false);
    
    cluster     = table2array(T);
    
    % open si
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'si.txt');
    T = readtable(filename, 'Delimiter','space', 'ReadVariableNames',false);
    si  = table2array(T);
    
    % plot silhouette
    col = cbrewer('qual','Set1',9);
    
    sub = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];

    
    x = 0;
    cpt = cpt + 1;
    nclu = 1;
    for ccc = 1:6
        subplot(4,5,sub(cpt))
        
        f = find(cluster(:,nclu) == ccc);
        Y = sort(si(f,nclu),'descend');
        X = [x:(x+length(Y)-1)];
        area(X,Y,'FaceColor',col(ccc,:),'LineStyle','none')
        hold on
        x = X(end)+1;
        set(gca,'Ylim',[-0.2 0.55],'Xlim',[-100 length(cluster)+100], 'fontsize',9)
        text(prctile(X,25),prctile(si(:,nclu),25),[num2str(round(sum(Y < mean(si(:,nclu))).*100./length(Y))),'%'], 'fontsize',9)
        text(prctile(X,25),prctile(si(:,nclu),75),[num2str(round(sum(Y >= mean(si(:,nclu))).*100./length(Y))),'%'], 'fontsize',9)
    
    end
    plot([0 length(cluster)],[mean(si(:,1)) mean(si(:,1))],'k-')
    ylabel('Silhouette value', 'fontsize',9)
    title(c_year)
end
%pic_name = strcat('C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/m_figure/all_silo.png');
%saveas(gcf,pic_name)
%% Plot All Maps (Individually)
lat = 60.125  : .25 : 80.875;
lon = -43.875 : .25 : 17.875;

[LON, LAT] = meshgrid(lon,lat);
years = [1998:2017];

for year = 1:20
    
    figure(year+30)
    figure()
    c_year = num2str(years(year));
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'metadata.txt');
    % open metadata
    T = readtable(filename);
    
    maxo        = table2array(T(:,1));
    xo          = table2array(T(:,2));
    yo          = table2array(T(:,3));
    good_points = table2array(T(:,4));
    
    % open clusters
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'clusters.txt');
    T = readtable(filename, 'Delimiter','space', 'ReadVariableNames',false);
    
    cluster     = table2array(T);
    
    % plot map
    map_ccc = nan(size(LON));
    
    for ii = 1:length(xo)
        map_ccc(yo(ii),xo(ii)) = cluster(ii,1);
    
    end
    set(0,'DefaultFigureRenderer','zbuffer')
    m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]);
    col = cbrewer('qual','Set1',max(cluster(:,1)));
    colormap(col)
    [X,Y] = m_ll2xy(LON,LAT);
    h = pcolor(X,Y,map_ccc);
    set(h,'edgecolor','none')
    
    tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
    set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])
    
    m_grid()
    m_gshhs_l('patch',[0.8 0.8 0.8],'edgecolor','k');
    title(c_year)
   
    pic_name = strcat('C:/Files/Work/Bigelow/Data/figures/Annual/',c_year,'map.png');
    saveas(gcf,pic_name)
end

%% Plot Silhouette Analysis (Individually)

clear all
close all
clc
cpt = 0;

for year = 1:20 
    
    years = [1998:2017];
    c_year = num2str(years(year));
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'metadata.txt');
    % open metadata
    T = readtable(filename);
    
    maxo        = table2array(T(:,1));
    xo          = table2array(T(:,2));
    yo          = table2array(T(:,3));
    good_points = table2array(T(:,4));
    
    % open clusters
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'clusters.txt');
    T = readtable(filename, 'Delimiter','space', 'ReadVariableNames',false);
    
    cluster     = table2array(T);
    
    % open si
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'si.txt');
    T = readtable(filename, 'Delimiter','space', 'ReadVariableNames',false);
    si  = table2array(T);
    
    % plot silhouette
    col = cbrewer('qual','Set1',9);
    
    sub = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];

    
    x = 0;
    cpt = cpt + 1;
    nclu = 1;
    
    figure(sub(cpt))
    for ccc = 1:6
        
        
        f = find(cluster(:,nclu) == ccc);
        Y = sort(si(f,nclu),'descend');
        X = [x:(x+length(Y)-1)];
        area(X,Y,'FaceColor',col(ccc,:),'LineStyle','none')
        hold on
        x = X(end)+1;
        set(gca,'Ylim',[-0.2 0.55],'Xlim',[-100 length(cluster)+100], 'fontsize',9)
        text(prctile(X,25),prctile(si(:,nclu),25),[num2str(round(sum(Y < mean(si(:,nclu))).*100./length(Y))),'%'], 'fontsize',9)
        text(prctile(X,25),prctile(si(:,nclu),75),[num2str(round(sum(Y >= mean(si(:,nclu))).*100./length(Y))),'%'], 'fontsize',9)
    
    end
    plot([0 length(cluster)],[mean(si(:,1)) mean(si(:,1))],'k-')
    ylabel('Silhouette value', 'fontsize',9)
    title(c_year)
    pic_name = strcat('C:/Files/Work/Bigelow/Data/figures/Annual/',c_year,'silho.png');
    saveas(gcf,pic_name)
end
    
%% Plot Mode Clusters Map

lat = 60.125  : .25 : 80.875;
lon = -43.875 : .25 : 17.875;

[LON, LAT] = meshgrid(lon,lat);
years = [1998:2017];

% open clusters
    T = readtable('C:/Files/Work/Bigelow/Data/txt_files/m_clusters.txt');
    
    cluster     = table2array(T(:,1));
    xo          = table2array(T(:,2));
    yo          = table2array(T(:,3));

figure(1)
% plot map
    map_ccc = nan(size(LON));
    
    for ii = 1:length(cluster)
        map_ccc(yo(ii),xo(ii)) = cluster(ii,1);
    
    end
    set(0,'DefaultFigureRenderer','zbuffer')
    m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]);
    col = cbrewer('qual','Set1',max(cluster(:,1)));
    colormap(col)
    [X,Y] = m_ll2xy(LON,LAT);
    h = pcolor(X,Y,map_ccc);
    set(h,'edgecolor','none')
    
    tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
    set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])
    
    m_grid()
    m_gshhs_l('patch',[0.8 0.8 0.8],'edgecolor','k');
    
   
    pic_name = strcat('C:/Files/Work/Bigelow/Data/figures/Summary/MODEmap.png');
    saveas(gcf,pic_name)
