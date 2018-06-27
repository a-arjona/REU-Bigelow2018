%% Build Average Time Series by Cluster
clc
clear all
close all

mean_ts_clu = zeros(120, 28);
years = [1998:2017];
cpt = 0;

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
    col = cbrewer('qual','Set1',max(cluster(:,1)));
    
    for ccc = 1:6
        f = find(cluster(:,1) == ccc); % index of all with cluster number ccc
        cpt = cpt + 1;
        mean_ts_clu(cpt,1) = ccc;
        mean_ts_clu(cpt,2:28) = mean(timeseries(f,:)*mean(maxo(f))); % avg time series of cluster for that year
        
    end
            
end

mean_ts_by_clu = zeros(6,27);
sd = zeros(6,27);

for clu = 1:6
    f = find(mean_ts_clu(:,1) == clu);
    mean_ts_by_clu(clu,:) = mean(mean_ts_clu(f,2:28)*mean(maxo(f)));
    sd(clu,:) = std(mean_ts_clu(f,2:28)*mean(maxo(f)));
end
    
M = [3:10];
D = ones(1,length(M));
Y = repmat(2000,1,length(M));
datem = datenum('01-01-2000')+ [3:8:365];
xtick = datenum(Y,M,D);
figure()
for clu = 1:6
    
    subplot(2,3,clu)
    plot(datem(9:35), mean_ts_by_clu(clu,:),'-','color',col(clu,:),'LineWidth',2)
    hold on
    plot(datem(9:35), mean_ts_by_clu(clu,:)+sd(clu,:),'k-')
    plot(datem(9:35), mean_ts_by_clu(clu,:)-sd(clu,:),'k-')
    set(gca,'Ylim',[0 5],'Xtick',xtick,'Xticklabel',datestr(xtick,'mmm'),'Xlim',[xtick(1) xtick(end)+10],...
        'Xgrid','on','Ygrid','on','fontsize',9)
    ylabel('[Chl-a] (mg m^-^3)','fontsize',9)
    title("Average Timeseries by Cluster over 20 years")
end 

%% Average Si

lat = 60.125  : .25 : 80.875;
lon = -43.875 : .25 : 17.875;

[LON, LAT] = meshgrid(lon,lat);
years = [1998:2017];

% open mean si data
    T = readtable('C:/Files/Work/Bigelow/Data/txt_files/mean_si.txt');
    
    m_si     = table2array(T(:,1));
    xo          = table2array(T(:,2));
    yo          = table2array(T(:,3));

figure(1)
% plot map
    map_ccc = nan(size(LON));
    
    for ii = 1:length(m_si)
        map_ccc(yo(ii),xo(ii)) = m_si(ii,1);
    
    end
    set(0,'DefaultFigureRenderer','zbuffer')
    m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]);
    colormap(jet)
    
    [X,Y] = m_ll2xy(LON,LAT);
    h = pcolor(X,Y,map_ccc);
    set(h,'edgecolor','none')
    
    tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
    set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])
    
    m_grid()
    
    %pic_name = strcat('C:/Users/Ade/Documents/Bigelow/DINEOF-2018/Interannual-Data/m_figure/MODEmap.png');
    %saveas(gcf,pic_name)

    %change so warm colors above mean(m_si)
%% Percent time associated with mode cluster
lat = 60.125  : .25 : 80.875;
lon = -43.875 : .25 : 17.875;

[LON, LAT] = meshgrid(lon,lat);
years = [1998:2017];

% open mean si data
    T = readtable('C:/Files/Work/Bigelow/Data/txt_files/perc_assoc.txt');
    
    perc_assoc     = table2array(T(:,1));
    xo          = table2array(T(:,2));
    yo          = table2array(T(:,3));
    
pa = perc_assoc .* 100;

% plot map
    map_ccc = nan(size(LON));
    
    for ii = 1:length(pa)
        map_ccc(yo(ii),xo(ii)) = pa(ii,1);
    
    end
    set(0,'DefaultFigureRenderer','zbuffer')
    m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]);
    colormap(jet)
    
    [X,Y] = m_ll2xy(LON,LAT);
    h = pcolor(X,Y,map_ccc);
    set(h,'edgecolor','none')
    
    tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
    set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])
    
    m_grid()
    title("Percent Associated with Mode Cluster")

