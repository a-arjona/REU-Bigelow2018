%% Plot results from Fuzzy Cluster analysis in Sub-Arctic Atlantic
%
% also plot ts for each cluster (centers data)
%

clear all
close all
clc

% open fuzzy clustering data
T = readtable('C:/Files/Work/Bigelow/Data/txt_files/19982007fuzzy_cl.txt');

xo          = table2array(T(:,7));
yo          = table2array(T(:,8));

nclu = 6;
nclu = nclu - 1;
len = 6;

white = [1, 1, 1];
figure()
col = cbrewer('qual','Set1',6);
for i = 1:6 
    k = subplot(2,3,i);
    colors_p = [linspace(white(1),col(i,1),len)',linspace(white(2),col(i,2),len)',linspace(white(3),col(i,3),len)'];
    clu = table2array(T(:,i));
    
    lat = 60.125  : .25 : 80.875;
    lon = -43.875 : .25 : 17.875;

    [LON, LAT] = meshgrid(lon,lat);

    map_ccc = nan(size(LON));
    for ii = 1:length(xo)
        map_ccc(yo(ii),xo(ii)) = clu(ii,1);

    end
    
    set(0,'DefaultFigureRenderer','zbuffer') % une ligne de code necessaire car sinon parfois ca bug au moment de faire la figure (#JeSaisPasPourquoi)
    m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]); % Projection a appliquer et bien specififier les latitudes limites

    
    %ncol = [col(1,:); col(2,:); col(4,:); col(3,:); col(5,:); col(6,:)];
    colormap(k, colors_p)
    [X,Y] = m_ll2xy(LON,LAT);
    h = pcolor(X,Y,map_ccc);
    set(h,'edgecolor','none')

    tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
    set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])

    m_grid()
    m_gshhs_l('patch',[0.8 0.8 0.8],'edgecolor','k');
    j = num2str(i);
    gt = strcat('Cluster ', j)
    title(gt);
    
    
end

%     
% figure(2)
% lat = 60.125  : .25 : 80.875;
% lon = -43.875 : .25 : 17.875;
% 
% [LON, LAT] = meshgrid(lon,lat);
% 
% map_ccc = nan(size(LON));
% for ii = 1:length(xo)
%     map_ccc(yo(ii),xo(ii)) = cluster(ii,1);
% 
% end
% 
% set(0,'DefaultFigureRenderer','zbuffer') % une ligne de code necessaire car sinon parfois ca bug au moment de faire la figure (#JeSaisPasPourquoi)
% m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]); % Projection a appliquer et bien specififier les latitudes limites
% 
%     
% %ncol = [col(1,:); col(2,:); col(4,:); col(3,:); col(5,:); col(6,:)];
% colormap(col)
% [X,Y] = m_ll2xy(LON,LAT);
% h = pcolor(X,Y,map_ccc);
% set(h,'edgecolor','none')
% 
% tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
% set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])
% 
% m_grid()
% m_gshhs_l('patch',[0.8 0.8 0.8],'edgecolor','k');
% title('20-year Mean Cluster (by FCM)')

%% FCM plots for annual data

clear all
close all
clc

years = [1998:2017];

for year = 1:20
    
    c_year = num2str(years(year));
    filename = strcat('C:/Files/Work/Bigelow/Data/txt_files/', c_year, 'FCM.txt');
    
    % open fuzzy clustering data
    T = readtable(filename);

    xo          = table2array(T(:,7));
    yo          = table2array(T(:,8));

    cluster     = table2array(T(:,9));


    nclu = 6;
    nclu = nclu - 1;
    len = 6;

    white = [1, 1, 1];
    figure(year)
    col = cbrewer('qual','Set1',max(cluster(:,1)));
    for i = 1:6 
        k = subplot(2,3,i);
        colors_p = [linspace(white(1),col(i,1),len)',linspace(white(2),col(i,2),len)',linspace(white(3),col(i,3),len)'];
        clu = table2array(T(:,i));

        lat = 60.125  : .25 : 80.875;
        lon = -43.875 : .25 : 17.875;

        [LON, LAT] = meshgrid(lon,lat);

        map_ccc = nan(size(LON));
        for ii = 1:length(xo)
            map_ccc(yo(ii),xo(ii)) = clu(ii,1);

        end

        set(0,'DefaultFigureRenderer','zbuffer') % une ligne de code necessaire car sinon parfois ca bug au moment de faire la figure (#JeSaisPasPourquoi)
        m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]); % Projection a appliquer et bien specififier les latitudes limites


        %ncol = [col(1,:); col(2,:); col(4,:); col(3,:); col(5,:); col(6,:)];
        colormap(k, colors_p)
        [X,Y] = m_ll2xy(LON,LAT);
        h = pcolor(X,Y,map_ccc);
        set(h,'edgecolor','none')

        tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
        set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])

        m_grid()
        m_gshhs_l('patch',[0.8 0.8 0.8],'edgecolor','k');
        j = num2str(i);
        gt = strcat('Cluster ', j);
        title(gt);


    end   
 
    % hard clusters map
    figure(20+year)
    lat = 60.125  : .25 : 80.875;
    lon = -43.875 : .25 : 17.875;

    [LON, LAT] = meshgrid(lon,lat);

    map_ccc = nan(size(LON));
    for ii = 1:length(xo)
        map_ccc(yo(ii),xo(ii)) = cluster(ii,1);

    end

    set(0,'DefaultFigureRenderer','zbuffer') % une ligne de code necessaire car sinon parfois ca bug au moment de faire la figure (#JeSaisPasPourquoi)
    m_proj('mercator','latitude',[lat(1) lat(end)],'longitude',[lon(1) lon(end)]); % Projection a appliquer et bien specififier les latitudes limites


    %ncol = [col(1,:); col(2,:); col(4,:); col(3,:); col(5,:); col(6,:)];
    colormap(col)
    [X,Y] = m_ll2xy(LON,LAT);
    h = pcolor(X,Y,map_ccc);
    set(h,'edgecolor','none')

    tb = m_etopo2('contour',[-5000:500:-500],'edgecolor',[0 0 0]);
    set(gca,'Clim',[min(map_ccc(:)) max(map_ccc(:))])

    m_grid()
    m_gshhs_l('patch',[0.8 0.8 0.8],'edgecolor','k');
    gt = strcat(c_year, ' FCM');
    title(gt)
    pic_name = strcat('C:/Files/Work/Bigelow/Data/figures/Annual/',c_year,'FCMmap.png');
    saveas(gcf,pic_name)
end