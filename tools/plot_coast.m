function plot_coast(bbox)
%bbox [minlon minlat;maxlon maxlat]

tools_dir = fileparts(mfilename('fullpath'));
coastline_file = fullfile(tools_dir, 'ne_10m_coastline', 'ne_10m_coastline.shp');
S = shaperead(coastline_file,'BoundingBox',bbox);

for k=1:length(S)
    
    plot(S(k).X,S(k).Y,'k')
    
end
