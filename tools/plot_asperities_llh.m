
for k=1:length(asp_centers)
    
    points = asp_centers{k} + asp_points{k};
    points_llh = local2llh(points',fliplr(origin))';
    points_llh(end+1,:) = points_llh(1,:);

    asp_llh{k} = points_llh;

    zshift = 100*ones(size(points_llh,1),1);
    plot3(points_llh(:,1),points_llh(:,2),zshift,'Color',col)

end


