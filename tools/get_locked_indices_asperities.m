function i_lock = get_locked_indices_asperities(asp_scale,asp_centers,asp_points,asp_xy,Rupture_Scale,centerx,centery)

i_lock = false(length(centerx),1);


for k=1:length(asp_points)

    asp_points_scale{k} = [asp_centers{k}(1)+Rupture_Scale(k)*asp_scale(k)*asp_points{k}(:,1) asp_centers{k}(2)+Rupture_Scale(k)*asp_scale(k)*asp_points{k}(:,2)];

    i_lock = i_lock | inpoly([centerx centery], asp_points_scale{k});
    %slip_index{k} = inpoly([centerx centery], asp_points_scale{k});



end


for k=1:size(asp_xy,1)
        
        dist = sqrt((asp_xy(k,1) - centerx).^2 + (asp_xy(k,2) - centery).^2);

        i_lock = i_lock | dist<Rupture_Scale(14+k);

      
end





