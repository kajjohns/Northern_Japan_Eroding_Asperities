function ring_tau = get_ringtau_asperities(ring_tau_asp,D,asp_scale,asp_centers,asp_points,asp_xy,Rupture_Scale,centerx,centery,i_lock)

ring_tau = zeros(size(i_lock));

for k=1:length(asp_points)


    asp_points_scale{k} = [asp_centers{k}(1)+(Rupture_Scale(k)*asp_scale(k)+D(k))*asp_points{k}(:,1) asp_centers{k}(2)+(Rupture_Scale(k)*asp_scale(k)+D(k))*asp_points{k}(:,2)];

    i_ring =  inpoly([centerx centery], asp_points_scale{k}) & ~i_lock;
 

    ring_tau(i_ring) = ring_tau_asp(k);

end
