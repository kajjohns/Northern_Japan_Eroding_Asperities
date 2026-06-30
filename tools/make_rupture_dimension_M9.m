function [date_mag, Rupture_Areas, Rupture_Scale] = make_rupture_dimension_M9(asp_points)
%MAKE_RUPTURE_DIMENSION_M9 Compute rupture areas and asperity scale factors.
%
% date_mag rows correspond to the 14 polygon asperities in asp_points,
% followed by the 5 circular asperities in asp_xy. Rupture areas are
% estimated from magnitude using the Wells and Coppersmith (1994) reverse
% or thrust rupture scaling relationship.

date_mag = [
    % Miyagi earthquakes. Rows 1 and 2 are subregions of the 1978 event.
    1978     7.5
    1978     7.5
    2005     7.2
    1936     7.5
    2003     7.0
    2003     7.0
    1981     7.1

    % Fukushima earthquakes.
    1938     7.7
    1938     7.7
    1938     7.7
    1938     7.7

    % Sanriku earthquakes.
    1994     7.9
    1994     7.7

    % Tohoku M9 polygon asperity.
    1454     9.1

    % Circular asperities.
    2008.55  6.9
    2008.35  6.8
    2005.8   6.3
    2010.5   6.4
    2003.38  7.0];

num_polygon_asperities = length(asp_points);
PolyAreas = zeros(num_polygon_asperities, 1);

for k = 1:num_polygon_asperities
    PolyAreas(k) = getPolygonArea(asp_points{k});
end

% Wells and Coppersmith (1994, BSSA):
% Area = 10^(-3.99 + 0.98*Mw), for reverse or thrust ruptures.
Rupture_Areas = 10.^(-3.99 + 0.98*date_mag(:, 2));  % km^2

Rupture_Scale = zeros(size(date_mag, 1), 1);

% Polygon asperities are kept at their digitized outlines.
Rupture_Scale(1:num_polygon_asperities) = 1;

% Circular asperity scale values are equivalent rupture radii in km.
circle_rows = (num_polygon_asperities + 1):size(date_mag, 1);
Rupture_Scale(circle_rows) = sqrt(Rupture_Areas(circle_rows) / pi);

end
