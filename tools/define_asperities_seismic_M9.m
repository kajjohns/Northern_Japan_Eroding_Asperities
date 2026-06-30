function [asp_centers, asp_points, asp_xy, centerx, centery, centers] = define_asperities_seismic_M9(origin, centroids)
%DEFINE_ASPERITIES_SEISMIC_M9 Define seismic asperity outlines and centers.
%
% Polygon asperities are digitized rupture outlines from seismic inversions.
% The polygon points are stored in local coordinates relative to each
% asperity center. Centers listed as lon/lat are converted into the same
% local coordinate system used by the mesh.
%
% Outputs:
%   asp_centers - 1xN cell array of polygon-center [x y] coordinates.
%   asp_points  - 1xN cell array of polygon boundary points around centers.
%   asp_xy      - Mx2 array of circular asperity centers in local [x y].
%   centerx,
%   centery,
%   centers     - mesh centroid coordinates used by locking/ring functions.

asp_centers = {};
asp_points = {};

% Polygon asperities 1-13. These are stored with positive y downward in the
% digitized source, so the y-coordinate is flipped before use.
polygon_specs = {
    'Miyagi 1', [142 38.47], [
        11.84 -15.54
        28.86 -10.36
        32.19 0.74
        28.86 12.95
        27.01 18.13
        25.16 21.09
        12.58 18.87
        2.22 16.65
        -10.36 15.17
        -19.61 11.1
        -24.79 9.25
        -35.15 4.07
        -34.78 -1.85
        -33.3 -7.03
        -28.12 -15.54
        -22.57 -20.35
        -17.76 -14.8
        -9.25 -15.91
        0.74 -19.24];

    'Miyagi 2', [141.73 38.16], [
        -5.18 -16.65
        2.59 -17.39
        15.54 -16.28
        17.76 -8.14
        19.24 0
        18.87 6.66
        12.21 11.47
        1.11 13.69
        -8.88 8.88
        -16.28 -0.37
        -10.73 -5.18
        -7.77 -11.84];

    'Miyagi 3', [142.25 38.17], [
        0 -11.84
        8.14 -7.77
        10.36 0.37
        3.7 10.36
        -2.22 11.1
        -11.47 8.88
        -15.91 3.7
        -14.8 -4.44
        -11.47 -8.51
        -4.81 -11.84];

    'Miyagi 4', [142.21 37.92], [
        2.22 -33.3
        12.58 -31.82
        22.2 -28.86
        24.79 -17.39
        20.35 -5.55
        17.02 3.7
        12.58 14.06
        6.66 19.98
        -2.96 19.61
        -13.69 15.54
        -28.49 8.14
        -34.78 -1.48
        -32.56 -8.88
        -31.08 -16.65
        -25.16 -30.71
        -18.5 -39.59
        -8.88 -38.85];

    'Miyagi 5a', [142.78 38.1], [
        0 -15.91
        10.36 -12.58
        20.35 -5.18
        21.46 4.81
        11.84 12.58
        -0.37 14.8
        -11.1 14.43
        -17.39 6.29
        -16.28 -2.59
        -12.95 -9.99
        -8.14 -12.95];

    'Miyagi 5b', [142.75 37.81], [
        0 -17.76
        12.21 -10.36
        18.5 1.85
        15.17 9.99
        8.51 15.17
        0 17.02
        -12.58 17.02
        -25.53 4.81
        -24.05 -6.29
        -21.83 -11.84
        -8.14 -17.02];

    'Miyagi 6', [142.88 38.38], [
        0 -22.2
        4.44 -17.39
        10.73 -17.02
        19.61 -21.83
        22.2 -14.8
        19.98 -4.81
        15.17 0.74
        19.24 8.88
        18.87 17.76
        14.43 25.16
        7.77 24.42
        -5.18 24.79
        -19.98 26.27
        -24.42 18.5
        -22.57 12.95
        -18.87 6.29
        -13.69 -7.4
        -12.95 -14.06
        -17.76 -21.09
        -16.65 -28.12
        -9.62 -30.71
        -2.22 -27.75
        -1.48 -25.9];

    'Fukushima 1', [141.78 37.6], [
        -1.027777778 -21.58333333
        14.38888889 -13.36111111
        25.69444444 -5.138888889
        24.66666667 11.30555556
        11.30555556 45.22222222
        2.055555556 55.5
        -1.027777778 43.16666667
        -5.138888889 27.75
        -16.44444444 9.25
        -30.83333333 -7.194444444
        -28.77777778 -20.55555556
        -17.47222222 -25.69444444];

    'Fukushima 2', [141.36 37.2], [
        0 -27.75
        15.41666667 -31.86111111
        27.75 -21.58333333
        30.83333333 -7.194444444
        26.72222222 8.222222222
        21.58333333 21.58333333
        7.194444444 23.63888889
        -11.30555556 17.47222222
        -18.5 6.166666667
        -14.38888889 -9.25
        -7.194444444 -21.58333333];

    'Fukushima 3', [141.51 36.57], [
        1.027777778 -28.77777778
        16.44444444 -22.61111111
        33.91666667 -18.5
        34.94444444 -1.027777778
        25.69444444 16.44444444
        9.25 32.88888889
        -17.47222222 23.63888889
        -43.16666667 12.33333333
        -41.11111111 3.083333333
        -35.97222222 -10.27777778
        -22.61111111 -33.91666667
        -11.30555556 -33.91666667];

    'Fukushima 4', [142.19 36.88], [
        -1.027777778 -46.25
        16.44444444 -39.05555556
        42.13888889 -25.69444444
        34.94444444 4.111111111
        22.61111111 32.88888889
        3.083333333 45.22222222
        -20.55555556 35.97222222
        -32.88888889 25.69444444
        -31.86111111 6.166666667
        -27.75 -16.44444444
        -19.52777778 -37
        -11.30555556 -49.33333333];

    'Sanriku 1', [142.94 40.35], [
        -5.138888889 -39.05555556
        10.27777778 -24.66666667
        40.08333333 -17.47222222
        49.33333333 12.33333333
        47.27777778 42.13888889
        28.77777778 75.02777778
        7.194444444 90.44444444
        -17.47222222 74
        -37 35.97222222
        -52.41666667 3.083333333
        -45.22222222 -16.44444444
        -25.69444444 -30.83333333];

    'Sanriku 2', [142.31 40.98], [
        -1.027777778 -66.80555556
        11.30555556 -71.94444444
        22.61111111 -44.19444444
        31.86111111 -41.11111111
        45.22222222 -49.33333333
        53.44444444 -33.91666667
        57.55555556 -5.138888889
        51.38888889 12.33333333
        48.30555556 26.72222222
        26.72222222 41.11111111
        -3.083333333 61.66666667
        -19.52777778 35.97222222
        -29.80555556 11.30555556
        -43.16666667 -37
        -52.41666667 -58.58333333
        -49.33333333 -80.16666667
        -30.83333333 -67.83333333
        -21.58333333 -60.63888889
        -11.30555556 -62.69444444];
};

for k = 1:size(polygon_specs, 1)
    [asp_centers{k}, asp_points{k}] = polygonAsperity( ...
        polygon_specs{k, 2}, polygon_specs{k, 3}, origin, true);
end

% Polygon asperity 14: the M9 outline is already stored in local
% coordinates in M9_points.mat, so it does not use lon/lat conversion.
load(fullfile('data', 'M9_points.mat'))
m9_x = mirror_m9_x - cent_m9_x;
m9_y = mirror_m9_y - cent_m9_y;
m9_points_a = [m9_x(1:33), m9_y(1:33)];
m9_points_b = [m9_x(35:end), m9_y(35:end)];
asp_centers{14} = [cent_m9_x cent_m9_y];
asp_points{14} = [
    m9_points_a(:, 1), m9_points_a(:, 2)
    m9_points_b(:, 1), m9_points_b(:, 2)];

% Mesh centroid coordinates used by the polygon/circular asperity tests.
centerx = centroids(:, 1);
centery = centroids(:, 2);
centers = centroids(:, 3);

% Circular asperities. The first center is already in local coordinates;
% the remaining centers are stored as lon/lat and converted below.
asp_xy = [
    -10 -85
    localxy([141.52 36.31], origin)   % southern asperity for 2008.35 eq
    localxy([141.52 36.31], origin)   % 2005.8 M6.3
    localxy([142.6528 39.5], origin)  % 2010.5 M6.4
    localxy([141.57 + 0.35 38.85], origin)]; % 2003.39

end

function [center_xy, points] = polygonAsperity(center_lonlat, points, origin, flip_y)
center_xy = localxy(center_lonlat, origin);

if flip_y
    points(:, 2) = -points(:, 2);
end
end

function xy = localxy(lonlat, origin)
xy = llh2local(lonlat', fliplr(origin))';
xy = [xy(1) xy(2)];
end
