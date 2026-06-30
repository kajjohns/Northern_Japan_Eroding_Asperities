function outArea = getPolygonArea(varargin)

% input is a 2-column vector of xy verteces
% optional input is 2-element vector of pixel grid dimensions
%
% this function was written to account for self-intersecting segments
% output area matches the active area used in MATLAB 'inpolygon' function

if nargin>2, error('Incorrect inputs'); end

verts = varargin{1};

if nargin==2
    xDim = varargin{2}(1);
    yDim = varargin{2}(2);
else
    xDim = 1;
    yDim = 1;
end

%close polygon if open
if ~all(verts(1,:)==verts(end,:))
    verts = [verts; verts(1,:)];
end

outArea = 0;
%original segments from verteces
segX = [verts(:,1) [verts(2:end,1); verts(1,1)]];
segY = [verts(:,2) [verts(2:end,2); verts(1,2)]];
%relative direction of each segment
w = sign(segX(:,2) - segX(:,1));
%call 'polyxpoly' to get all x intersection verteces for sorting into vertical slices
[x,~] = polyxpoly(verts(:,1),verts(:,2),verts(:,1),verts(:,2),'unique');
[~, xsrtidx] = sort(x);
for n=1:length(x)-1
    %find segments traversing this slice
    idx = find((segX(:,1)<=x(xsrtidx(n)) & segX(:,2)>=x(xsrtidx(n+1))) | ...
        (segX(:,1)>=x(xsrtidx(n+1)) & segX(:,2)<=x(xsrtidx(n))));
    nsegs = numel(idx); %must be even
    %slice segments
    sliceSegY = zeros(nsegs,2);
    for m=1:nsegs
        %find intersects of line segments with slice edges
        s = (segY(idx(m),2)-segY(idx(m),1))/...
            (segX(idx(m),2)-segX(idx(m),1));
        tempY = [s*(x(xsrtidx(n))-segX(idx(m),1))+segY(idx(m),1); ...
            s*(x(xsrtidx(n+1))-segX(idx(m),1))+segY(idx(m),1)];
        sliceSegY(m,:) = round(tempY'*1e6)/1e6;
    end
    [~, ysrtidx] = sort(mean(sliceSegY,2));
    ysrtidx = ysrtidx(numel(ysrtidx):-1:1);
    %get "winding" of first encountered segment in slice
    weff = w(idx(ysrtidx(1)));
    %calulate area by integration
    xint = (x(xsrtidx(n+1)) - x(xsrtidx(n))) * xDim;
    if xint~=0
        for m=1:nsegs
            if m==1
                sliceArea = xint * w(idx(ysrtidx(m))) * weff * ...
                    mean(sliceSegY(ysrtidx(m),:)) * yDim;
                wcurr = weff;
                linescrossed = 0;
                inpoly = true;
            else
                if w(idx(ysrtidx(m)))==(-1*wcurr)
                    if linescrossed~=0
                        linescrossed = linescrossed - 1;
                    else
                        temp = (xint * w(idx(ysrtidx(m))) * weff * ...
                            mean(sliceSegY(ysrtidx(m),:)) * yDim);
                        sliceArea = sliceArea + temp;
                        wcurr = wcurr * -1;
                        if inpoly
                            inpoly = false;
                        else
                            inpoly = true;
                        end
                    end
                else
                    if ~inpoly
                        temp = (xint * ...
                            mean(sliceSegY(ysrtidx(m),:)) * yDim);
                        sliceArea = sliceArea + temp;
                        weff = weff * -1;
                        inpoly = true;
                    else
                        linescrossed = linescrossed + 1;
                    end
                end
            end
        end
        outArea = outArea + sliceArea;
    end
end

end