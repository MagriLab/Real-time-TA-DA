function [p z] = zoomPlot(x,y,xbounds,pos,varargin)
% Please retain the following:
% 
% Original Author: 
% Kelsey Bower, kelsey.bower@case.edu

if nargin > 5
    printf('Too many arguments. zoomPlot(x,y,xbounds,pos,vertex)\n')
elseif nargin < 5
    vertex = [1 4];
elseif nargin == 5
    vertex = varargin{1};
end
% Get current axis position and limits
p = gca;

% Calculate x,y points of zoomPlot
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = (pos(1)+pos(3)-p.Position(1))/p.Position(3)*diff(p.XLim)+(p.XLim(1));
y1 = (pos(2)-p.Position(2))/p.Position(4)*diff(p.YLim)+p.YLim(1);
y2 = ((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(p.YLim)+p.YLim(1);

% Plot lines connecting zoomPlot to original plot points
index = find(x>=xbounds(1) & x<=xbounds(2)); % Find indexes of points in zoomPlot
rectangle('Position',[xbounds(1) min(y(index)) diff(xbounds) max(y(index))-min(y(index))]);
hold on
if any(vertex==1)
    plot([xbounds(1) x1], [max(y(index)) y2], 'k'); % Line to vertex 1
end
if any(vertex==2)
    plot([xbounds(2) x2], [max(y(index)) y2], 'k'); % Line to vertex 2
end
if any(vertex==3)
    plot([xbounds(2) x2], [min(y(index)) y1], 'k'); % Line to vertex 4
end
if any(vertex==4)
    plot([xbounds(1) x1], [min(y(index)) y1], 'k'); % Line to vertex 3
end

% Plot zoomPlot and change axis
z = axes('position',pos);
box on 
plot(x,y,'c','linewidth',1,'color',[.8 .1 0])
axis([xbounds(1) xbounds(2) min(y(index)) max(y(index))]);
