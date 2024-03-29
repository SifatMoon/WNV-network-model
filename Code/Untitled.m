figure
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
faceColors = makesymbolspec('Polygon',...
    {'INDEX', [1 numel(states)], 'FaceColor', ... 
    polcmap(numel(states))}); %NOTE - colors are random
geoshow(ax, states, 'DisplayType', 'polygon', ...
   'SymbolSpec', faceColors)
framem off; gridm off; mlabel off; plabel off;
colorbar('southoutside')
sifu=1;