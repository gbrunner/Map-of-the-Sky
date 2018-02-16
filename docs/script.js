require([
  'esri/core/promiseUtils',
  'esri/core/urlUtils',
  'esri/Map',
  'esri/views/MapView',
  'esri/views/SceneView',
  'esri/layers/BaseTileLayer',
  'esri/layers/TileLayer',
  'esri/widgets/Compass',
  'esri/widgets/Home',
  'dojo/domReady!'
], function(
  promiseUtils,
  urlUtils,
  Map,
  MapView,
  SceneView,
  BaseTileLayer,
  TileLayer,
  Compass,
  Home
) {
  var is2DView = (urlUtils.urlToObject(window.location.href).query && !!(Number(urlUtils.urlToObject(window.location.href).query.is2DView)));

  var ViewModule = is2DView ? MapView : SceneView;

  // black base layer for SceneView and MapView adapted from @ycabon's
  // https://codepen.io/ycabon/pen/gvXqqj?editors=1000
  var BlackLayer = BaseTileLayer.createSubclass({
    constructor: function() {
      var canvas = this.canvas = document.createElement('canvas');
      canvas.width = canvas.height = 256;
      var ctx = canvas.getContext('2d');
      ctx.fillRect(0, 0, 256, 256);
    },
    fetchTile: function(level, row, col) {
      return promiseUtils.resolve(this.canvas);
    }
  });

  var view = new ViewModule({
    container: 'viewDiv',
    map: new Map({
      basemap: {
        baseLayers: [
          new BlackLayer()
        ]
      },
      layers: [
        new TileLayer({
          url: 'https://tiles.arcgis.com/tiles/g2TonOxuRkIqSOFx/arcgis/rest/services/Galactic_Plane_Bright/MapServer'
        })
      ]
    }),
    center: [-14.5, -0.36],
    zoom: 10,
    constraints: {
      minZoom: 4
    }
  });

  if (!is2DView) {
    view.environment.atmosphereEnabled = false;
  }

  if (is2DView) {
    view.ui.add(new Compass({
      view: view
    }), 'top-left');
  }

  view.ui.add(new Home({
    view: view
  }), 'top-left');

  view.when(function() {
    var switchViewNode = is2DView ? document.getElementById('switchViewTo3D') : document.getElementById('switchViewTo2D');
    view.ui.add(switchViewNode, 'top-left');
    switchViewNode.style.display = 'flex';
  });
});
