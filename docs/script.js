require([
  'esri/core/promiseUtils',
  'esri/core/urlUtils',

  'esri/Basemap',
  'esri/Map',
  'esri/views/MapView',
  'esri/views/SceneView',
  'esri/layers/BaseTileLayer',
  'esri/layers/TileLayer',

  'esri/widgets/BasemapGallery',
  'esri/widgets/Compass',
  'esri/widgets/Expand',
  'esri/widgets/Home',

  'dojo/domReady!'
], function(
  promiseUtils,
  urlUtils,

  Basemap,
  Map,
  MapView,
  SceneView,
  BaseTileLayer,
  TileLayer,

  BasemapGallery,
  Compass,
  Expand,
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

  var blackLayer = new BlackLayer();

  var basemapsCollection = [
    new Basemap({
      title: 'Galactic Plane',
      baseLayers: [
        new TileLayer({
          url: 'https://tiles.arcgis.com/tiles/g2TonOxuRkIqSOFx/arcgis/rest/services/Galactic_Plane/MapServer'
        })
      ],
      thumbnailUrl: 'https://tiles.arcgis.com/tiles/g2TonOxuRkIqSOFx/arcgis/rest/services/Galactic_Plane/MapServer/tile/10/512/471'
    }),
    new Basemap({
      title: 'Galactic Plane - Faint',
      baseLayers: [
        new TileLayer({
          url: 'https://tiles.arcgis.com/tiles/g2TonOxuRkIqSOFx/arcgis/rest/services/Galactic_Plane_Faint/MapServer'
        })
      ],
      thumbnailUrl: 'https://tiles.arcgis.com/tiles/g2TonOxuRkIqSOFx/arcgis/rest/services/Galactic_Plane_Faint/MapServer/tile/10/512/471'
    }),
    new Basemap({
      title: 'Galactic Plane - Bright',
      baseLayers: [
        new TileLayer({
          url: 'https://tiles.arcgis.com/tiles/g2TonOxuRkIqSOFx/arcgis/rest/services/Galactic_Plane_Bright/MapServer'
        })
      ],
      thumbnailUrl: 'https://tiles.arcgis.com/tiles/g2TonOxuRkIqSOFx/arcgis/rest/services/Galactic_Plane_Bright/MapServer/tile/10/512/471'
    })
  ];

  var view = new ViewModule({
    container: 'viewDiv',
    map: new Map(), // let the basemapGallery widget control the map's basemap value
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

  var basemapGallery = new BasemapGallery({
    source: basemapsCollection,
    view: view
  });

  // the custom blackLayer isn't fully valid with the basemapGallery widget
  // so we re-add it to the bottom of the map's basemap baseLayers
  // each time the user makes a new selection
  basemapGallery.watch('activeBasemap', function() {
    view.map.basemap.baseLayers.add(blackLayer, 0);
  });

  // set the initial basemap value, which will also set it on the map
  basemapGallery.activeBasemap = basemapsCollection[0];

  view.ui.add(new Expand({
    autoCollapse: true,
    content: basemapGallery,
    expandIconClass: 'esri-icon-basemap',
    view: view
  }), {
    position: 'top-right'
  });

  view.when(function() {
    var switchViewNode = is2DView ? document.getElementById('switchViewTo3D') : document.getElementById('switchViewTo2D');
    view.ui.add(switchViewNode, 'top-left');
    switchViewNode.style.display = 'flex';
  });
});
