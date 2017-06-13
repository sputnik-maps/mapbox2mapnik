## mapbox2mapnik

mapbox2mapnik is a plugin for [Mapnik library](https://github.com/mapnik/mapnik)
which allows it to use [Mapbox style format](https://www.mapbox.com/mapbox-gl-js/style-spec/).

## Dependencies

* libmapnik >= 3.0.x
* libglog
* libjsoncpp
* Boost
    - optional
    - filesystem

## Building

```bash
mkdir build
cd build
cmake ..
make
```

## Installation

```bash
sudo make install
```
