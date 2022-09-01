# Qmesh_extent_high_resolution
Hi! This is an extension for Qmesh. While Qmesh takes the shoreline as the region of interest,
this package provides function to take another rectangular extent as the chosen place. Comparing to Qmesh,
this package mainly applying the PCA smoothing method to the shoreline.
# Requirement

  * Python 3.x
  * Ubuntu 
  * Qgis installed in Ubuntu
  * Gmsh installed in Ubuntu
  * Qmesh3 installed in Ubuntu
# Installation 
To install Qmesh_extent_high_resolution, please make sure Qmesh3 is running properly in your environment, and the package is usually running on the ubuntu.
If you happen to not have Qmesh3 installed in your package, you can use command line: <br>
```
$ pip install qmesh3 <br>
```
To successfully install the Qmesh package, Qgis and Gmsh are required. Command line below can be used.
```
$ sudo apt install qgis
$ sudo apt install gmsh
```
Then just download this repo to your local computer. 
# How to use this repo
After downloading this repo, you have your local environment like:
```
- your working repo
-- Qmesh_extent_high_resolution
-- mesh_generation.py
-- enclosed.shp
-- enclosed.cpg
-- enclosed.prj
-- enclosed.shx
-- enclosed.dbf
-- unenclosed.shp
....
```
Where mesh_generation.py and example shapefiles can be find in the folder example/case1.
Afterwards, you can just simply run
```
python3 mesh_generation.py
```
A .msh file will be generated, which can be opened by gmsh to see the mesh generated.
# parameter description

```
# read the name
enclosed_filename = 'complex_shorelin_with_manual_boundary.shp'
unencloed_filename = 'complex_shoreline.shp'
# output filename for check for self-intersection
target_enclosed = "target_enclosed"
target_unenclosed = "target_unenclosed"
contour = qmesh.contour.Contour(enclosed_filename, unencloed_filename, target_enclosed, target_unenclosed)
extent_of_interest = (-2.4846, -2.3198, 59.3140, 59.4364)
```
Filename is what the input shape filename, and target is the filename for output. Extent of interest is the extent where you put the place of interest,
not the extent for whole map. 
```
generate_mesh(mesh_name='mesh_generation_complex', # Name the generated files
enclosed_name='target_enclosed.shp', # Enclosed shapefile’s name
unenclosed_name='target_unenclosed.shp', # Unenclosed shapefile’s name
EPSG=4326, # Set CRS
extent=(-4.2673085305105776, -1.8445937110834305, 58.2864165276842456, 59.742403168),
resolution=(5000, 5000), # Set raster resolution, not the mesh resolution
# Set (min size, max size, gradation distance, buffer)
# of mesh elements
# The units are degree for the CRS 'EPSG:4326'
gradation=(0.0008, 0.12, 0.3, 0.0), #min, max, gradation distance, buffer
extent_of_interest = extent_of_interest, 
multiple_times = 6) 
```
where multiple times represents the difference between the extent of interest and outside. Explanations for other parameters can be found in 

  * Semi-automated Mesh Generation of the Coastal Oceans with Engineering Structures Using Satellite Data.pdf
Which is the tutorial for Qmesh.

