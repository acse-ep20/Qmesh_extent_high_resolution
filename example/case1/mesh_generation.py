import Qmesh_extent_high_resolution as qmesh
# read the name
enclosed_filename = 'complex_shorelin_with_manual_boundary.shp'
unencloed_filename = 'complex_shoreline.shp'
# output filename for check for self-intersection
target_enclosed = "target_enclosed"
target_unenclosed = "target_unenclosed"

contour = qmesh.contour.Contour(enclosed_filename, unencloed_filename, target_enclosed, target_unenclosed)
extent_of_interest = (-2.4846, -2.3198, 59.3140, 59.4364)
contour.simplify_the_feature(extent_of_interest, points_per_group=20, overlapping_points=19, num_modes=1)

def generate_mesh(mesh_name, enclosed_name, unenclosed_name, EPSG, extent,
resolution, gradation, extent_of_interest, multiple_times):
    print('generating mesh '+mesh_name)
    qmesh.initialise()
    enclosed = qmesh.vector.Shapes()
    enclosed.fromFile(enclosed_name)
    loopShapes = qmesh.vector.identifyLoops(enclosed,
    isGlobal=False, defaultPhysID=1000, fixOpenLoops=True)
    polygonShapes = qmesh.vector.identifyPolygons(loopShapes,
    smallestNotMeshedArea=100, meshedAreaPhysID=99)
    # Apply one gradation
    unenclosed = qmesh.vector.Shapes()

    unenclosed.fromFile(unenclosed_name)
    grad = qmesh.raster.meshMetricTools.gradationToShapes()
    grad.setShapes(unenclosed)
    grad.setRasterBounds(extent[0], extent[1], extent[2], extent[3])
    grad.setRasterResolution(resolution[0], resolution[1])
    grad.setGradationParameters(
    gradation[0], gradation[1], gradation[2], gradation[3])
    grad.calculateLinearGradationWithExtent(extent, multiple_times)
    # grad.calculateLinearGradation()
    domain = qmesh.mesh.Domain()
    domain.setGeometry(loopShapes, polygonShapes)
    domain.setMeshMetricField(grad)
    domain.setTargetCoordRefSystem('EPSG:'+str(EPSG), fldFillValue=1000.0)
    domain.gmsh(geoFilename=mesh_name+'.geo',
        fldFilename=mesh_name+'.fld',
        mshFilename=mesh_name+'.msh')
    mesh = qmesh.mesh.Mesh()
    mesh.readGmsh(mesh_name+'.msh', 'EPSG:'+str(EPSG))
    mesh.writeShapefile(mesh_name)

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
