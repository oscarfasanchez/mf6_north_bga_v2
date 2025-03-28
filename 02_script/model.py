# code to build the groundwtaer model for the north of bucaramanga


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import flopy
import geopandas as gpd
# import gmshflow v0.1.0
import gmshflow

def create_voronoi_grid(wd_shp):
    #load the shapefiles

    #load domain
    dom = gpd.read_file(os.path.join(wd_shp, 'Domin_Mod.shp'))
    #load rivers
    riv = gpd.read_file(os.path.join(wd_shp, 'Rios.shp'))
    #load creeks
    crks = gpd.read_file(os.path.join(wd_shp, 'queb_corr.shp'))
    #load observations and filter pizometer observations with the column DEPTH_MEA > 0
    obs = gpd.read_file(os.path.join(wd_shp, 'INV_PAS_V5_DEM.shp'))
    obs = obs[obs['DEPTH_MEA'] > 0]
    #get a geoseries with the unique points
    obs_geo = obs.remove_repeated_points(tolerance=1.0).normalize().drop_duplicates()
    obs = gpd.GeoDataFrame(geometry=obs_geo)
    #load Faults
    faults = gpd.read_file(os.path.join(wd_shp, 'Fallas_V2_EditPau.shp'))
    # filter Sinistral Faults
    faults = faults.loc[faults['cinematica'].str.contains('Sinistral')]
    #load Gallery
    gal = gpd.read_file(os.path.join(wd_shp, 'Alineamiento_Galeria.shp'))
    

    # Create the mesh

    # Define cell_size in each area
    cs_dom = 200.0
    cs_crk = 50.0
    cs_obs = 10.0
    cs_fault = 50
    cs_gal = 3

    # Add cell size column
    crks['cs'] = cs_crk
    obs['cs'] = cs_obs
    faults['cs'] = cs_fault
    gal['cs'] = cs_gal


    # set creek only if its less than 300 meters from gal
    crks = crks.loc[crks.distance(gal.geometry[0]) < 300]

    # Initialize GmshModel and GmshMeshDomain
    gmsh_model = gmshflow.GmshModel("north_bga_gmsh")
    mesh_domain = gmshflow.GmshMeshDomain("domain", dom, cs_dom)

    # Prepare mesh domain
    mesh_domain.prepare_mesh_domain(mesh_area=1, min_overlap=0.5)

    # Create domain loop
    c_ind = mesh_domain.create_domain_loop_from_poly()

    # add features to create fields
    crks_handler = gmshflow.LineGeometryHandler()
    crks_handler.set_gpd_line(crks)
    crks_ind_list = crks_handler.create_line_from_line()

    obs_handler = gmshflow.PointGeometryHandler()
    obs_handler.set_gdf_point(obs)
    obs_ind_list = obs_handler.create_point_from_point(df_coord=True)

    faults_handler = gmshflow.LineGeometryHandler()
    faults_handler.set_gpd_line(faults)
    faults_ind_list = faults_handler.create_line_from_line()

    gal_handler = gmshflow.LineGeometryHandler()
    gal_handler.set_gpd_line(gal)
    gal_ind_list = gal_handler.create_line_from_line()

    # Create domain surface
    ind_s_dom = mesh_domain.create_domain_surface()

    # #convert to points for exponential fields
    gdf_crk_coord = crks_handler.convert_to_points_for_size_fields()
    gdf_faults_coord = faults_handler.convert_to_points_for_size_fields()
    gdf_gal_coord = gal_handler.convert_to_points_for_size_fields()

    df_coords = pd.concat([ obs_handler.gdf_coord])#gdf_crk_coord, gdf_faults_coord, gdf_gal_coord,

    # #create exponential fields
    mesh_domain.create_exponential_field(df_coords, fac=1.1)

    # # # Set exponential field as background mesh
    mesh_domain.set_field()

    # # Deactivate the mesh size from the geometry points
    # mesh_domain.set_mesh_size_from_geometries(use_boundaries=False, use_points=False)

    ind_s_dom = mesh_domain.create_domain_surface()
    # Run GUI 
    gmsh_model.run_gui()
    # Generate mesh
    print('Generating the mesh...')
    gmsh_model.generate_mesh()
    print('Mesh ready')

      # Run GUI 
    gmsh_model.run_gui()

    #get some quality stats of the triangular mesh
    df_quality = gmsh_model.get_triangular_quality()
    print(df_quality.gamma.min(), df_quality.gamma.mean())
    df_quality.hist(column='gamma', bins=50)

    # Export to Voronoi
    shp_mesh_name = "gdf_voro_bga"
    gmsh_model.export_to_voronoi(wd_shp, "gdf_voro_bga", [ind_s_dom])


    # Generate mesh
    print('Generating the mesh...')
    gmsh_model.generate_mesh()
    print('Mesh ready')

    #convert shapefile to cvfd to be imported in modflow6
    print('Converting shapefile to cvfd...')
    verts, iverts = flopy.utils.cvfdutil.shapefile_to_cvfd(
                    os.path.join(wd_shp, f'{shp_mesh_name}.shp'), verbose=True, duplicate_decimals=3, skip_hanging_node_check=True)
    gridprops = flopy.utils.cvfdutil.get_disv_gridprops(verts, iverts, xcyc=None)

    vertices = gridprops["vertices"]
    cell2d = gridprops["cell2d"]
    nlay = 1
    ncpl = gridprops["ncpl"]
    idomain = np.ones((nlay, ncpl), dtype=int)

    modelgrid = flopy.discretization.VertexGrid(
        vertices=vertices,
        cell2d=cell2d,
        idomain=idomain,
        nlay=nlay,
        ncpl=ncpl,
    )

    return modelgrid


def define_vertical_discretizations(wd_raster):
    pass   
    
    

def create_mf6_model(ws):
    # Create the Flopy simulation object
    sim = flopy.mf6.MFSimulation(sim_name='north_bga_sim', sim_ws=ws, exe_name='mf6', version='mf6')
    gwf = flopy.mf6.ModflowGwf(sim, modelname='north_bga_model', save_flows=True)

    wd_shp = os.path.join('..', '01_GIS',  'shp') 
    wd_raster = os.path.join('..', '01_GIS',  'raster')
    # Define the spatial discretization
    modelgrid = create_voronoi_grid(wd_shp)
    define_vertical_discretizations(wd_raster)
    #define the temporal discretization

    #define the initial conditions

    #define the boundary conditions

    #define the RIVER package

    #define the ghb package

    #define the drain package for creeks

    #define the drn package for the gallery

    #define the drn package for horizontal drains

    


    return sim


if __name__ == '__main__':
    ws = os.path.join('..','03_models', 'base')
    create_mf6_model(ws)
