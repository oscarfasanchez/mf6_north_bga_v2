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

def create_disv_voronoi_grid(sim, wd_shp, wd_raster, coarse_cs=200.0, medium_cs=30, fine_cs=5,
                                    smooth_factor=1.1, bottom_model_offset=5,
                                      geol_layer_division=[2, 2, 2, 2, 1], plot=False):
    gwf = sim.get_model()
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
    cs_dom = coarse_cs
    cs_crk = medium_cs
    cs_obs = 10.0
    cs_fault = 50
    cs_gal = fine_cs

    # Add cell size column
    crks['cs'] = cs_crk
    obs['cs'] = cs_obs
    faults['cs'] = cs_fault
    gal['cs'] = cs_gal


    # set creek only if its less than 500 meters from gal
    crks = crks.loc[crks.distance(gal.geometry[0]) < 500]

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

    df_coords = pd.concat([
        obs_handler.gdf_coord, gdf_faults_coord, gdf_crk_coord,
        gdf_gal_coord])#obs_handler.gdf_coord ,gdf_crk_coord, gdf_faults_coord, gdf_gal_coord,

    # #create exponential fields
    mesh_domain.create_exponential_field(df_coords, fac=smooth_factor)

    # # # Set exponential field as background mesh
    mesh_domain.set_field()

    # # Deactivate the mesh size from the geometry points
    # mesh_domain.set_mesh_size_from_geometries(use_boundaries=False, use_points=False)

    ind_s_dom = mesh_domain.create_domain_surface()
    
    # Generate mesh
    print('Generating the mesh...')
    gmsh_model.generate_mesh()
    print('Mesh ready')

    
    if plot:
        # Run GUI 
        gmsh_model.run_gui()

        #get some quality stats of the triangular mesh
        df_quality = gmsh_model.get_triangular_quality()
        print('Min and mean gamma:',df_quality.gamma.min(), df_quality.gamma.mean())
        df_quality.hist(column='gamma', bins=50)

    # Export to Voronoi
    shp_mesh_name = "gdf_voro_bga"
    mesh_domain.export_to_voronoi(wd_shp, "gdf_voro_bga", [ind_s_dom])


    # Generate mesh
    print('Generating the mesh...')
    gmsh_model.generate_mesh()
    print('Mesh ready')

    #convert shapefile to cvfd to be imported in modflow6
    print('Converting shapefile to cvfd...')
    verts, iverts = flopy.utils.cvfdutil.shapefile_to_cvfd(
                    os.path.join(wd_shp, f'{shp_mesh_name}.shp'), verbose=True,
                    duplicate_decimals=3, skip_hanging_node_check=True)
    gridprops = flopy.utils.cvfdutil.get_disv_gridprops(verts, iverts, xcyc=None)


    ncpl = gridprops["ncpl"]
    nlay = len(geol_layer_division)#sum(geol_layer_division)  # Number of layers based on geol_layer_division
    # Create a simple idomain array with all cells active
    # This can be modified later to include inactive cells if needed
    idomain = np.ones((nlay, ncpl), dtype=int)
    
    modelgrid = flopy.discretization.vertexgrid.VertexGrid(
        vertices=gridprops["vertices"],
        cell2d=gridprops["cell2d"],
        idomain=idomain,
        nlay=nlay,
        ncpl=ncpl,)

    #now lets define the layers

    raster_names = ["R_Topo_resamp.tif",'R_Qd.tif', 'R_Qbg.tif', "R_Qbo2.tif", "R_Qbo1.tif"]
    raster_data = {}
    # Load the raster data for topography and other properties
    for raster_name in raster_names:
        raster_path = os.path.join(wd_raster, raster_name)
        if not os.path.exists(raster_path):
            raise FileNotFoundError(f"Raster file {raster_path} does not exist.")
        # Load the raster data
        raster = flopy.utils.Raster.load(raster_path)
        raster_data[raster_name] = raster.resample_to_grid(modelgrid, method='linear', band=raster.bands[0])
    
    assert geol_layer_division== [2, 2, 2, 2, 1], "geol_layer_division must be [2, 2, 2, 2, 1] or the code should be modified to handle different divisions"
    # Qd
    # botm_ly1 = raster_data['R_Qd.tif'] + (raster_data["R_Topo_resamp.tif"] - raster_data['R_Qd.tif'])*(1 / geol_layer_division[0])
    botm_ly2 = raster_data['R_Qd.tif'] 
    
    #Qbg
    # botm_ly3 = raster_data['R_Qbg.tif'] + (raster_data["R_Qd.tif"] - raster_data['R_Qbg.tif'])*(1 / geol_layer_division[1])
    botm_ly4 = raster_data['R_Qbg.tif']

    #Qbo2
    # botm_ly5 = raster_data['R_Qbo2.tif'] + (raster_data["R_Qbg.tif"] - raster_data['R_Qbo2.tif'])*(1 / geol_layer_division[2])
    botm_ly6 = raster_data['R_Qbo2.tif']

    #Qbo1
    # botm_ly7 = raster_data['R_Qbo1.tif'] + (raster_data["R_Qbo2.tif"] - raster_data['R_Qbo1.tif'])*(1 / geol_layer_division[3])
    botm_ly8 = raster_data['R_Qbo1.tif'] 

    # Bottom layer -Rock
    min_value = np.minimum.reduce([botm_ly2, botm_ly4, botm_ly6, botm_ly8]).min()
    botm_ly9 = np.ones(ncpl) * (min_value - bottom_model_offset)  # Adjusting the bottom

    print(f"Bottom model elevation: {botm_ly9}")

    

    #create the botm array based on the geol_layer_division
    botm = np.array([
         botm_ly2, botm_ly4, botm_ly6, botm_ly8, botm_ly9    #botm_ly1, , botm_ly3, botm_ly5 , botm_ly7
    ])

    top=raster_data["R_Topo_resamp.tif"]
    
    # get the area of each raster where its not equal to raster.nodatavals
    active = np.where(botm < 1e30, 1, -1)#raster.nodatavals[0] there are two???

    #IDOMAIN required to handle negative thickness,raster holes and inactive cells, also where the raster is not defined
    # calculate the thickness array
    thickness = np.zeros((nlay, ncpl), dtype=float)

    botm[0, active[0] == -1] = top[active[0] == -1]  # Adjust bottom elevation for inactive cells 
    thickness[0, :] = top[:] - botm[0, :]
    idomain[0, thickness[0, :] <= 0] = -1  # Set pass through cells only in first layer
    botm[0, thickness[0, :] <= 0] = top[thickness[0, :] <= 0]  # Adjust bottom elevation for pass through cells, that are
    #recalculate thickness 
    thickness[0, thickness[0, :] <= 0] = 0.0  # Set thickness to zero for pass through cells

    for i in range(1, nlay):
        botm[i, active[i] == -1] = botm[i - 1, active[i] == -1]  # Adjust bottom elevation for inactive cells
        thickness[i, :] = botm[i - 1, :] - botm[i, :] 
        # Set idomain based on active cells
        idomain[i , thickness[i, :] <= 0] = -1  # Set pass through cells to inactive
        # Set the last layer's thickness
        botm[i, thickness[i, :] <= 0] = botm[i - 1, thickness[i, :] <= 0]  # Adjust bottom elevation for pass through cells
        # Set the last layer's thickness to a minimum value
        thickness[i, thickness[i, :] <= 0] = 0.0  # Set thickness to zero for pass through cells

    # now lets set the missing layers by diving the original geological layers
    #and duplicating the idomain
    nlay = sum(geol_layer_division)  # Total number of layers after division
    botm_div = np.zeros((nlay, ncpl), dtype=int)
    idomain_div = np.zeros((nlay, ncpl), dtype=int)

    botm_div[0, :] = botm[0, :]  + (top - botm[0, :]) * (1 / geol_layer_division[0])
    idomain_div[0, :] = idomain[0, :]  

    botm_div[1, :] = botm[0, :]
    idomain_div[1, :] = idomain[0, :]

    botm_div[2, :] = botm[1, :]  + (botm[0, :] - botm[1, :]) * (1 / geol_layer_division[1])
    idomain_div[2, :] = idomain[1, :]

    botm_div[3, :] = botm[1, :]
    idomain_div[3, :] = idomain[1, :]

    botm_div[4, :] = botm[2, :]  + (botm[1, :] - botm[2, :]) * (1 / geol_layer_division[2])
    idomain_div[4, :] = idomain[2, :]

    botm_div[5, :] = botm[2, :]
    idomain_div[5, :] = idomain[2, :]

    botm_div[6, :] = botm[3, :]  + (botm[2, :] - botm[3, :]) * (1 / geol_layer_division[3])
    idomain_div[6, :] = idomain[3, :]

    botm_div[7, :] = botm[3, :]
    idomain_div[7, :] = idomain[3, :]

    botm_div[8, :] = botm[4, :]
    idomain_div[8, :] = idomain[4, :]






    # # get the cells in the active domain, its not necessary when the cells has been excluded already from the model
    # shp_dom = dom.geometry[0]
    # ix = flopy.utils.GridIntersect(modelgrid, method='vertex')
    # domain_cells = ix.intersect(shp_dom).cellids
 

    disv = flopy.mf6.ModflowGwfdisv(
        gwf,
        nlay=nlay,
        **gridprops, # Unpack grid properties
        top=top,  # Top elevation of the model
        botm=botm_div,  # Bottom elevation for each layer
        idomain=idomain_div,  # Domain indicator array
    )
    if plot:
        disv.plot()

    return disv

def define_temporal_discretization(sim, totsim_days=365):
    # time in seconds

    nsper = totsim_days  # number of stress periods
    perlen = 86400.0  # 1 day in seconds
    nstp = 1  # number of time steps per period 
    tsmult = 1.0  # time step multiplier

    time_disc = [(perlen, nstp, tsmult) for _ in range(nsper-1)]#[perlen, nstp, tsmult]
    time_disc.insert(0,(1,1,1.0))#inserting the steady stress period at the beginning of list
    tdis= flopy.mf6.ModflowTdis(sim, pname="tdis",
                              time_units="SECONDS", 
                              nper=nsper, perioddata=time_disc)
    period_ats =[(i, 86400, 1.0e-5, 86400, 2.0, 5.0) for i in range(1, nsper)]
    ats = flopy.mf6.ModflowUtlats(tdis, maxats=len(period_ats), perioddata= period_ats, pname="ats"   )

    return tdis, ats




def define_solver(sim):   
    ims = flopy.mf6.ModflowIms(sim, pname="ims",
                                          complexity= "MODERATE",print_option="SUMMARY",
                                          outer_maximum=350, outer_dvclose=0.1,
                                          under_relaxation="DBD", under_relaxation_gamma=0.1,
                                          under_relaxation_theta=0.7, under_relaxation_kappa=0.1,
                                          backtracking_number=20, backtracking_tolerance=20,
                                          backtracking_reduction_factor=0.1, backtracking_residual_limit=0.002, 
                                          inner_maximum=500, inner_dvclose=0.01,
                                          linear_acceleration="bicgstab", preconditioner_levels=20,
                                          preconditioner_drop_tolerance=0.0001,ats_outer_maximum_fraction=0.05,
                                          scaling_method="DIAGONAL", reordering_method="MD",
                                          number_orthogonalizations=2)
    return ims

def create_npf_package(sim, kh, kv, geol_layer_division, plot=False):
    gwf = sim.get_model()
    # Define the conductivity for each layer
    # Ensure that the length of kh and kv matches the number of layers
    kh_div=[]
    kv_div=[]
    for glay in range(len(geol_layer_division)):
        for div in range(geol_layer_division[glay]):
            kh_div.append(kh[glay])
            kv_div.append(kv[glay])

    kh_div = np.array(kh_div)
    kv_div = np.array(kv_div)
    # Create the NPF package with specified parameters
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=kh_div, k33overk=True, k33=kv_div,
                                save_flows=True, save_specific_discharge=True)
    if plot:
        npf.plot()
    return npf

def create_river_chd_package(sim, wd_shp):
    gwf = sim.get_model()
    # lets get the model top elevation
    top = gwf.dis.top.array
    nlay = gwf.dis.nlay.array

    # Load the river shapefile
    riv_shp = gpd.read_file(os.path.join(wd_shp, 'Rios.shp'))
    # lets merge all the element in the river shapefile into one geometry,
    # because we dont have special data on it
    ix = flopy.utils.GridIntersect(gwf.modelgrid, method='vertex')
    # Intersect the river shapefile with the model grid
    riv_geo = riv_shp.geometry.union_all()
    riv_intersect = ix.intersect(riv_geo)
    chd_spd = []
    for i in range(riv_intersect.shape[0]):
        # Get the cell ids and the top elevation for each intersected cell
        cellid = riv_intersect["cellids"][i]
        #let's check if the cellid is active if its -1 we apply to the cell below, or skip if inactive
        # Look for the upper cell active
        
        for lay in range(nlay):
            if gwf.dis.idomain.array[lay, cellid] == 1:
                chd_spd.append([lay, *cellid, top[cellid] + 1])#add one meter of water to terrain
                break
    
    chd_riv  = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd,
                                       pname="chd", save_flows=True,
                                       print_input=False, print_flows=True,)

    return chd_riv
def create_ghb_package_lat_inflow(sim, wd_shp):
    gwf = sim.get_model()
    # Load the GHB shapefile
    ghb_shp = gpd.read_file(os.path.join(wd_shp, 'Chd_In.shp'))
    ghb_spd = []
    ghb_geo = ghb_shp.geometry#there is only one geometry in the shapefile
    ix = flopy.utils.GridIntersect(gwf.modelgrid, method='vertex')
    # Get the cell ids and the elevation for each GHB
    ghb_intersect = ix.intersect(ghb_geo)
    for i in range(ghb_intersect.shape[0]):
        # Get the cell ids and the top elevation for each intersected cell
        cellid = ghb_intersect["cellids"][i]
        #let's check if the cellid is active
        for lay in range(gwf.dis.nlay.array):
            if gwf.dis.idomain.array[lay, cellid] == 1:
                ghb_spd.append([lay, *cellid, ghb_intersect["elev"][i]]) 


    ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=ghb_spd,
                                   pname="ghb", save_flows=True,
                                   print_input=False, print_flows=True)
    return ghb

def create_mf6_model(ws):
    # Create the Flopy simulation object
    sim = flopy.mf6.MFSimulation(sim_name='north_bga_sim', sim_ws=ws,
                                  exe_name='mf6', version='mf6')
    gwf = flopy.mf6.ModflowGwf(sim, modelname='north_bga_model',
                            newtonoptions="UNDER_RELAXATION",
                              save_flows=True)

    wd_shp = os.path.join('..', '01_GIS',  'shp') 
    wd_raster = os.path.join('..', '01_GIS',  'raster')


    ims = define_solver(sim)
    
    #define the temporal discretization
    tdis, ats = define_temporal_discretization(sim)
    # Define the spatial discretization1
    #define the cell sizes inside the voronoi function
    #TODO define minimum thickness for each layer
    geol_layer_division = [2, 2, 2, 2, 1]  # Geological layer division
    disv = create_disv_voronoi_grid(sim, wd_shp, wd_raster, coarse_cs=200.0, medium_cs=30, fine_cs=5,
                                    smooth_factor=1.1, bottom_model_offset=50,
                                      geol_layer_division=geol_layer_division, plot=False)

    
    
    # define node property flow package
    k_qd = 1e-7
    k_qbg = 1e-5
    k_qbo2 = 1e-7
    k_qbo1 = 1e-8
    k_rock = 1e-8

    kv_qd = 0.1  # Example horizontal hydraulic conductivity for Qd
    kv_qbg = 0.1  # Example vertical hydraulic conductivity for Qbg
    kv_qbo2 = 0.1  # Example vertical hydraulic conductivity for Qbo2
    kv_qbo1 = 0.1  # Example vertical hydraulic conductivity for Qbo
    kv_rock = 0.1  # Example vertical hydraulic conductivity for rock

    kh = np.array([
        k_qd, k_qbg, k_qbo2, k_qbo1, k_rock
    ])  # Horizontal hydraulic conductivity for each layer

    kv = np.array([
        kv_qd, kv_qbg, kv_qbo2, kv_qbo1, kv_rock
    ])  # Vertical hydraulic conductivity for each layer

    
    npf = create_npf_package(sim, kh=kh, kv=kv, geol_layer_division=geol_layer_division,
                             plot=False)

    #define the initial conditions
    start=np.empty((disv.nlay.array, disv.ncpl.array), dtype=float)
    start [:] = gwf.dis.top.array
    ic = flopy.mf6.ModflowGwfic(gwf, strt=start, pname="ic")

    #define the boundary conditions

    #define the RIVER package
    riv_chd = create_river_chd_package(sim, wd_shp)

    #define the ghb package
    ghb = create_ghb_package_lat_inflow(sim, wd_shp)

    #define the drain package for creeks

    #define the drn package for the gallery

    #define the drn package for horizontal drains

    


    return sim


if __name__ == '__main__':
    ws = os.path.join('..','03_models', 'base')
    sim = create_mf6_model(ws)
