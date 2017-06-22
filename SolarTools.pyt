import sys, os, time
import arcpy


class Toolbox(object):

    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "SolarTools"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Extract_Lidar_Building_Footprints, Build_raster_from_las, Calibrate_solar_potential]


class Extract_Lidar_Building_Footprints(object):

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Building Footprints from Lidar"
        self.description = "This tool uses classified airborne lidar las files to estimate building footprints. It " \
                           "has options to also extract trees, output a shapefile (if file size < 2Gb) and compare " \
                           "lidar-based output building footprints to an independant dataset covering the same area "
        self.canRunInBackground = False


    def getParameterInfo(self):
        """Define parameter definitions"""
        # Input Features parameter
        in_folder = arcpy.Parameter(
            displayName="Input Folder",
            name="in_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        out_csrs = arcpy.Parameter(
            displayName="Las files Coordinate System",
            name="out_csrs",
            datatype="GPCoordinateSystem",
            parameterType="Required",
            direction="Input")

        out_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        out_geodatabase = arcpy.Parameter(
            displayName="Output Geodatabase Name",
            name="out_geodatabase",
            datatype="GPString",
            parameterType="Required",
            direction="Output")

        ref_feature = arcpy.Parameter(
            displayName="Reference Feature Class for Comparison",
            name="ref_feature",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input")
        ref_feature.filter.list = ["Polygon"]

        trees_rast = arcpy.Parameter(
            displayName="Generate Trees Raster Mask",
            name="trees_rast",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        trees_rast.value = False

        out_shp = arcpy.Parameter(
            displayName="Output Shapefile",
            name="out_shp",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Output")
        out_shp.value = False

        res = arcpy.Parameter(
            displayName="Resolution used when converting lidar point cloud to raster (in squared units of coordinate system)",
            name="res",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        res.value = 0.5

        hole_area = arcpy.Parameter(
            displayName="Minimum hole area to be removed in output polygons",
            name="hole_area",
            datatype="GPArealUnit",
            parameterType="Optional",
            direction="Input")
        hole_area.value = "6 SquareMeters"

        hole_percentTotArea = arcpy.Parameter(
            displayName="Percentage of minimum area to be removed in output polygons (in percent of total building area)",
            name="hole_percentTotArea",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        hole_percentTotArea.value = 1

        simplify_tol = arcpy.Parameter(
            displayName="Tolerance used in polygon simplification step",
            name="simplify_tol",
            datatype="GPLinearUnit",
            parameterType="Optional",
            direction="Input")
        simplify_tol.value = "0.8 Meters"

        parameters = [in_folder, out_csrs, out_folder, out_geodatabase, res, hole_area, hole_percentTotArea, simplify_tol, trees_rast, ref_feature, out_shp]

        return parameters


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True


    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return


    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return


    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Get parameters
        in_folder = parameters[0].valueAsText
        out_csrs = parameters[1].valueAsText
        out_folder = parameters[2].valueAsText
        out_geodatabase = parameters[3].valueAsText
        res = parameters[4].valueAsText
        hole_area = parameters[5].valueAsText
        hole_percentTotArea = parameters[6].valueAsText
        simplify_tol = parameters[7].valueAsText
        trees_rast = parameters[8].valueAsText
        ref_feature = parameters[9].valueAsText
        ref_feature_altered = parameters[9].altered
        out_shp = parameters[10].valueAsText

        # Clean-up parameters
        res = res.replace(',','.')

        # Create output geodatabase
        arcpy.CreateFileGDB_management(out_folder, out_geodatabase)

        # Define workspace environment
        workspace = os.path.join(out_folder,out_geodatabase+'.gdb')  # workspace is the geodatabase
        arcpy.env.workspace = workspace

        # Extract building points as raster
        arcpy.AddMessage("Extracting building points as raster...")
        arcpy.management.CreateLasDataset(in_folder,os.path.join(out_folder,'montreal_las.lasd'), spatial_reference=out_csrs)
        arcpy.management.MakeLasDatasetLayer(os.path.join(out_folder,'montreal_las.lasd'), 'buildings_lyr', class_code=6)  ## extract building class (6), all returns
        arcpy.management.LasPointStatsAsRaster('buildings_lyr', 'buildings_raster', 'PREDOMINANT_CLASS', 'CELLSIZE', float(res))

        if trees_rast=='true':
            # Extract tree points as raster
            arcpy.AddMessage("Extracting tree points as raster...")
            arcpy.management.MakeLasDatasetLayer(os.path.join(out_folder,'montreal_las.lasd'), 'veg_lyr', class_code=[3,4,5])  ## extract vegetation classes (3,4,5), all returns
            arcpy.management.LasPointStatsAsRaster('veg_lyr', 'vegetation_raster', 'PREDOMINANT_CLASS', 'CELLSIZE', res)

        # Raster to polygon
        arcpy.AddMessage("Transforming raster to polygon...")
        arcpy.RasterToPolygon_conversion('buildings_raster', 'buildings_poly', 'NO_SIMPLIFY')

        # Regularize building footprints (arcpy.ddd.RegularizeBuildingFootprint() is not available before ArcGIS 10.4)
        arcpy.AddMessage("Regularizing building footprints...")
        ## Eliminate polygon holes smaaller than X m2 or Y% of polygon area
        arcpy.management.EliminatePolygonPart('buildings_poly', 'buildings_elim', condition='AREA_OR_PERCENT', part_area=hole_area, part_area_percent=hole_percentTotArea, part_option="CONTAINED_ONLY")
        ## Simplify polygon
        arcpy.cartography.SimplifyPolygon('buildings_elim', 'buildings_simplify', algorithm='POINT_REMOVE', tolerance=simplify_tol, minimum_area=hole_area, error_option='NO_CHECK')

        # Convert to shapefile
        if out_shp=='true':
            try:
                shp_path = os.path.join(out_folder,'buildings_simplify.shp')
                arcpy.AddMessage('Exporting to shapefile: {}'.format(os.path.join(out_folder,'buildings_simplify.shp')))
                arcpy.FeatureClassToShapefile_conversion('buildings_simplify', out_folder)
            except arcpy.ExecuteError:
                arcpy.AddMessage(arcpy.GetMessages())

        # Verify against reference data
        if ref_feature_altered:
            ## Generate Jaccard: Intersection/Union
            arcpy.AddMessage('Comparing to reference building footprint data...')
            arcpy.Intersect_analysis(['buildings_simplify',ref_feature], 'intersection')
            arcpy.Union_analysis(['buildings_simplify',ref_feature], 'union')
            arcpy.Statistics_analysis('intersection', 'stats_intersection', [["Shape_Area", "SUM"]])
            arcpy.Statistics_analysis('union', 'stats_union', [["Shape_Area", "SUM"]])
            with arcpy.da.SearchCursor('stats_union', ['SUM_Shape_Area']) as cursor_union:
                with arcpy.da.SearchCursor('stats_intersection', ['SUM_Shape_Area']) as cursor_intersect:
                    for i,u in zip(cursor_intersect,cursor_union):
                        jaccard = i[0]/u[0]
                        arcpy.AddMessage('The jaccard index (Intersection/Union) is: {}'.format(str(round(jaccard,3))))

        # Delete temporary files
        arcpy.AddMessage("Deleting temporary files...")
        arcpy.Delete_management(os.path.join(out_folder,'montreal_las.lasd'))
        arcpy.Delete_management('buildings_lyr')
        arcpy.Delete_management('veg_lyr')
        arcpy.Delete_management('stats_intersection')
        arcpy.Delete_management('stats_union')






class Build_raster_from_las(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Build Raster from Lidar"
        self.description = " "
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        # Input Features parameter
        in_folder = arcpy.Parameter(
            displayName="Input Folder",
            name="in_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        out_csrs = arcpy.Parameter(
            displayName="Las files Coordinate System",
            name="out_csrs",
            datatype="GPCoordinateSystem",
            parameterType="Required",
            direction="Input")

        res = arcpy.Parameter(
            displayName="Resolution used when converting lidar point cloud to raster (in squared units of coordinate system)",
            name="res",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        res.value = 0.5

        out_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Output Raster Name",
            name="out_raster",
            datatype="GPString",
            parameterType="Required",
            direction="Output")

        parameters = [in_folder, out_csrs, res, out_folder, out_raster]

        return parameters


    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Get parameters
        in_folder = parameters[0].valueAsText
        out_csrs = parameters[1].valueAsText
        res = parameters[2].valueAsText
        out_folder = parameters[3].valueAsText
        out_raster= parameters[4].valueAsText

        # Clean-up parameters
        res = res.replace(',','.')

        # Define workspace environment
        workspace = out_folder
        arcpy.env.workspace = workspace

        # Extract building points as raster
        arcpy.AddMessage("Extracting las points as raster...")
        arcpy.management.CreateLasDataset(in_folder,'montreal_las.lasd', spatial_reference=out_csrs)
        arcpy.AddMessage("Making layer...")
        arcpy.management.MakeLasDatasetLayer('montreal_las.lasd', 'lasLyr')
        arcpy.AddMessage("Las dataset to raster...")
        arcpy.conversion.LasDatasetToRaster('lasLyr', out_raster+'.tif', 'ELEVATION', 'BINNING MAXIMUM LINEAR', 'FLOAT', 'CELLSIZE', float(res))

        # Delete temporary files
        arcpy.AddMessage("Deleting temporary files...")
        arcpy.Delete_management(os.path.join(out_folder,'montreal_las.lasd'))
        arcpy.Delete_management('lasLyr')




class Calibrate_solar_potential(object):

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Calibrate Solar Potential"
        self.description = " "
        self.canRunInBackground = False


    def getParameterInfo(self):
        """Define parameter definitions"""
        # Input Features parameter
        slope_in = arcpy.Parameter(
            displayName="Input slope raster",
            name="slope_in",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        rast_in = arcpy.Parameter(
            displayName="Input radiation raster",
            name="rast_in",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        point_ref = arcpy.Parameter(
            displayName="Input reference point coordinates where data reference data was taken (zero-slope)",
            name="point_ref",
            datatype="GPPoint",
            parameterType="Optional",
            direction="Input")
        point_ref.value = "285406.224 5036270.935"

        rast_ref = arcpy.Parameter(
            displayName="Reference radiation raster",
            name="rast_ref",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        rad_month_hor = arcpy.Parameter(
            displayName="Monthly incident solar radiation on flat surface (MJ m-2) (jan;feb;mar;apr;mai;jun;jul;aug;sep;oct;nov;dec)",
            name="rad_month_hor",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        rad_month_hor.value = '5.958;9.8;14.281;16.515;20.12;21.705;21.711;18.512;13.593;8.835;5.231;4.605'

        rad_month_incl = arcpy.Parameter(
            displayName="Monthly incident solar radiation an inclined (30deg) surface (MJ m-2) (jan;feb;mar;apr;mai;jun;jul;aug;sep;oct;nov;dec)",
            name="rad_month_incl",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        rad_month_incl.value = '10.295;14.861;18.395;18.199;20.13;20.751;21.182;19.64;16.258;12.358;8.281;8.191'

        out_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        out_geodatabase = arcpy.Parameter(
            displayName="Output Geodatabase Name",
            name="out_geodatabase",
            datatype="GPString",
            parameterType="Required",
            direction="Output")

        parameters = [slope_in, rast_in, point_ref, rast_ref, rad_month_hor, rad_month_incl, out_folder, out_geodatabase]

        return parameters


    def execute(self, parameters, messages):
        """The source code of the tool."""

        ### Define functions
        # def get_bands(path_to_raster):
        #     """ Get a list of bands as Raster objects from a multiband raster 
        #     https://gis.stackexchange.com/a/150207"""
        #     oldws = arcpy.env.workspace #Save previous workspace
        #     #Get raster objects from band names
        #     arcpy.env.workspace = path_to_raster
        #     bands = [arcpy.sa.Raster(os.path.join(path_to_raster, b)) for b in arcpy.ListRasters()]
        #     #Restore previous workspace
        #     arcpy.env.workspace = oldws
        #     return bands
        ###

        # Get parameters
        slope_in = parameters[0].valueAsText
        rast_in = parameters[1].valueAsText
        point_ref = parameters[2].valueAsText
        rast_ref = parameters[3].valueAsText
        rad_month_hor = parameters[4].valueAsText
        rad_month_incl = parameters[5].valueAsText
        out_folder = parameters[6].valueAsText
        out_geodatabase = parameters[7].valueAsText

        # Create output geodatabase
        arcpy.CreateFileGDB_management(out_folder, out_geodatabase)

        # Define workspace environment
        workspace = os.path.join(out_folder,out_geodatabase+'.gdb')  # workspace is the geodatabase
        arcpy.env.workspace = workspace


        # Get lists from monthly strings
        rad_month_hor, rad_month_incl = rad_month_hor.split(";"), rad_month_incl.split(";")

        # Convert to floats and change units from MJ m-2 to WH m-2
        hor, incl = [], []
        for i,j in zip(rad_month_hor,rad_month_incl):
            hor.append(float(i)*277.777777778*365/12)
            incl.append(float(j)*277.777777778*365/12)

        # Get values of reference raster bands (months) at reference point
        arcpy.AddMessage("Obtaining value of reference raster at point...")
        cellValues = arcpy.GetCellValue_management(rast_ref, point_ref)
        values = cellValues.getOutput(0)
        values = [float(x) for x in values.split('\\n')]

        # For every month,  calculate offset in hor rad between data and calculated (offset defined as fraction of reference)
        arcpy.AddMessage("Calculating monthly offset between reference and calculated horizontal radiation...")
        hor_offset = []
        for i in range(0,12):
            hor_offset.append(values[i]/hor[i])
        arcpy.AddMessage("Offset fractions: {}".format(hor_offset))


        # For every month, apply offset to input raster
        arcpy.AddMessage("Applying monthly offset to input raster and calculating annual sum...")
        #bands = get_bands(rast_in)
        oldws = arcpy.env.workspace #Save previous workspace
        #Get raster objects from band names
        arcpy.env.workspace = rast_in
        bands = [arcpy.sa.Raster(os.path.join(rast_in, b)) for b in arcpy.ListRasters()]
        #Restore previous workspace
        arcpy.env.workspace = oldws
        months = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
        monthly_values, monthly_incl_addedFrac = [], []
        for month, i in zip(months, range(0,12)):
            arcpy.AddMessage("Correcting {} values...".format(month))
            monthly_values.append(arcpy.sa.Float(bands[i]) * hor_offset[i])
            monthly_incl_addedFrac.append(incl[i]/hor[i] * hor_offset[i])

        # Calculate annual sum
        arcpy.AddMessage("Calculating monthly and annual sums...")
        summed = arcpy.sa.Float(bands[0]) * 0
        summed_wIncl = arcpy.sa.Float(bands[0]) * 0
        for i in range(0,12):
            arcpy.AddMessage("Processing month {}".format(i+1))
            monthly_values[i].save("month{}".format(i+1))
            monthly_incl = arcpy.sa.Con(slope_in < 2, monthly_values[i] * monthly_incl_addedFrac[i] , monthly_values[i])
            monthly_incl.save("month{}_wIncl".format(i+1))
            summed += monthly_values[i]
            summed_wIncl += monthly_incl
        arcpy.AddMessage("Saving annual sums...")
        summed.save("calibrated_solar_radiation_annualSum")
        summed_wIncl.save("calibrated_solar_radiation_annualSum_wIncl")


