# acse-ep20

from qgis.core import *
import numpy as np
import numpy.linalg as la
import os
a = 1
class Contour(object):
    def __init__(self, enclosed_file, unenclosed_file, target_enclosed_filename, target_unenclosed_filename):
        """
        enclosed file is the contour vector file with manual boundary
        unencloed_file is contour vector file with only geometries
        """
        self.enclosed_file = enclosed_file
        self.unenclosed_file = unenclosed_file
        self.enclosed_layer = self.get_layer_information(enclosed_file)
        self.unenclosed_layer = self.get_layer_information(unenclosed_file)
        self.target_enclosed_filename = target_enclosed_filename
        self.target_unenclosed_filename = target_unenclosed_filename
        self.enclosed_features_list = []
        self.unenclosed_features_list = []
        self.classify_features_by_enclosed_or_not(self.unenclosed_layer)

    def get_layer_information(self, file_name):
        """
        get the temp layer for vector manipulation
        return as temp layer
        """
        vector_layer = QgsVectorLayer(file_name, 'Shapefile', 'ogr')
        features = [feat for feat in vector_layer.getFeatures()]
        mem_layer = QgsVectorLayer('LineString?crs=epsg:4326', "duplicated","memory")
        mem_layer_data = mem_layer.dataProvider()
        attr = vector_layer.dataProvider().fields().toList()
        mem_layer_data.addAttributes(attr)
        mem_layer.updateFields()
        mem_layer_data.addFeatures(features)
        return mem_layer
    
    def deepcopy_layer(self, layer):
        """
        return a deepcopy of a layer
        """
        features = [feat for feat in layer.getFeatures()]
        mem_layer = QgsVectorLayer('LineString?crs=epsg:4326', "duplicated","memory")
        mem_layer_data = mem_layer.dataProvider()
        attr = layer.dataProvider().fields().toList()
        mem_layer_data.addAttributes(attr)
        mem_layer.updateFields()
        mem_layer_data.addFeatures(features)
        return mem_layer

    def split_feature_by_extent(self, line, extent):
        """
        giving and extent, and a line, split them by indices
        """
        xiMin, xiMax, etaMin, etaMax = extent
        # two flag to determin whether two index found
        search_first = False
        counter = 0
        for index, point in enumerate(line):
            x = point[0]
            y = point[1]
            if xiMin <= x and x <= xiMax and etaMin <= y and y <= etaMax and not search_first:
                end_first = index - 1
                search_first = True
            if search_first and not(xiMin <= x and x <= xiMax and etaMin <= y and y <= etaMax):
                end_sec = index - 1
                return (end_first, end_sec)
            if search_first:
                counter += 1

    def smooth_a_part(self, points, points_per_group=20, overlapping_points=19, num_modes=2):
        """
        this function provides the detail of the PCA smoothing algorithm
        """
        num_points = len(points)
        num_groups = (num_points - overlapping_points) // (points_per_group-overlapping_points)
        x_average = np.zeros(num_groups)
        y_average = np.zeros(num_groups)
        Sx = np.zeros((points_per_group, num_groups))
        Sy = np.zeros((points_per_group, num_groups))
        # get the matrix ready from overlapping
        for i in range(num_groups):
            # this one is ok
            start_index = i * (points_per_group-overlapping_points)
            for j in range(points_per_group):
                Sx[j, i] = points[start_index + j][0]
                Sy[j, i] = points[start_index + j][1]
                x_average[i] += points[start_index + j][0]
                y_average[i] += points[start_index + j][1]
            x_average[i] /= points_per_group
            y_average[i] /= points_per_group
            for j in range(points_per_group):
                Sx[j, i] = points[start_index + j][0] - x_average[i]
                Sy[j, i] = points[start_index + j][1] - y_average[i]
    # the construction of matrix is perfectly OK
        Cx = 1/num_groups * Sx @ Sx.T
        Cy = 1/num_groups * Sy @ Sy.T
        eigenvalue_x, Ax= la.eig(Cx)

        eigenvalue_y, Ay= la.eig(Cy)
        PCA_modes_x = Sx.T @ Ax

        PCA_modes_y = Sy.T @ Ay

        # parameter that could be adjusted
        reconstruct_Sx =  Ax[:,:num_modes]@PCA_modes_x[:,:num_modes].T
        reconstruct_Sy = Ay[:,:num_modes]@PCA_modes_y[:,:num_modes].T

        new_points_x = np.zeros(num_points)
        new_points_y = np.zeros(num_points)
        # this one is used to count the number of appearance
        new_points_counter = np.zeros(num_points)
        new_points = []
        # first group length
        # get the sum and count of xs and ys
        for i in range(num_groups):
            start_index = i * (points_per_group-overlapping_points)
            for j in range(points_per_group):
                new_points_x[start_index + j]
                reconstruct_Sx[j, i]
                new_points_x[start_index + j] += reconstruct_Sx[j, i] + x_average[i]
                new_points_y[start_index + j] += reconstruct_Sy[j, i] + y_average[i]
                new_points_counter[start_index + j] += 1
        # get the average
        for i in range(num_points):
            new_points_x[i] /= new_points_counter[i]
            new_points_y[i] /= new_points_counter[i]
            new_points.append(QgsPointXY(new_points_x[i], new_points_y[i]))
        return new_points

    def check_unenclosed_feature_numbers(self):

        if len(self.unenclosed_features_list) > 1:
            raise("more than one unenclosed shoreline(feature) is included")


    def check_whether_feature_all_in_extent(self, feature, extent):
        """
        check whetehr all points in a feature in an extent
        """
        xiMin, xiMax, etaMin, etaMax = extent 
        area_points = []
        geom = feature.geometry()
        point_counter = 0
        if geom.wkbType() == QgsWkbTypes.MultiLineString:
            temp_lines = geom.asMultiPolyline()
            for t in temp_lines:
                t = geom.fromPolylineXY(t)
                pointList = t.asPolyline()
                for point in pointList:
                    x = point.x()
                    y = point.y()
                    point_counter += 1
                    if xiMin <= x and x <= xiMax and etaMin <= y and y <= etaMax:
                        area_points.append(point)
            if point_counter == len(area_points):
                return True
            else:
                return False
    
    def check_whether_a_line_continous_in_one_extent(self, manipulated_features, extent):
        """given an extent and features,
        check whether all manipulated features are continuous in the extent
        """
        xiMin, xiMax, etaMin, etaMax = extent 
        for feature in manipulated_features:
            area_points = []
            geom = feature.geometry()
            point_counter = 0
            if geom.wkbType() == QgsWkbTypes.MultiLineString:
                temp_lines = geom.asMultiPolyline()
                for t in temp_lines:
                    t = geom.fromPolylineXY(t)
                    pointList = t.asPolyline()
                    for point in pointList:
                        x = point.x()
                        y = point.y()
                        if xiMin <= x and x <= xiMax and etaMin <= y and y <= etaMax:
                            area_points.append(point)
                counter = 0
                counting_flag = False
                end_flag = False
                if len(area_points):
                    for t in temp_lines:
                        t = geom.fromPolylineXY(t)
                        pointList = t.asPolyline()
                        for point in pointList:
                            if point.compare(area_points[0]) and not counting_flag:
                                counting_flag = True
                            if counting_flag:
                                if point.compare(area_points[counter]):
                                    counter += 1
                                    if counter == len(area_points):
                                        break
                                else:
                                    end_flag = True
                                    break
                        if end_flag:
                            break
                    if not counter == len(area_points):
                        return False
        return True

    def saving_a_layer(self, vector_layer, shapefile_name):
        """
        save a layer to a file
        """
        save_options = QgsVectorFileWriter.SaveVectorOptions()
        save_options.driverName = "ESRI Shapefile"
        save_options.fileEncoding = "UTF-8"
        transform_context = QgsProject.instance().transformContext()
        error = QgsVectorFileWriter.writeAsVectorFormatV2(vector_layer,
                                                        shapefile_name,
                                                        transform_context,
                                                        save_options)

    def classify_features_by_enclosed_or_not(self, vector_layer):
        """
        enclosed  features need to be adjust to
        """
        manipulated_features = vector_layer.getFeatures()
        
        for feature in manipulated_features:
            geom = feature.geometry()
            
            if geom.wkbType() == QgsWkbTypes.MultiLineString:
                    temp_lines = geom.asMultiPolyline()
                    first = temp_lines[0][0]
                    last = temp_lines[-1][-1]
                    if first.compare(last):
                        self.enclosed_features_list.append(feature.id())
                    else:
                        self.unenclosed_features_list.append(feature.id())

    def delete_unnecessary_feature_for_unenclose(self, vector_layer, uplimit_number, extent):
        """
        delete the small geometries for
        """
        temp_layer = self.deepcopy_layer(vector_layer)
        manipulated_features = temp_layer.getFeatures()
        counter = 0
        short_features = []
        for index, feature in enumerate(manipulated_features):
            geom = feature.geometry()
            point_counter = 0
            if not self.check_whether_feature_all_in_extent(feature, extent):
                if geom.wkbType() == QgsWkbTypes.MultiLineString:
                    temp_lines = geom.asMultiPolyline()
                    for t in temp_lines:
                        t = geom.fromPolylineXY(t)
                        pointList = t.asPolyline()
                        point_counter += len(pointList)
                    if point_counter < uplimit_number:
                        counter += 1
                        short_features.append(feature.id())
        caps = temp_layer.dataProvider().capabilities()
        # Check if a particular capability is supported:
        if caps & QgsVectorDataProvider.DeleteFeatures:
            res = temp_layer.dataProvider().deleteFeatures(short_features)
        self.saving_a_layer(temp_layer, self.target_unenclosed_filename)
        return short_features
    
    def check_connect_manual_boundary_to_natural_boundary(self):
        """
        this function connects the manual boundary to natural boundary
        0 first first
        1 represents natural first manual last
        2 represents natural last and manual last
        3 represents natural last and manual first
        """
        feature_ids = []
        for feature in self.unenclosed_layer.getFeatures():
            feature_ids.append(feature.id())
        for feature in self.enclosed_layer.getFeatures():
            if feature.id() not in feature_ids:
                # usually, manual boundary is a connected feature
                manual_boundary_id = feature.id()
                break
        request = QgsFeatureRequest()
        request.setFilterFid(manual_boundary_id)
        features = self.enclosed_layer.getFeatures(request)
        for feature in features:
            geom = feature.geometry()
            if geom.wkbType() == QgsWkbTypes.MultiLineString:
                temp_lines = geom.asMultiPolyline()
                first = temp_lines[0][0]
                last = temp_lines[-1][-1]
        
        # now find the corresponding point
        request = QgsFeatureRequest()
        request.setFilterFids(self.unenclosed_features_list)
        features = self.enclosed_layer.getFeatures(request)
        for feature in features:
            natural_boundary = feature
            geom = feature.geometry()
            if geom.wkbType() == QgsWkbTypes.MultiLineString:
                temp_lines = geom.asMultiPolyline()
                first_natural = temp_lines[0][0]
                last_natural = temp_lines[-1][-1]
        # 
        if first_natural.compare(first):
            return 0, 0, manual_boundary_id
        elif first_natural.compare(last):
            return 0, -1, manual_boundary_id
        else:
            raise("Boundary is not closed!!")
            
    def connect_the_boundary(self, index_natural, index_manual, manual_id):
        """
        this function connects the boundary
        """
        request = QgsFeatureRequest()
        request.setFilterFids(self.unenclosed_features_list)
        features = self.enclosed_layer.getFeatures(request)
        for feature in features:
            natural_boundary = feature
            geom = feature.geometry()
            if geom.wkbType() == QgsWkbTypes.MultiLineString:
                print(1)
                temp_lines = geom.asMultiPolyline()
                first_natural = temp_lines[0][0]
                last_natural = temp_lines[-1][-1]
            elif geom.wkbType()==QgsWkbTypes.LineString:
                temp_line = geom.asPolyline()
                first_natural = temp_line[0]
                last_natural = temp_line[-1]
        request = QgsFeatureRequest()
        request.setFilterFid(manual_id)
        features = self.enclosed_layer.getFeatures(request)
        for feature in features:
            points = []
            natural_boundary = feature
            geom = feature.geometry()
            if geom.wkbType() == QgsWkbTypes.MultiLineString:
                temp_lines = geom.asMultiPolyline()
                for t in temp_lines:
                    points += t
        if index_natural == 0 and index_manual == 0:
            points[0] = first_natural
            points[-1] = last_natural
        if index_natural == 0 and index_manual == -1:
            points[0] = last_natural
            points[-1] = first_natural
        self.enclosed_layer.dataProvider().changeGeometryValues({manual_id: QgsGeometry.fromPolylineXY(points)})
        

        
    def simplify_the_feature(self, extent, points_per_group=20, overlapping_points=19, num_modes=2):
        index_natural, index_manual, manual_id = self.check_connect_manual_boundary_to_natural_boundary()
        xiMin, xiMax, etaMin, etaMax = extent
        unenclosed_layer = self.unenclosed_layer
        manipulated_features = unenclosed_layer.getFeatures()
        if not self.check_whether_a_line_continous_in_one_extent(manipulated_features, extent):
            raise("the part is discountinous in the extent")
        # refresh the iterator
        manipulated_features = unenclosed_layer.getFeatures()
        for feature in manipulated_features:
            area_points = []
            geom = feature.geometry()
            point_counter = 0
            points = []
            # if the a feature all in one extent, then it will not be processed
            if not self.check_whether_feature_all_in_extent(feature, extent):
                if geom.wkbType() == QgsWkbTypes.MultiLineString:
                    temp_lines = geom.asMultiPolyline()
                    for t in temp_lines:
                        t = geom.fromPolylineXY(t)
                        pointList = t.asPolyline()
                        for point in pointList:
                            x = point.x()
                            y = point.y()
                            points.append((x, y))
                    split_indices = self.split_feature_by_extent(points, extent)
                    if split_indices:
                        line_one = points[:split_indices[0] + 1]
                        line_middle = points[split_indices[0] + 1 : split_indices[1] + 1]
                        line_end = points[split_indices[1] + 1:]
                        line_middle_transformed = []
                        for point in line_middle:
                            line_middle_transformed.append(QgsPointXY(point[0], point[1]))
                        line_one_smoothed = self.smooth_a_part(line_one, points_per_group, overlapping_points, num_modes)
                        line_end_smoothed = self.smooth_a_part(line_end, points_per_group, overlapping_points, num_modes)
                        if feature.id() in self.enclosed_features_list:
                            line_end_smoothed[-1] = line_one_smoothed[0]
                        complete_line = line_one_smoothed + line_middle_transformed + line_end_smoothed
                        smoothed_line = QgsGeometry.fromPolylineXY(complete_line)
                        unenclosed_layer.dataProvider().changeGeometryValues({feature.id(): smoothed_line})
                    else:
                        if len(points) > points_per_group * 3:
                            line_smoothed = self.smooth_a_part(points, points_per_group, overlapping_points, num_modes)
                            if feature.id() in self.enclosed_features_list:
                                line_smoothed[-1] = line_smoothed[0] 
                            smoothed_line = QgsGeometry.fromPolylineXY(line_smoothed)
                            unenclosed_layer.dataProvider().changeGeometryValues({feature.id(): smoothed_line})


        manipulated_features = unenclosed_layer.getFeatures()
        for feature in manipulated_features:
            self.enclosed_layer.dataProvider().changeGeometryValues({feature.id():feature.geometry()})

        deleted_id = self.delete_unnecessary_feature_for_unenclose(unenclosed_layer,  points_per_group * 3, extent)
        caps = self.enclosed_layer.dataProvider().capabilities()
        # Check if a particular capability is supported:
        if caps & QgsVectorDataProvider.DeleteFeatures:
            res = self.enclosed_layer.dataProvider().deleteFeatures(deleted_id)
        self.connect_the_boundary(index_natural, index_manual, manual_id)
        self.saving_a_layer(self.enclosed_layer, self.target_enclosed_filename) 

