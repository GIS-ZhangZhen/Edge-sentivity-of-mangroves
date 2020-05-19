# -*- coding: utf-8 -*-
import arcpy as ap
import numpy as np
import scipy
import scipy.stats
from collections import deque 
import matplotlib.pyplot as plt
import math
import datetime
import sys

ap.env.overwriteOutput = True
if ap.CheckExtension("spatial") == "Available":
    ap.CheckOutExtension("spatial")
else:
    print("License error")


def FindRadius(polygon,center):
    """
    返回一个多边形的最大半径
    """
    X_list,Y_list=[],[]
    for row in ap.da.SearchCursor(polygon,["SHAPE@"]):
        for part in row[0]:
            for point in part:
                if point:
                    X_list.append("{0},{1}".format(point.X,point.Y))
                    Y_list.append("{0},{1}".format(point.Y,point.X))
        a,b,c,d=max(X_list).split(","), min(X_list).split(","), max(Y_list).split(","), min(Y_list).split(",")
        max_X, min_X=[float(a[0]),float(a[1])], [float(b[0]),float(b[1])]
        max_Y, min_Y=[float(c[1]),float(c[0])], [float(d[1]),float(d[0])]
        
    distance=[]
    x,y=center[0],center[1]
    distance.append(((max_X[0]-x)**2.0+(max_X[1]-y)**2.0)**(1/2.0))
    distance.append(((min_X[0]-x)**2.0+(min_X[1]-y)**2.0)**(1/2.0))
    distance.append(((max_Y[0]-x)**2.0+(max_Y[1]-y)**2.0)**(1/2.0))
    distance.append(((min_Y[0]-x)**2.0+(min_Y[1]-y)**2.0)**(1/2.0))
    return max(distance)

def CreateLine(center,point):
    points = [center,point]
    array = ap.Array()
    point_set = ap.Point()
    for i in range(2):
        point_set.ID= i
        point_set.X = points[i][0]
        point_set.Y = points[i][1]
        array.add(point_set)
    polyline=ap.Polyline(array)
    return polyline


def ShapeIndex(polygon,id,n=32):
    ap.env.workspace='F:/DATA/TemporalFile/'
    point="point{0}.shp".format(id)
    point_result=ap.FeatureToPoint_management(polygon,point, "INSIDE")
    center=[row[0] for row in ap.da.SearchCursor(point_result,['SHAPE@XY'])]
    center=[center[0][0],center[0][1]]
    ap.Delete_management(point)
    r = FindRadius(polygon,center)
    
    alpha = 360.0/n
    fc="NewLine{0}.shp".format(id)
    line_result=ap.CreateFeatureclass_management(ap.env.workspace,
                                     fc,"Polyline")
    cursor = ap.da.InsertCursor(line_result,["SHAPE@"])
    for num in range(n):
        x = center[0] + math.sin(math.radians(alpha*num))*r
        y = center[1] + math.cos(math.radians(alpha*num))*r
        point = [x,y]
        line = CreateLine(center,point)
        cursor.insertRow([line])
    del cursor
    ap.AddField_management(fc,'ID1','SHORT')
    ap.CalculateField_management(fc,'ID1','!FID!',"PYTHON_9.3")
    clip = "clip{0}.shp".format(id)
    clip_result=ap.Clip_analysis(line_result,polygon,clip)
    length = [row[1] for row in sorted(ap.da.SearchCursor(clip_result,['ID1','SHAPE@LENGTH']))]
    ap.Delete_management(fc)
    ap.Delete_management(clip)
    return np.array(length)

def normalization(data):
    range = np.max(data) - np.min(data)
    if range != 0:
        return (data - np.min(data)) / range
    else :
        return np.ones(len(data))

def Shed(List, n):
    List = deque([x for x in List])
    for i in range(n):
        b = List.popleft()
        List.append(b)
    return List

def FillPolygon(infc,n):
    data=[]
    spatial_ref = ap.Describe(infc).spatialReference
    for row in ap.da.SearchCursor(infc, ["OID@", "SHAPE@"]):
        data.append([])
        for part in row[1]:
            for pnt in part:
                if not pnt :break
                data[-1].append([pnt.X,pnt.Y])
        features = []

    for feature in data:
        features.append(
            ap.Polygon(
                ap.Array([ap.Point(*coords) for coords in feature]),spatial_ref))
    outfc = 'F:/DATA/TemporalFile/outfc{0}.shp'.format(n)
    fill_result=ap.CopyFeatures_management(features, outfc )
    return fill_result

def CursorSBC(row):
    outfeature = "{0}.shp".format(row[0])
    ToPolygon_result=ap.FeatureToPolygon_management(row[1],outfeature)
    Fill_polygon = FillPolygon(ToPolygon_result,row[0])
    q = ShapeIndex(Fill_polygon,row[0])
    q=normalization(np.around(q,3))
    
    ap.Delete_management(Fill_polygon)
    ap.Delete_management(outfeature)
    return q

#---------------------------------得到标准图像的信号序列--------------------------------------------
start_time = datetime.datetime.now()
print "Start Time:",start_time
SBC_list = []
model = 'Z:/home/项目1_红树林与互花米草关系/Data/test/标准形状.shp'
for row in sorted(ap.da.SearchCursor(model, ["FID_match","Class_Name","SHAPE@"])):
    polygon=row[2]
    outfeature = "F:/DATA/TemporalFile/{0}.shp".format(row[0])
    geometry=ap.FeatureToPolygon_management(polygon,outfeature)
    p = ShapeIndex(geometry,row[0])
    ap.Delete_management(outfeature)
    name = row[1]
    SBC_list.append(name)
    globals()[name] = normalization(p)

#---------------------------------班块形状识别--------------------------------------------
mangrove = 'Z:/home/项目1_红树林与互花米草关系/Data/mangrove1501/mangrove1501.shp'
for row in ap.da.SearchCursor(mangrove, ["OID@","SHAPE@","ID"]):
    q=CursorSBC(row)
    js_list1,shape_list1=[],[]
    for i in range(len(q)):
        a = Shed(q,i)
        js_list2 = []
        for shape in SBC_list:
            p=globals()[shape]
            M=(p+a)/2
            js2=0.5*scipy.stats.entropy(p, M)+0.5*scipy.stats.entropy(a, M)
            js_list2.append(js2)
        js1 = min(js_list2)
        index1 = js_list2.index(js1)
        js_list1.append(js1)
        shape_list1.append(SBC_list[index1])
    js = min(js_list1)
    index = js_list1.index(js)
    print(row[2],js,shape_list1[index])
    txt1 = "mangrove_1501_test.txt"
    f=open(txt1,'a') 
    f.writelines('{0},{1},{2}\n'.format(row[2],js,shape_list1[index]))
end_time = datetime.datetime.now()
print "Succeeded at",end_time
print "Elapsed Time:",end_time-start_time





















