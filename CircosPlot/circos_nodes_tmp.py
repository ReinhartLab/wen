
import pycircos
import matplotlib.pyplot as plt
import pandas as pd

Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle

networkR = [680,730]
lineR1 = [600,630]
lineR2 = [635,665]
nodeR = [570, 590]
linkR = 585

import collections
color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
               "acen":"#D82322"}
#age = [18,35]
# age = [35,50]
# age = [50,65]
age = [65,88]
freqList = [4,8,12,20]

for freq in freqList:

    #Set chromosomes
    circle = Gcircle() 
    df = pd.read_csv("circos_data/NodeNetwork_general.csv")

    for i in range(len(df)):
        name   = df.network[i] 
        length = df.end[i]
        arc    = Garc(arc_id=name, size=length, interspace=2, raxis_range=networkR, labelposition=70, labelsize= 11,label_visible=True)
        circle.add_garc(arc) 
        circle.set_garcs(1,359) 


    #scatter plot
    # values_all   = [] 
    arcdata_dict = collections.defaultdict(dict)

    df = pd.read_csv("circos_data/NodeNetwork_dots.csv")

    for i in range(len(df)):
        name   = df.network[i] 
        mid = df.center[i]
        value = df.value[i]
        if name not in arcdata_dict:
                arcdata_dict[name]["positions"] = []
                arcdata_dict[name]["values"] = []
        arcdata_dict[name]["positions"].append(mid) 
        arcdata_dict[name]["values"].append(value)

    for key in arcdata_dict:
        circle.scatterplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                        raxis_range=nodeR, facecolor="orangered", markersize=8,spine=True) 


    values_all   = [] 
    #arcdata_dict = collections.defaultdict(dict)
    
    #linkcsv = f'circos_data/NodeConn_linkplot_{freq}Hz_{age[0]}to{age[1]}_hipp.csv'
    linkcsv = f'circos_data/NodeConn_linkplot_{freq}Hz_{age[0]}to{age[1]}_hipp.csv'
    df = pd.read_csv(linkcsv)

    for i in range(len(df)):    
    
        name1  = df.network[i]    
        start1 = df.start1[i]
        end1   = df.end1[i]
        name2  = df.network2[i]    
        start2 = df.start2[i]
        end2   = df.end2[i]
        source = (name1, start1, end1, linkR)
        destination = (name2, start2, end2, linkR)
        circle.chord_plot(source, destination, facecolor=circle.garc_dict[name1].facecolor)


    #line plot
    values_all   = [] 
    values_all2 = []
    arcdata_dict = collections.defaultdict(dict)
    linecsv = f'circos_data/NodeConn_line_{freq}Hz_{age[0]}to{age[1]}_hipp.csv'

    df = pd.read_csv(linecsv)
    for i in range(len(df)):    
        
            name  = df.label[i]    
            start = df.start[i]-1
            end   = df.end[i]
            mid   = (start+end)/2
            value   = df.value1[i]
            values_all.append(value) 
            value2 = df.value2[i]
            values_all2.append(value2) 
            if name not in arcdata_dict:
                arcdata_dict[name]["positions"] = []
                arcdata_dict[name]["values"]    = [] 
                arcdata_dict[name]["values2"]    = [] 
            arcdata_dict[name]["positions"].append(mid) 
            arcdata_dict[name]["values"].append(value)
            arcdata_dict[name]["values2"].append(value2)

    vmin, vmax = min(values_all), max(values_all) 
    for key in arcdata_dict:
        circle.lineplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                        rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)], raxis_range=lineR1, linecolor="royalblue", spine=False) 

    vmin, vmax = min(values_all2), max(values_all2) 
    for key in arcdata_dict:
        circle.lineplot(key, data=arcdata_dict[key]["values2"], positions=arcdata_dict[key]["positions"], 
                        rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)], raxis_range=lineR2, linecolor="royalblue", spine=False) 

    fig = plt.figure()
    fig.suptitle('test title', fontsize=20)
    circle.figure.suptitle('Test title', fontsize=20)
    picName = f'circos_pics/linkplot_{freq}Hz_{age[0]}to{age[1]}.pdf'
    circle.figure.savefig(picName)
