


import pycircos
import matplotlib.pyplot as plt
Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle

import collections
# color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
#                "acen":"#D82322"}
color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
               "acen":"#D82322"}

#Set chromosomes
circle = Gcircle() 
#with open("sample_data/example_data_chromosome_general.csv") as f:
with open("circos_data/NodeNetwork_general.csv") as f:
    f.readline()
    for line in f:
        line   = line.rstrip().split(",") 
        name   = line[0]
        length = int(line[-1]) 
        arc    = Garc(arc_id=name, size=length, interspace=2, raxis_range=(950,1000), labelposition=60, label_visible=True)
        circle.add_garc(arc) 
circle.set_garcs(0,359) 

# #cytoband
import collections
color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
               "acen":"#D82322"}

# arcdata_dict = collections.defaultdict(dict)
# with open("sample_data/example_data_chromosome_cytoband.csv") as f:
#     f.readline()
#     for line in f:
#         line  = line.rstrip().split(",")
#         name  = line[0]     
#         start = int(line[1])-1 
#         width = int(line[2])-(int(line[1])-1) 
#         if name not in arcdata_dict:
#             arcdata_dict[name]["positions"] = []
#             arcdata_dict[name]["widths"]    = [] 
#             arcdata_dict[name]["colors"]    = [] 
#         arcdata_dict[name]["positions"].append(start) 
#         arcdata_dict[name]["widths"].append(width)
#         arcdata_dict[name]["colors"].append(color_dict[line[-1]])

# for key in arcdata_dict:
#     circle.barplot(key, data=[1]*len(arcdata_dict[key]["positions"]), positions=arcdata_dict[key]["positions"], 
#                    width=arcdata_dict[key]["widths"], raxis_range=[950,1000], facecolor=arcdata_dict[key]["colors"])    
    

values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("circos_data/NodeConnUni_linkplot_Theta_27to30.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name1  = line[0]     
        start1 = int(line[1])-1
        end1   = int(line[2])
        name2  = line[3]     
        start2 = int(line[4])-1
        end2   = int(line[5])
        source = (name1, start1, end1, 830)
        destination = (name2, start2, end2, 830)
        circle.chord_plot(source, destination, facecolor=circle.garc_dict[name1].facecolor)






circle.figure

circle.figure.savefig("tutotial.pdf")