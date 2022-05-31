
import pycircos
import matplotlib.pyplot as plt
import pandas as pd

Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle

networkR = [680,730]
lineR1 = [600,630]
lineR2 = [635,665]
nodeR = [560, 580]
linkR = 575

import collections
color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
               "acen":"#D82322"}

colorlist = ["#ff8a80","#ff80ab","#ea80fc","#b388ff","#8c9eff","#84ffff","#b9f6ca","#ccff90",
"#ffff8d","#ffd180","#ff9e80","#bcaaa4","#ff5252","#e040fb","#7c4dff","#536dfe","#448aff","#18ffff","#64ffda","#69f0ae","#b2ff59","#eeff41","#ffff00","#ffd740","#ffab40","#ff6e40","#a1887f","#e0e0e0","#90a4ae"]

colorlist = ['#ff7f00', '#fdbf6f', '#e31a1c', '#fb9a99', '#33a02c', '#b2df8a', '#1f78b4', '#a6cee3',
"#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00"]

colorlist = ["#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3",
'#b3b3b3', '#e5c494', '#ffd92f', '#a6d854', '#e78ac3', '#8da0cb', '#fc8d62', '#66c2a5']

colorlist = ['#7C1D6F', '#DC3977', '#F0746E', '#FCDE9C', '#b3b3b3', '#7CCBA2', '#089099', '#045275',
"#045275","#089099","#7CCBA2","#b3b3b3","#FCDE9C","#F0746E","#DC3977","#7C1D6F"]

colorlist = ['#7C1D6F', '#b3b3b3','#DC3977', '#089099','#F0746E', '#045275', '#FCDE9C',  '#7CCBA2', 
'#7CCBA2', '#FCDE9C', '#045275', '#F0746E', '#089099', '#DC3977', '#b3b3b3', '#7C1D6F']

# colorlist = ["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf",
# '#f781bf', '#a65628', '#ffff33', '#ff7f00', '#984ea3', '#4daf4a', '#377eb8', '#e41a1c']



#colorlist = ["#ff8a80","#ff80ab","#ea80fc","#b388ff","#8c9eff","#82b1ff","#84ffff","#a7ffeb","#b9f6ca","#ccff90","#f4ff81","#ffff8d","#ffe57f","#ffd180","#ff9e80","#bcaaa4","#eeeeee","#b0bec5","#ff5252","#ff4081","#e040fb","#7c4dff","#536dfe","#448aff","#18ffff","#64ffda","#69f0ae","#b2ff59","#eeff41","#ffff00","#ffd740","#ffab40","#ff6e40","#a1887f","#e0e0e0","#90a4ae"]

ageList = [[18,35],[35,50],[50,65],[65,88]]

freqList = [4,8,12,20]
#freqList = [4]
for age in ageList:
    for freq in freqList:

        #Set chromosomes
        circle = Gcircle() 
        df = pd.read_csv("circos_data/NodeNetwork_general.csv")

        for i in range(len(df)):
            name   = df.network[i] 
            length = df.end[i]
            arc    = Garc(arc_id=name, size=length, interspace=2, raxis_range=networkR, labelposition=70, labelsize= 11,label_visible=True,facecolor=colorlist[i])
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
                            raxis_range=nodeR, facecolor=circle.garc_dict[key].facecolor, markersize=8,spine=True) 


        values_all   = [] 
        #arcdata_dict = collections.defaultdict(dict)
        
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
            #circle.chord_plot(source, destination, facecolor=circle.garc_dict[name1].facecolor)
            circle.chord_plot(source, destination, facecolor="#808080")


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
        titleStr = f'{freq}Hz_{age[0]}to{age[1]}'
        fig.suptitle(titleStr, fontsize=20)
        circle.figure.suptitle(titleStr, fontsize=20)
        picName = f'circos_pics/linkplot_{freq}Hz_{age[0]}to{age[1]}.pdf'
        circle.figure.savefig(picName)