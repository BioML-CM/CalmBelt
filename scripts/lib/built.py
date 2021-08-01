import pandas as pd
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt

import plotly.express as px
import plotly.graph_objects as go

from lib import alignment

from plotly.figure_factory import create_dendrogram
from scipy.spatial.distance import pdist

import plotly.figure_factory as ff

import itertools


#----------------------------------------------------------------
def get_tick(scale):
    genome_size = 29903
    tickvals = [scale, scale*2, scale*3, scale*4, scale*5, scale*6, scale*7, scale*8,
                scale*9, scale*10, scale*11, scale*12, scale*13, scale*14, scale*15, scale*16, 
                scale*17, scale*18, scale*19, scale*20, scale*21]
    ticktext = ['Jan-2020', 'Feb-2020', 'Mar-2020', 'Apr-2020', 'May-2020', 
                'Jun-2020', 'July-2020', 'Aug-2020', 'Sep-2020', 'Oct-2020', 'Nov-2020', 
                'Dec-2020', 'Jan-2021', 'Feb-2021', 'Mar-2021', 'Apr-2021', 'May-2021', 
                'Jun-2021', 'July-2021', 'Aug-2021', 'Sep-2021']
    return tickvals, ticktext
#----------------------------------------------------------------

def find_insertion_position(id_insertion,snps_insertion_df):
    insertion_list = []
    for sample in id_insertion:
        se_list = list(set(snps_insertion_df[snps_insertion_df['id']==sample]['start_end']))

        for se in se_list: #แยกคิดแต่ละสาย
            list_temp = sorted(list(snps_insertion_df[(snps_insertion_df['id']==sample)&(snps_insertion_df['start_end']==se)]['position']))
            list_temp+=[0]

            count_insertion = [] #จำนวนตำแหน่งที่่อยู่ติดกัน     
            pos_temp = []
            sum_count = 1
            for i in range(len(list_temp)-1):
                if list_temp[i] == list_temp[i+1]-1:
                    sum_count += 1
                else:
                    count_insertion += [sum_count]   #จำนวน
                    pos_temp += [list_temp[i]]     #ตำแหน่ง
                    sum_count = 1



            #ตำแหน่งตาม ref
            pos_insertion = [i-j for i,j in zip(pos_temp,np.cumsum(count_insertion))]
            insertion_list += [[sample,i,j,se] for i,j in zip(pos_insertion,count_insertion)]
    return insertion_list

#----------------------------------------------------------------------

def fig_dendro_built_tree(X_dendro, X_id_list, clade_drop_dendro, snps_max_df, change_protein_df):

    n = len(X_id_list)

    fig_dendro = create_dendrogram(X_dendro, orientation='right', labels=X_id_list,
                           distfun=lambda x: pdist(x,metric='hamming'))
    
    labels_list = fig_dendro.layout['yaxis']['ticktext'] #label เรียงตามลำดับจากล่างขึ้นบน
    
    
    # สร้าง heatmap
    m = np.zeros([n,n])
    vals = pdist(X_dendro,metric='hamming')

    xs,ys = np.triu_indices(n,k=1)
    m[xs,ys] = vals
    m[ys,xs] = vals
    m[ np.diag_indices(n) ] = 0
    m = m*29903  #จำนวน dna
    m = np.round(m, 0).astype(int)        
    
    
    new_index = [X_id_list.index(k) for k in list(labels_list)] 
    
    m = m[new_index] #row=yaxis
    m = m[:,new_index[::-1]]  #column=xaxis
    
    

    fig_heatmap = ff.create_annotated_heatmap(m,x=list(labels_list[::-1]),y=list(labels_list), colorscale='cividis', hoverinfo='none')
    
    
    heatmap_df = pd.DataFrame(m[::-1],index=list(labels_list[::-1]), columns=list(labels_list[::-1]))
    
    # Make text size smaller
    for i in range(len(fig_heatmap.layout.annotations)):
        fig_heatmap.layout.annotations[i].font.size = 9
        
    h = max((int(n/20)+1)*300,600)
    fig_heatmap['layout'].update({'template' :'plotly_white','height':h,'width':h,
                        'xaxis': {'title': '', 'showticklabels' : True, 'ticks' :''}, 
                        'yaxis': {'title': '', 'showticklabels' : True,'ticks' :''},
                         })

    #-----------------------------
    
    

    position_dendro = []

    size = len(fig_dendro.data)
    
    min_x=0   
    for i in range(size):           
        x = fig_dendro.data[i]['x']
        y = fig_dendro.data[i]['y']
        fig_dendro.data[i]['marker']['color']='grey'
        fig_dendro.data[i]['marker']['size']=0.2
        fig_dendro.data[i]['showlegend']=False
        
        min_x = min(min(x),min_x)
#         position_dendro += [[x,y]]  
    
    position_dendro_df = pd.DataFrame(labels_list, columns=['id'])
    position_dendro_df = position_dendro_df.sort_values(['id'])
    position_dendro_df = position_dendro_df.reset_index()

    position_dendro_df = pd.concat([position_dendro_df, clade_drop_dendro[['year','month','day','group','lineage','result']]], axis=1)

    
    #วาดเส้นแนวนอน
    x_new = []

    ratio = min(max(3,(int(n/5)+2)),5)
    scale = np.abs(min_x)/ratio
#     scale = 0.0005
    for _,row in position_dendro_df.iterrows():          
        y = int(row['year'])-2019
        m = int(row['month']) + (12*y)-12
        t = m*scale
        x_new += [t+(int(row['day'])/30)*scale]


    position_dendro_df['xnew'] = x_new
    position_dendro_df = position_dendro_df.sort_values(['index'])

    indications = sorted(list(set(clade_drop_dendro['group'])))
    group_dict = dict(zip(indications,range(len(indications))))
#     clade_color = {0:'lightcoral', 1:'lightseagreen', 2:'goldenrod', 3:'cornflowerblue', 4:'pink', 5:'plum', 6:'aquamarine', 7:'sandybrown', 8:'deepskyblue', 9:'slateblue', 10:'hotpink', 11:'gold', 12:'slategray', 13:'palegreen', 14:'tan', 15:'navy', 16:'darkorange', 17:'crimson', 18:'sadddlebrown', 19:'mediumseagreen', 20:'blueviolet', 21:'rosybrown'
# }

    clade_color = {0:'lightcoral',
                    1:'darkcyan',
                    2:'goldenrod', 
                    3:'olive',
                    4:'turquoise',
                    5:'deepskyblue',
                    6:'slategrey', 
                    7:'yellowgreen',
                    8:'hotpink',
                    9:'indianred',
                    10:'sandybrown',
                    11:'darkorange',
                    12:'plum',
                    13:'salmon',
                    14:'steelblue',
                    15:'cornflowerblue',
                    16:'slateblue',
                    17:'lightpink',
                    18:'lightgreen',
                    19:'dodgerblue',
                    20:'darkorchid',}

    color_list  = []
    label_list = []
    for i in X_id_list:
        c = clade_drop_dendro[clade_drop_dendro['id']==i]['group'].values[0]
        label_list += ['{}'.format(i)]
        color_list += [clade_color[group_dict[c]]]

    color_label = zip(label_list, color_list)
    dict_color_label = dict(color_label)

    color_list_sample = [dict_color_label[str(i)] for i in labels_list]
    position_dendro_df['color']=color_list_sample

    y_new = []

    for i in range(0,len(labels_list)):
        color_line = dict_color_label[str(labels_list[i])]
    #     print(color_line)
        xnew = position_dendro_df[position_dendro_df['index']==i]['xnew'].values[0]
        fig_dendro.add_trace(go.Scatter(x=[0,xnew], y=[5+i*10,5+i*10],mode='lines',hoverinfo='none',opacity=0.5,
                                    showlegend=False,
                                marker=dict(color=color_line, line=dict(width=0.01))
                            ))
        y_new += [5+i*10]

    #วาดจุด
    position_dendro_df['ynew'] = y_new

    #     summary_change_protein_df = pd.read_csv('../preprocessed/{}_summary_change_protein_df.csv'.format('asia')) 
    change_protein_df = change_protein_df[change_protein_df['id'].isin(X_id_list)][['id','change_protein','gene','check']]

    for i in range(n):
        x=position_dendro_df[position_dendro_df['id']==labels_list[i]]['xnew']#[position_dendro_df['xnew'][i]]
        y=position_dendro_df[position_dendro_df['id']==labels_list[i]]['ynew']#[position_dendro_df['ynew'][i]]
        id_sample = str(labels_list[i])
        
        

        snps_max_df = snps_max_df.sort_values('position')
        change_neucleotide_list = list(snps_max_df[snps_max_df['id']==labels_list[i]]['change'])
        change_neucleotide = alignment.get_mutation_text(change_neucleotide_list)


        change_protein_label = change_protein_df[(change_protein_df['id']==labels_list[i])&(change_protein_df['check']==False)]
        if len(list(change_protein_label['change_protein']))==0:
            change_protein = '-'
        else:
            change_protein = ''
            gene = list(set(change_protein_label['gene']))
            for g in ['S']: #gene:
                list_gene = list(set(change_protein_label[change_protein_label['gene']==g]['change_protein']))

                change_protein += f'<i>{g}</i>' + ' : ' + alignment.get_mutation_text(list_gene) + '</br> '



        color_marker = dict_color_label[str(id_sample)]
    #     print(color_marker)
        clade = position_dendro_df[position_dendro_df['id']==labels_list[i]]['result'].values[0]
        lin = position_dendro_df[position_dendro_df['id']==labels_list[i]]['lineage'].values[0]
        text_lin = id_sample +' ('+ lin + ')'
        g = position_dendro_df[position_dendro_df['id']==labels_list[i]]['group'].values[0]

        fig_dendro.add_trace(go.Scatter( x=x, y=y, text=text_lin,
                                mode='markers+text', name='', textposition="middle right",
                                        textfont=dict(size=8,),
                                hovertemplate = '<b>Id</b> : {}'.format(id_sample)+
                                        '<br><b>Group</b> : {}'.format(g)+
                                        '<br><b>Lineage</b> : {}'.format(lin)+
                                        '<br><b>Clade</b> : {}'.format(clade)+
                                        '<br><b>Mutation</b>:<br> {}'.format(change_neucleotide)+
                                        '<br><br><b>Amino acid change (Spike) </b>:<br> {}'.format(change_protein),
                                hoverlabel=dict(bgcolor='#F2F2F3'),
                                marker=dict(color=color_marker,size=8, line_width=1),
                                legendgroup=g, showlegend=False))

    for k,v in group_dict.items():
    #         print(k,v)
            fig_dendro.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                               marker=dict(size=8, color=clade_color[v], line_width=1),
                               legendgroup=k, showlegend=True, name=k))

    tickvals, ticktext = get_tick(scale)
    fig_dendro.update_layout(
        xaxis = dict(
            tickmode = 'array',
            tickvals = tickvals,
            ticktext = ticktext,
            tickangle = -45,
            ),
        margin=dict(l=0,r=0,b=0,t=50,pad=4)
    )
    fig_dendro.update_xaxes(showgrid=True)
    
    h = max(500,(int(n/30)+1)*300)

    fig_dendro['layout'].update({'template' :'plotly_white',
                             'height':h,
                        'xaxis': {'title': '', 'showticklabels' : True, 'ticks' :''}, 
                        'yaxis': {'title': '', 'showticklabels' : False,'ticks' :''},
                                })
    fig_dendro.update_layout(legend_font_size=10)
    fig_dendro.update_xaxes(tickfont_size=10)

#     fig_dendro.update_xaxes(range=(-0.001,position_dendro_df['xnew'].max()+scale*3))
    fig_dendro.update_xaxes(range=(1.5*min_x,position_dendro_df['xnew'].max()+np.abs(min_x)*1.5))
    fig_dendro.update_yaxes(range=(-20,position_dendro_df['ynew'].max()+25))
    fig_dendro.update_layout(showlegend=True)
    fig_dendro['layout'].update({'legend' : {'traceorder':'normal'}})
    
    return fig_dendro,fig_heatmap,position_dendro_df,dict_color_label,heatmap_df

#-------------------------------------------------------
