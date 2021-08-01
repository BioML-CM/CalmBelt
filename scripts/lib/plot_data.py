import pandas as pd
import numpy as np
import pickle

import plotly.express as px
import plotly.graph_objects as go
from lib import alignment,utils

clade_color = {0: 'cornflowerblue',
              1:'darkorange',
               2:'deepskyblue',
               3: 'pink',
               4:'plum',
               5:'lightcoral',
               6:'goldenrod',               
              7:'lightseagreen',
              8:'mediumseagreen'}


clade_who = {'Alpha': 'cornflowerblue',
              'Beta':'darkorange',
               'Gamma':'deepskyblue',
               'Delta': 'pink',
               'Eta':'plum',
               'Iota':'lightcoral',
               'Kappa':'goldenrod',               
              'Lambda':'lightseagreen',
              'Others':'mediumseagreen'}


def get_tick(scale):
    tickvals = [scale, scale*2, scale*3, scale*4, scale*5, scale*6, scale*7, scale*8,
                scale*9, scale*10, scale*11, scale*12, scale*13, scale*14, scale*15, scale*16, 
                scale*17,scale*18,scale*19,scale*20]
    ticktext = ['Jan-2020', 'Feb-2020', 'Mar-2020', 'Apr-2020', 'May-2020', 
                'Jun-2020', 'July-2020', 'Aug-2020', 'Sep-2020', 'Oct-2020', 'Nov-2020', 
                'Dec-2020', 'Jan-2021', 'Feb-2021', 'Mar-2021', 'Apr-2021', 'May-2021',
                'Jun-2021', 'July-2021', 'Aug-2021']
    return tickvals, ticktext



def fig_clade(sample_type):
    data_df = pd.read_csv('../preprocessed/{}_data_clade_df.csv'.format(sample_type), index_col=0)
    clade_name = sorted(list(data_df.columns[:-1]))
    
    fig_clade = go.Figure() 
    
    for i,l in enumerate(clade_name):
        fig_clade.add_trace(go.Scatter( 
            name = l, 
            x = data_df.index, 
            y = data_df[l], 
            mode='lines',
            line=dict(color=clade_color[i]),
            hovertemplate = '<b>Clade </b> {}'.format(l)+
                            '<br>%{x}'+
                            '<br>%{y:.2f}%',
            stackgroup='one'
           )) 
        
    fig_clade['layout'].update({'template' :'simple_white','width':800, 'height':400, 
#                           'title': '{} Trends (COVID-19)'.format(sample_type.upper()), 
                          'xaxis': {'title': 'Month'}, 
                          'yaxis': {'title': '%'},
                          'legend' : {'traceorder':'normal'}})

    
    fig_clade.update_layout(margin=dict(l=0,r=0,b=0,t=50,pad=4))

    return fig_clade



def fig_who(sample_type):
    data_df = pd.read_csv('../preprocessed/{}_data_who_df.csv'.format(sample_type), index_col=0)
    lin_name = list(data_df.columns[:-1])
    
    fig_lin = go.Figure() 
    
    for l in lin_name:
        fig_lin.add_trace(go.Scatter( 
            name = l, 
            x = data_df.index, 
            y = data_df[l], 
            line=dict(color=clade_who[l]),
            hovertemplate = '<b>Lineage </b> {}'.format(l)+
                            '<br>%{x}'+
                            '<br>%{y:.2f}%',
            stackgroup='one'
           )) 



    fig_lin['layout'].update({'template' :'simple_white','width':800, 'height':400, 
#                           'title': '{} Trends (COVID-19)'.format(sample_type.upper()), 
                          'xaxis': {'title': 'Month'}, 
                          'yaxis': {'title': '%'},
                           'legend' : {'traceorder':'normal'}})
    
    
    fig_lin.update_layout(margin=dict(l=0,r=0,b=0,t=50,pad=4))

    return fig_lin




def fig_sum(sample_type):
    data_df = pd.read_csv('../preprocessed/{}_data_clade_df.csv'.format(sample_type), index_col=0)
    fig_sum = go.Figure() 

    fig_sum.add_trace(go.Scatter( 
        name = '',
        x = data_df.index, 
        y = data_df['sum'], 
        hovertemplate = '<br>%{x}'+
                        '<br>%{y:.0f}',
        stackgroup='one'
       )) 


    fig_sum['layout'].update({'template' :'simple_white','width':800, 'height':400, 
                          'xaxis': {'title': 'Date'}, 
                          'yaxis': {'title': '# Submitted Genome'}})
    fig_sum.update_layout(margin=dict(l=0,r=0,b=0,t=50,pad=4))

    return fig_sum




from plotly.subplots import make_subplots

def fig_enp(sample_type):
    count_align_df = pd.read_csv('../preprocessed/{}_count_align_df.csv'.format(sample_type))
    count_align_df['pos'] = range(1,count_align_df.shape[0]+1)
    
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, shared_yaxes=True,
                        row_heights=[0.8, 0.2])

    fig.add_traces(go.Scatter(
                    x=count_align_df['pos'],
                    y=count_align_df['entropy'],name=''
                ), rows=1, cols=1)

    fig.add_traces(go.Scatter(name='', x=[265,265,21555,21555], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[21562,21562,25384,25384], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[25392,25392,26220,26220], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[26244,26244,26472,26472], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[26522,26522,27191,27191], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[27201,27201,27387,27387], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[27393,27393,27887,27887], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[27893,27893,28259,28259], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[28273,28273,29533,29533], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[29557,29557,29674,29674], y=[0,1,1,0], 
                              fill="toself",mode='lines',line_width=0)
                   , rows=2, cols=1)
    fig.add_traces(go.Scatter(name='', x=[29557,29557,29902,29902], y=[0,1,1,0], 
                               fill='toself', fillcolor='white',mode='lines',line_width=0)
                   , rows=2, cols=1)




    fig.add_annotation(dict(x=265+(21555-265)/2, y=0.4, 
                                text='ORF1ab', showarrow=False),row=2, col=1)
    fig.add_annotation(dict(x=21562+(25384-21562)/2, y=0.4, 
                                text='S', showarrow=False),row=2, col=1)
    fig.add_annotation(dict(x=25392+(26220-25392)/2, y=0.4, 
                                text='ORF3a', showarrow=False),row=2, col=1)
    fig.add_annotation(dict(x=26244+(26472-26244)/2, y=0.4,
                                text='E', showarrow=False),row=2, col=1)
    fig.add_annotation(dict(x=26522+(27191-26522)/2, y=0.4,
                                text='M', showarrow=False),row=2, col=1)
    fig.add_annotation(dict(x=27201+(27387-27201)/2, y=0.4, 
                                text='ORF6', showarrow=False),row=2, col=1)
    fig.add_annotation(dict(x=27393+(27887-27393)/2, y=0.4,
                                text='ORF7ab', showarrow=False),row=2, col=1)
    fig.add_annotation(dict(x=27893+(28259-27893)/2, y=0.4,
                                text='ORF8', showarrow=False),row=2, col=1)
    fig.add_annotation(dict(x=28273+(29533-28273)/2, y=0.4, 
                                text='N', showarrow=False),row=2, col=1)
    fig.add_annotation(dict(x=29557+(29674-29557)/2, y=0.4,
                                text='ORF10', showarrow=False),row=2, col=1)


    fig.update_layout(title_text=None,
    #                 title_font_size=16 ,
                    template='simple_white',height=400,
                    xaxis1=dict(showticklabels=False,ticks='',
                              hoverformat="{n}f" ),
                     yaxis1=dict(
                         title_text= 'Diversity',
                               hoverformat=".4f"),
                     showlegend=False)

    fig.update_layout(xaxis2=dict(
                            rangeslider=dict(visible=True,),
                            tickmode = 'array',
                            ticks='outside',
                            tickvals = [2000, 4000, 6000, 8000, 10000, 12000,14000, 16000, 18000, 20000, 22000, 24000, 26000, 28000],
                            ticktext = ['2000', '4000', '6000', '8000', '10000', '12000','14000', '16000', '18000', '20000', '22000',
                            '24000', '26000', '28000'],
                             title_text= 'Position',
                              hoverformat="{n}f" ),
                      yaxis2=dict(showticklabels=False, ticks='', visible=False),
                      margin=dict(l=0,r=0,b=0,t=50,pad=4),
                      showlegend=False)
    return fig







def fig_cluster(sample_type, cluster_type):
    Y_tse = pickle.load(open("../preprocessed/{}_Y_tse.pickle".format(sample_type), "rb"))
    Y_tse = Y_tse.sort_values([cluster_type])
    
    color_list = ['cornflowerblue','darkorange','deepskyblue','pink','plum','lightcoral','goldenrod','lightseagreen','mediumseagreen']
    
    if cluster_type == 'k_mean':
        Y_tse['k-mean'] = 'C'+ Y_tse['k_mean']
        cluster_type = 'k-mean'
    elif cluster_type == 'result':
        Y_tse['Clade'] = Y_tse['result']
        cluster_type = 'Clade'
    elif cluster_type == 'who_name':
        Y_tse['WHO'] = Y_tse['who_name']
        cluster_type = 'WHO'
        color_list = ['cornflowerblue','darkorange','pink','plum','deepskyblue','lightcoral','goldenrod','lightseagreen','mediumseagreen']
        
    clade = Y_tse['{}'.format(cluster_type)].values
     
    
    symbol_list = ['circle','cross','x','triangle-up','triangle-down','square','diamond']
    symbol_open = ['circle-open','cross-open','x-open','triangle-up-open','triangle-down-open','square-open','diamond-open']

    fig_cluster = px.scatter(data_frame=Y_tse, x='x', y='y', color=cluster_type, symbol=cluster_type, 
                             symbol_sequence=symbol_open, opacity=0.75,
                             color_discrete_sequence = color_list,
                             hover_data={'x': False, 'y': False, cluster_type: False, 
                                        'Cluster': (clade)})
    fig_cluster.update_traces(marker=dict(size=6,
                   line=dict(width=1
                )),
                      selector=dict(mode='markers'))

#     fig_cluster.update_layout(title='Clustering : color = {}'.format(cluster_type))
    fig_cluster['layout'].update({'template' :'simple_white', 
                                  'width':450, 'height':400,
                            'xaxis': {'title': 't-SNE 1', 'showticklabels' : False, 'mirror': True}, 
                            'yaxis': {'title': 't-SNE 2', 'showticklabels' : False, 'mirror': True},
                            'legend' : {'traceorder':'normal'}})
    fig_cluster.update_layout(margin=dict(l=0,r=0,b=0,t=50,pad=4))
    
    return fig_cluster




from plotly.figure_factory import create_dendrogram
from scipy.spatial.distance import pdist

def fig_dendro(sample_type,color_style):
    
    X_dendro = pickle.load(open('../preprocessed/{}_X_dendro.pickle'.format(sample_type), "rb"))
    X_id_list = pickle.load(open('../preprocessed/{}_X_id_list_dendro.pickle'.format(sample_type), "rb"))
    change_protein_df = pd.read_csv('../preprocessed/{}_summary_change_protein_df.csv'.format(sample_type)) 
    snps_max_df = pd.read_csv('../preprocessed/{}_snps_max_df.csv'.format(sample_type)) 

    clade_drop_dendro = pickle.load(open('../preprocessed/{}_clade_drop_dendro.pickle'.format(sample_type), "rb"))
    
    change_protein_df = change_protein_df[change_protein_df['id'].isin(X_id_list)][['id','change_protein','gene','check']]
    snps_max_df = snps_max_df[snps_max_df['id'].isin(X_id_list)]

    n = len(X_id_list)
    
    fig_dendro = create_dendrogram(X_dendro, orientation='right', labels=X_id_list,
                           distfun=lambda x: pdist(x,metric='hamming'))
    
    labels_list = fig_dendro.layout['yaxis']['ticktext'] #label เรียงตามลำดับ
    
    position_dendro = []

    size = len(fig_dendro.data)

    for i in range(size):    
        x = fig_dendro.data[i]['x']
        y = fig_dendro.data[i]['y']
        fig_dendro.data[i]['marker']['color']='grey'
        fig_dendro.data[i]['marker']['size']=0.2
        fig_dendro.data[i]['showlegend']=False
        position_dendro += [[x,y]]  


    position_dendro_df = pd.DataFrame(labels_list, columns=['id'])
    position_dendro_df = position_dendro_df.sort_values(['id'])
    position_dendro_df = position_dendro_df.reset_index()

    position_dendro_df = pd.concat([position_dendro_df, clade_drop_dendro[['date','year','month','day',color_style]]], axis=1)

    #วาดเส้นแนวนอน
    x_new = []
    scale = 0.5
    for _,row in position_dendro_df.iterrows():          
        y = int(row['year'])-2019
        m = int(row['month']) + (12*y)-12
        t = m*scale
        x_new += [t+(int(row['day'])/30)*scale]
        


    position_dendro_df['xnew'] = x_new
    position_dendro_df = position_dendro_df.sort_values(['index'])

    indications = sorted(list(set(clade_drop_dendro[color_style])))
#     print(indications)
    group_dict = dict(zip(indications,range(len(indications))))
    


    color_list  = []
    label_list = []
    if color_style == 'who_name':
        for i in X_id_list:
            c = clade_drop_dendro[clade_drop_dendro['id']==i][color_style].values[0]
            label_list += ['{}'.format(i)]
            color_list += [clade_who[c]]
    else:
        for i in X_id_list:
            c = clade_drop_dendro[clade_drop_dendro['id']==i][color_style].values[0]
            label_list += ['{}'.format(i)]
            color_list += [clade_color[group_dict[c]]]


    color_label = zip(label_list, color_list)
    dict_color_label = dict(color_label)

    color_list_sample = [dict_color_label[str(i)] for i in labels_list]
    position_dendro_df['color'] = color_list_sample
    
    opacity_dict={'asia':0.03,'thailand':0.1,'singapore':0.1}
    y_new = []

    for i in range(0,len(labels_list)):
        color_line = str(dict_color_label[str(labels_list[i])])
        xnew = position_dendro_df[position_dendro_df['index']==i]['xnew'].values[0]
        fig_dendro.add_trace(go.Scatter(x=[0,xnew], y=[5+i*10,5+i*10],mode='lines',hoverinfo='none',opacity=opacity_dict[sample_type],
                                    showlegend=False,
                                marker=dict(color=color_line, line=dict(
                    width=0.01
                ))))
        y_new += [5+i*10]


    #วาดจุด
    position_dendro_df['ynew'] = y_new
    
    

    for i in range(n):
        x=position_dendro_df[position_dendro_df['id']==labels_list[i]]['xnew']#[position_dendro_df['xnew'][i]]
        y=position_dendro_df[position_dendro_df['id']==labels_list[i]]['ynew']#[position_dendro_df['ynew'][i]]
        id_sample = str(labels_list[i])
        
        snps_max_df = snps_max_df.sort_values('position')
        change_neucleotide_list = list(snps_max_df[snps_max_df['id']==labels_list[i]]['change'])
        change_neucleotide = utils.get_mutation_text(change_neucleotide_list)
 
            
        change_protein_label = change_protein_df[(change_protein_df['id']==labels_list[i])&(change_protein_df['check']==False)]
        if len(list(change_protein_label['change_protein']))==0:
            change_protein = '-'
        else:
            change_protein = ''
            gene = list(set(change_protein_label['gene']))
            for g in gene:
                list_gene = list(set(change_protein_label[change_protein_label['gene']==g]['change_protein']))

                change_protein += f'<i>{g}</i>' + ' : ' + utils.get_mutation_text(list_gene) + '</br> '
            
            
        
        color = str(dict_color_label[str(id_sample)])
        if color_style == 'k_mean':
            clade = position_dendro_df[position_dendro_df['id']==labels_list[i]][color_style].values[0].astype(str)
        else:
            clade = position_dendro_df[position_dendro_df['id']==labels_list[i]][color_style].values[0]
            
        fig_dendro.add_trace(go.Scatter( x=x, y=y,
                                mode='markers', name='', 
                                hovertemplate = '<b>Id</b> : {}'.format(id_sample)+
                                        '<br><b>Clade</b> : {}'.format(clade)+
                                        '<br><b>Mutation</b>:<br> {}'.format(change_neucleotide)+
                                        '<br><br><b>Amino acid change</b>:<br> {}'.format(change_protein),
                                hoverlabel=dict(bgcolor='#F2F2F3'),
                                marker=dict(color=dict_color_label[str(id_sample)],size=8, line_width=1),
                                legendgroup=clade, showlegend=False))
    
    if color_style == 'who_name':
        for k,v in clade_who.items():
            fig_dendro.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                               marker=dict(size=8, color=clade_who[k], line_width=1),
                               legendgroup=k, showlegend=True, name=k))                    
    else:        
        for k,v in group_dict.items():
            fig_dendro.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                               marker=dict(size=8, color=clade_color[v], line_width=1),
                               legendgroup=k, showlegend=True, name=k))

    
    tickvals, ticktext = get_tick(scale)
    fig_dendro.update_layout(
        xaxis = dict(
            tickmode = 'array',
            tickvals = tickvals,
            ticktext = ticktext,
            tickangle = -45
    ))
    fig_dendro.update_xaxes(showgrid=True)

    fig_dendro['layout'].update({'template' :'plotly_white','height':800,
                            'xaxis': {'title': '', 'showticklabels' : True, 'ticks' :''}, 
                            'yaxis': {'title': '', 'showticklabels' : False,'ticks' :''},
                             'legend' : {'traceorder':'normal'}})
    fig_dendro.update_layout(margin=dict(l=0,r=0,b=0,t=50,pad=4))
    
    if sample_type == 'asia':
        fig_dendro.update_yaxes(showticklabels=False, range=(position_dendro_df['ynew'].min()-1200,position_dendro_df['ynew'].max()+1200))
        fig_dendro.update_xaxes(range=(position_dendro_df['xnew'].min()-2,position_dendro_df['xnew'].max()+0.3))
    else:   
        fig_dendro.update_yaxes(showticklabels=False, range=(position_dendro_df['ynew'].min()-300,position_dendro_df['ynew'].max()+300))
        fig_dendro.update_xaxes(range=(position_dendro_df['xnew'].min()-2,position_dendro_df['xnew'].max()+0.3))
    fig_dendro.update_layout(showlegend=True)

    return fig_dendro,position_dendro_df,dict_color_label






def alarm_plot(df,position,column):
    fig = go.Figure()

    for i in position:
        
        plot_df = df[df[column]==int(i)]
        
        sum_plot_df = plot_df.groupby(['month','total',column]).sum()
        sum_change = []
        for m in sum_plot_df.index:
            sum_change += [plot_df[plot_df['month']==m[0]]['graph']]
        sum_plot_df['sum_change'] = sum_change
        sum_plot_df = sum_plot_df.reset_index()

        if sum_plot_df.shape[0]>1:
            fig.add_trace(go.Scatter(x=sum_plot_df['month'], y=sum_plot_df['id'], customdata=sum_plot_df, hovertemplate =
                                        '<b>Position</b>: {}'.format(i)+
        #                                 '<br><b>Gene</b>: {}'.format(gene)+
                                        '<br><b>Month</b>: %{x}'+
                                        '<br><b>#pateint</b>: %{y}'+
                                        ' <b> of </b> %{customdata[1]}'+
                                         '<br><b>Change = </b> %{customdata[4]}'+
                                        '<br><b>occured in </b>: {} month'.format(sum_plot_df.shape[0]),
                                    name='{}'.format(i), mode='lines+markers',
        #                             marker=dict(color=clade_color[str(clade)],size=5, line_width=1),
        #                             legendgroup=clade, 
                                     showlegend=True
                                    ))


    fig['layout'].update({'template' :'simple_white','width':600, 'height':400,
                            'xaxis': {'title': 'Month', 'showticklabels' : True, 'ticks' :''}, 
                            'yaxis': {'title': 'Number of cases', 'showticklabels' : True,'ticks' :''}})
    
    tickvals, ticktext = get_tick(1)
    fig.update_layout(
            xaxis = dict(
                tickmode = 'array',
                tickvals = tickvals,
                ticktext = ticktext,
                tickangle = -45,),
            yaxis = dict( type='log'),
            margin=dict(l=0,r=0,b=0,t=50,pad=4),
    )

    return fig


#for insertion deletion
def alarm_plot2(df,position,column):
    fig = go.Figure()

    for i in position:
        
        plot_df = df[df[column]==int(i)]
        fig.add_trace(go.Scatter(x=plot_df['month'], y=plot_df['id'],text=plot_df['total'],
                            hovertemplate =
                                    '<b>Position</b>: {}'.format(i)+
    #                                 '<br><b>Gene</b>: {}'.format(gene)+
                                    '<br><b>Month</b>: %{x}'+
                                    '<br><b>#pateint</b>: %{y}'+
                                 ' <b> of </b> %{text}'+
    #                                 '<br><b>Percent</b>: %{y:.2f}'+
                                    '<br><b>occured in </b>: {} month'.format(plot_df.shape[0]),
                                name='{}'.format(i), mode='lines+markers',
    #                             marker=dict(color=clade_color[str(clade)],size=5, line_width=1),
    #                             legendgroup=clade, 
                                 showlegend=True
                                ))


    fig['layout'].update({'template' :'simple_white','width':600, 'height':400,
                            'xaxis': {'title': 'Month', 'showticklabels' : True, 'ticks' :''}, 
                            'yaxis': {'title': 'Number of cases', 'showticklabels' : True,'ticks' :''}})
    
    tickvals, ticktext = get_tick(1)
    
    fig.update_layout(
            xaxis = dict(
                tickmode = 'array',
                tickvals = tickvals,
                ticktext = ticktext,
                tickangle = -45,
                ),
        yaxis = dict( type='log'
                ),
    )
    return fig


def alarm_excel_plot(gene_list,pos_list):
    df = pickle.load(open('../preprocessed/{}_protein_alarm.pickle'.format('singapore'), "rb"))
    
    fig = go.Figure()
    for i in range(len(gene_list)):
        protein_plot = df[str(gene_list[i])]
        plot_df = protein_plot[protein_plot['change_protein']==pos_list[i]]
#         print(plot_df)
        fig.add_trace(go.Scatter(x=plot_df['month'], y=plot_df['id'],customdata=plot_df,
                                hovertemplate =
                                        '<b>Position</b>: %{customdata[1]}'+
                                        '<br><b>Gene</b>: %{customdata[4]}'+
                                        '<br><b>Month</b>: %{x}'+
                                        '<br><b>#pateint</b>: %{y}'+
                                         ' <b> of </b> %{customdata[5]}'+
        #                                 '<br><b>Percent</b>: %{y:.2f}'+
                                        '<br><b>occured in </b>: {} month'.format(plot_df.shape[0]),
                                    name='{}'.format(pos_list[i]), mode='lines+markers',
        #                             marker=dict(color=clade_color[str(clade)],size=5, line_width=1),
        #                             legendgroup=clade, 
                                     showlegend=True

                                    ))
        fig['layout'].update({'template' :'simple_white','width':600, 'height':500,
                                'xaxis': {'title': 'Month', 'showticklabels' : True, 'ticks' :''}, 
                                'yaxis': {'title': 'Number of cases', 'showticklabels' : True,'ticks' :''}})
        
        tickvals, ticktext = get_tick(1)
        fig.update_layout(
                xaxis = dict(
                    tickmode = 'array',
                    tickvals = tickvals,
                    ticktext = ticktext,
                    tickangle = -45,
                    ),
            yaxis = dict( type='log'
                    ),
        )
    return fig

def alarm_specific_protein(gene_list,pos_list):
    #load data
    id_month = pickle.load(open("../preprocessed/{}_id_month.pickle".format('singapore'), "rb"))
    id_month_key = sorted(list(id_month.keys()))

    summary_change_protein_df = pd.read_csv('../preprocessed/{}_summary_change_protein_df.csv'.format('singapore')) 
    
    
    change_protein_df = summary_change_protein_df.copy()
    change_protein_df = change_protein_df[change_protein_df['check']==False]
    protein_m_df = {}
    plot_list = []
    i=0
    for m in id_month_key:
        protein_m_df[m] = change_protein_df[change_protein_df['id'].isin(id_month[m].astype(int))]
        protein_m_df[m] = protein_m_df[m][['id','position_protein','change_protein','gene']]

        df= protein_m_df[m]
        for k in range(len(gene_list)):
            id_match = df[(df['change_protein']==pos_list[k].upper()) & (df['gene']==gene_list[k])]['id'].values
            df = df[df['id'].isin(id_match)]

        number_id_match = len(set(df['id']))
        total = len(id_month[m])
        plot_list += [[i+1,number_id_match,total]]
        i+=1

    plot_df = pd.DataFrame(plot_list, columns=['month', 'number', 'total'])
    
    
    #plot
    fig = go.Figure()

    fig.add_trace(go.Bar(x=plot_df['month'], y=plot_df['number'],text=plot_df['total'],
                            hovertemplate =
                                    '<br><b>Month</b>: %{x}'+
                                    '<br><b>#pateint</b>: %{y}'+
                                  ' <b> of </b> %{text}',
                                name='',
                                 showlegend=False

                                ))
    fig['layout'].update({'template' :'simple_white','width':600, 'height':400,
                                'xaxis': {'title': 'Month', 'showticklabels' : True, 'ticks' :''}, 
                                'yaxis': {'title': 'Number of cases', 'showticklabels' : True,'ticks' :''}})
    
    tickvals, ticktext = get_tick(1)
    fig.update_layout(
                xaxis = dict(
                    tickmode = 'array',
                    tickvals = tickvals,
                    ticktext = ticktext,
                    tickangle = -45,
                    ),
            yaxis = dict( type='log'
                    ),
        )
    return fig