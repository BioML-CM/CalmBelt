#!/usr/bin/env python
# coding: utf-8

# In[1]:


import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

import dash_table as dt

from flask import Flask, send_from_directory

import time
import datetime

import base64
import os
from urllib.parse import quote as urlquote


UPLOAD_DIRECTORY = "../upload/"
DOWNLOAD_DIRECTORY = "../download"

if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)


# In[2]:


server = Flask(__name__)
app = dash.Dash(__name__, suppress_callback_exceptions=True, server=server)
app.title = 'SARS-CoV-2 Analysis'


# In[3]:


import shutil #for copy file + delete folder
import time, os, fnmatch, shutil
import subprocess


# In[4]:


pango_version = subprocess.check_output(['pangolin', '--version']).decode().strip()
# pango_version


# In[5]:


@server.route("/upload/<path:path>")
def download(path):
    """Serve a file from the download directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)

@server.route("/download/<path:path>")
def download2(path):
    """Serve a file from the download directory."""
    return send_from_directory(DOWNLOAD_DIRECTORY, path, as_attachment=True)

def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "../upload/{}".format(urlquote(filename))
    return html.A(filename.split('/')[2], href=location)


# In[9]:


import pandas as pd
import numpy as np
import pickle

import plotly.express as px
import plotly.graph_objects as go

from lib import plot_data, built, utils


# # Top  lineage

# In[10]:


count_pango_df = pd.read_csv('../preprocessed/count_pango_df.csv', index_col=0)
count_pango_df['%'] = count_pango_df['count']/count_pango_df['total']*100


# In[11]:


region_lineage = {}

for r in ['North America', 'Asia', 'Europe', 'South America', 'Oceania', 'Africa']:
    temp_df = count_pango_df.loc[count_pango_df['region']==r,['lineage','%']].values.tolist()
    region_lineage[r] = ['{} : {}%'.format(l,round(c, 2)) for l,c in temp_df]    


# In[12]:


count_pango_country_df = pd.read_csv('../preprocessed/count_pango_country_df.csv', index_col=0)
count_pango_country_df['%'] = count_pango_country_df['count']/count_pango_country_df['total']*100


# In[13]:


country_lineage = {}

for r in ['Singapore']:
    temp_df = count_pango_country_df.loc[count_pango_country_df['country']==r,['lineage','%']].values.tolist()
    country_lineage[r] = ['{} : {}%'.format(l,round(c, 2)) for l,c in temp_df]


# # Plot figure

# In[14]:


# ['North America', 'Asia', 'Europe', 'South America', 'Oceania', 'Africa']
fig_clade_north_america = plot_data.fig_clade('North America')
fig_clade_asia = plot_data.fig_clade('Asia')
fig_clade_europe = plot_data.fig_clade('Europe')
fig_clade_south_america = plot_data.fig_clade('South America')
fig_clade_oceania = plot_data.fig_clade('Oceania')
fig_clade_africa = plot_data.fig_clade('Africa')

fig_who_north_america = plot_data.fig_who('North America')
fig_who_asia = plot_data.fig_who('Asia')
fig_who_europe = plot_data.fig_who('Europe')
fig_who_south_america = plot_data.fig_who('South America')
fig_who_oceania = plot_data.fig_who('Oceania')
fig_who_africa = plot_data.fig_who('Africa')


fig_sum_north_america = plot_data.fig_sum('North America')
fig_sum_asia = plot_data.fig_sum('Asia')
fig_sum_europe = plot_data.fig_sum('Europe')
fig_sum_south_america = plot_data.fig_sum('South America')
fig_sum_oceania = plot_data.fig_sum('Oceania')
fig_sum_africa = plot_data.fig_sum('Africa')


# In[15]:


fig_clade_sin = plot_data.fig_clade('Singapore')
fig_sum_sin = plot_data.fig_sum('Singapore')
fig_who_sin = plot_data.fig_who('Singapore')


# In[16]:


fig_enp_sin = plot_data.fig_enp('singapore')
fig_cluster_k_mean_sin = plot_data.fig_cluster('singapore','k_mean')
# fig_cluster_result_sin = plot_data.fig_cluster('singapore','result')
fig_cluster_who_sin = plot_data.fig_cluster('singapore','who_name')


# In[17]:


fig_cluster_k_mean_sin.write_image("../download/fig_cluster_k_mean_sin.png", scale=3)
# fig_cluster_result_sin.write_image("../download/fig_cluster_result_sin.png", scale=3)
fig_cluster_who_sin.write_image("../download/fig_cluster_who_sin.png", scale=3)


# In[18]:


# fig_dendro_sin,_,_ = plot_data.fig_dendro('singapore','result')
fig_dendro_sin,_,_ = plot_data.fig_dendro('singapore','who_name')


# # About

# In[19]:


def call_image(image_filename):
    encoded_image = base64.b64encode(open(image_filename, 'rb').read())
    return encoded_image.decode()


about_layout = html.Div([
                        html.Br(), 
                        html.H1('SARS-CoV-2 genomic characterization'),
                        html.P(['''The current pandemic caused by SARS-CoV-2 virus raises great concern around the world. 
                                Analysis of SARS-CoV-2 genomic information suggests rapid mutation rates since the start of the outbreak in December 2019. 
                                The very first genomic data from Wuhan is typically used as ''',
                                html.A(children='a standard reference genome.', href='https://www.ncbi.nlm.nih.gov/sars-cov-2/', target='_blank'), 
                                ''' GISAID has initially classified SARS-CoV-2 genomes into different clades (i.e. a group of genomes with a similar genetic profile). 
                                These include S (considered as an original clade), L, V and G, as well as other subclades. Clade information plays an important role in outbreak tracking and infectious properties. ''',
                                html.A(children='A recent study', href='https://www.nature.com/articles/d41586-020-02544-6', target='_blank'), 
                                ''' has reported that patients infected with the clade G virus show a larger amount of virus at the upper neck, allowing the virus to spread easier.''']),
                        html.Br(),            
    
                        html.P(['''Although qPCR has typically been used for rapid virus detection, 
                        it does not allow us to study the RNA sequence of the virus to determine clades. 
                        The advance of sequencing technology and the increasing accessibility allow scientists to generate more than tens of thousands of virus genome sequences. 
                        An intuitive graphical user interface for analysing multiple genomes is essential for end-users in hospitals, facilitating rapid detection of variants of concerns and outbreak tracking.''']),
                        html.Br(),            
    
                        html.P([                            
                        ''' CalmBelt’s interactive interface allows users to analyse multiple SARS-CoV-2 genomes by utilising whole genome information, collection date, and additional information such as predefined clusters based on contact tracing. 
                        CalmBelt also integrated external SARS-CoV-2 nomenclature assignments, GISAID and PANGO, allowing users to visualise whole genome similarity and nomenclatures at the same time.''']),  
                        
                        html.Hr(), 
                        html.H1(['''SARS-CoV-2 nomenclatures''']), 
                        html.P([html.A(children='SARS-CoV-2  GISAID Clades', href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5388101/', target='_blank')]),
                                html.Br(),

                            html.P(html.Table([
                                    html.Tr([html.Th('Clades'),html.Th('Mutations') ]),
                                    html.Tr([html.Td('L'),html.Td('C241, C3037, C8782, G11083, C22227, A23403, G25563, G26144, T28144, G28882')]),
                                    html.Tr([html.Td('S'),html.Td('C8782T, T28144C')]),
                                    html.Tr([html.Td('V'),html.Td('G11083T, G26144T')]),
                                    html.Tr([html.Td('G'),html.Td('C241T, C3037T, A23403G')]),
                                    html.Tr([html.Td('GH'),html.Td('C241T, C3037T, A23403G, G25563T')]),
                                    html.Tr([html.Td('GR'),html.Td('C241T, C3037T, A23403G, G28882A')]),
                                    html.Tr([html.Td('GRY'),html.Td('C241T, C3037T, A23063T, A23403G, G28882A, 21765-21770del, 21991-21993del')]),
                                    html.Tr([html.Td('GV'),html.Td('C241T, C3037T, A23403G, C22227T')]),
                                    
                                    html.Tr([html.Td('O'),html.Td('Others')])
                                ])),
                        html.Br(), 
                        html.P([html.A(children='Pango Lineages', href='https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/', target='_blank', )])
                                ,html.Br(),   
                            html.P(html.Table([
                                html.Tr([html.Th('WHO label'),html.Th('Lineage'),html.Th('GISAID clade')]),
                                html.Tr([html.Td('Alpha'),html.Td('B.1.1.7'),html.Td('GRY')]),
                                html.Tr([html.Td('Beta'),html.Td('B.1.351, B.1.351.2, B.1.351.3'),html.Td('GH')]),
                                html.Tr([html.Td('Gamma'),html.Td('P.1, P.1.1, P.1.2'),html.Td('GR')]),
                                html.Tr([html.Td('Delta'),html.Td('B.1.617.2, AY.1, AY.2, AY.3'),html.Td('G')]),
                                html.Tr([html.Td('Eta'),html.Td('B.1.525'),html.Td('G')]),
                                html.Tr([html.Td('Iota'),html.Td('B.1.526'),html.Td('GH')]),
                                html.Tr([html.Td('Kappa'),html.Td('B.1.617.1'),html.Td('G')]),
                                html.Tr([html.Td('Lambda'),html.Td('C.37'),html.Td('GR')]),
                            ])),
    
                        html.Br(), 
                        html.Br(),
                        html.Br(), 
                        html.Br(),
    ])


# # All Region

# In[20]:


# ['North America', 'Asia', 'Europe', 'South America', 'Oceania', 'Africa']

region_layout = html.Div([
                        html.H3('Asia Trends (COVID-19)'),
                        html.P(['''Overall statistics of COVID-19 cases in Asia. A time-series plot shows the percentage and number of cases for each clade or lineage (7-day average)''']),

                         html.Div(children=[
                                            dcc.Loading(
                                                id="loading-plot-asia",
                                                type="default",
                                                children=html.Div(id="total-data-plot")
                                                       )
                                      ]), 
                    
                       html.Div(children=[
                                            dcc.Loading(
                                                id="who-plot-asia",
                                                type="default",
                                                children=html.Div(id="who-data-plot")
                                                       )
                                      ]),
                        html.Div(children=[dcc.RadioItems(id='total-data',
                                    options=[{'label': 'WHO naming (%)', 'value': 'who'},
                                             {'label': 'GISAID clade (%)', 'value': 'clade'},
                                             {'label': 'Number of cases', 'value': 'number'},],
                                    value='who',
                                    labelStyle={'display': 'inline-block'}
                                    ) 
                                 ]),
                        
                        html.P(['''Top 5 PANGO lineages in last 30 days = ''',', '.join(region_lineage['Asia'])]),
                        html.Hr(),
                        
    
                        html.H3('Europe Trends (COVID-19)'),
                        html.P(['''Overall statistics of COVID-19 cases in Europe. A time-series plot shows the percentage and number of cases for each clade or lineage (7-day average)''']),

                         html.Div(children=[
                                            dcc.Loading(
                                                id="loading-plot-Europe",
                                                type="default",
                                                children=html.Div(id="total-data-plot-Europe")
                                                       )
                                      ]), 
                        html.Div(children=[
                                            dcc.Loading(
                                                id="who-plot-Europe",
                                                type="default",
                                                children=html.Div(id="who-data-plot-Europe")
                                                       )
                                      ]), 
                        html.Div(children=[dcc.RadioItems(id='total-data-Europe',
                                    options=[{'label': 'WHO naming (%)', 'value': 'who'},
                                             {'label': 'GISAID clade (%)', 'value': 'clade'},
                                             {'label': 'Number of cases', 'value': 'number'},],
                                    value='who',
                                    labelStyle={'display': 'inline-block'}
                                    ) 
                                 ]),
                        
                        html.P(['''Top 5 PANGO lineages in last 30 days = ''',', '.join(region_lineage['Europe'])]),
                        html.Hr(),
    
    
                        html.H3('North America Trends (COVID-19)'),
                        html.P(['''Overall statistics of COVID-19 cases in North America. A time-series plot shows the percentage and number of cases for each clade or lineage (7-day average)''']),

                         html.Div(children=[
                                            dcc.Loading(
                                                id="loading-plot-North-America",
                                                type="default",
                                                children=html.Div(id="total-data-plot-North-America")
                                                       )
                                      ]), 
                    
                        html.Div(children=[
                                            dcc.Loading(
                                                id="who-plot-North-America",
                                                type="default",
                                                children=html.Div(id="who-data-plot-North-America")
                                                       )
                                      ]), 
                        html.Div(children=[dcc.RadioItems(id='total-data-North-America',
                                    options=[{'label': 'WHO naming (%)', 'value': 'who'},
                                             {'label': 'GISAID clade (%)', 'value': 'clade'},
                                             {'label': 'Number of cases', 'value': 'number'},],
                                    value='who',
                                    labelStyle={'display': 'inline-block'}
                                    ) 
                                 ]),
                        
                        html.P(['''Top 5 PANGO lineages in last 30 days = ''',', '.join(region_lineage['North America'])]),
                        html.Hr(),
    
    
                        html.H3('South America Trends (COVID-19)'),
                        html.P(['''Overall statistics of COVID-19 cases in South America. A time-series plot shows the percentage and number of cases for each clade or lineage (7-day average)''']),

                         html.Div(children=[
                                            dcc.Loading(
                                                id="loading-plot-South-America",
                                                type="default",
                                                children=html.Div(id="total-data-plot-South-America")
                                                       )
                                      ]), 
                        html.Div(children=[
                                            dcc.Loading(
                                                id="who-plot-South-America",
                                                type="default",
                                                children=html.Div(id="who-data-plot-South-America")
                                                       )
                                      ]), 

                        html.Div(children=[dcc.RadioItems(id='total-data-South-America',
                                    options=[{'label': 'WHO naming (%)', 'value': 'who'},
                                             {'label': 'GISAID clade (%)', 'value': 'clade'},
                                             {'label': 'Number of cases', 'value': 'number'},],
                                    value='who',
                                    labelStyle={'display': 'inline-block'}
                                    ) 
                                 ]),
                        
                        html.P(['''Top 5 PANGO lineages in last 30 days = ''',', '.join(region_lineage['South America'])]),
                        html.Hr(),

    
                        html.H3('Africa Trends (COVID-19)'),
                        html.P(['''Overall statistics of COVID-19 cases in Africa. A time-series plot shows the percentage and number of cases for each clade or lineage (7-day average)''']),

                         html.Div(children=[
                                            dcc.Loading(
                                                id="loading-plot-Africa",
                                                type="default",
                                                children=html.Div(id="total-data-plot-Africa")
                                                       )
                                      ]), 
                        html.Div(children=[
                                            dcc.Loading(
                                                id="who-plot-Africa",
                                                type="default",
                                                children=html.Div(id="who-data-plot-Africa")
                                                       )
                                      ]), 
                        html.Div(children=[dcc.RadioItems(id='total-data-Africa',
                                    options=[{'label': 'WHO naming (%)', 'value': 'who'},
                                             {'label': 'GISAID clade (%)', 'value': 'clade'},
                                             {'label': 'Number of cases', 'value': 'number'},],
                                    value='who',
                                    labelStyle={'display': 'inline-block'}
                                    ) 
                                 ]),
                        
                        html.P(['''Top 5 PANGO lineages in last 30 days = ''',', '.join(region_lineage['Africa'])]),
                        html.Hr(),
    
    
                        html.H3('Oceania Trends (COVID-19)'),
                        html.P(['''Overall statistics of COVID-19 cases in Oceania. A time-series plot shows the percentage and number of cases for each clade or lineage (7-day average)''']),

                         html.Div(children=[
                                            dcc.Loading(
                                                id="loading-plot-Oceania",
                                                type="default",
                                                children=html.Div(id="total-data-plot-Oceania")
                                                       )
                                      ]), 
                        html.Div(children=[
                                            dcc.Loading(
                                                id="who-plot-Oceania",
                                                type="default",
                                                children=html.Div(id="who-data-plot-Oceania")
                                                       )
                                      ]), 
                    
                        html.Div(children=[dcc.RadioItems(id='total-data-Oceania',
                                    options=[{'label': 'WHO naming (%)', 'value': 'who'},
                                             {'label': 'GISAID clade (%)', 'value': 'clade'},
                                             {'label': 'Number of cases', 'value': 'number'},],
                                    value='who',
                                    labelStyle={'display': 'inline-block'}
                                    ) 
                                 ]),
                        
                        html.P(['''Top 5 PANGO lineages in last 30 days = ''',',  '.join(region_lineage['Oceania'])]),
                        html.Hr(),
                        
                        html.P(['''Last update : 30 June 2021''']),
                        html.Br(), 
                        html.Br(),
                        

    ])

                                

@app.callback(Output('total-data-plot', 'children'),
[Input('total-data', 'value')])
def plot_totaldata(value):
    time.sleep(1)
    if value == 'clade':
        return dcc.Graph(figure=fig_clade_asia)
    elif value == 'who':
        return dcc.Graph(figure=fig_who_asia)
    elif value == 'number':
        return dcc.Graph(figure=fig_sum_asia)

@app.callback(Output('total-data-plot-Europe', 'children'),
[Input('total-data-Europe', 'value')])
def plot_totaldata(value):
    time.sleep(1)
    if value == 'clade':
        return dcc.Graph(figure=fig_clade_europe)
    elif value == 'who':
        return dcc.Graph(figure=fig_who_europe)
    elif value == 'number':
        return dcc.Graph(figure=fig_sum_europe)

@app.callback(Output('total-data-plot-North-America', 'children'),
[Input('total-data-North-America', 'value')])
def plot_totaldata(value):
    time.sleep(1)
    if value == 'clade':
        return dcc.Graph(figure=fig_clade_north_america)
    elif value == 'who':
        return dcc.Graph(figure=fig_who_north_america)
    elif value == 'number':
        return dcc.Graph(figure=fig_sum_north_america)

@app.callback(Output('total-data-plot-South-America', 'children'),
[Input('total-data-South-America', 'value')])
def plot_totaldata(value):
    time.sleep(1)
    if value == 'clade':
        return dcc.Graph(figure=fig_clade_south_america)    
    elif value == 'who':
        return dcc.Graph(figure=fig_who_south_america)
    elif value == 'number':
        return dcc.Graph(figure=fig_sum_south_america)

@app.callback(Output('total-data-plot-Africa', 'children'),
[Input('total-data-Africa', 'value')])
def plot_totaldata(value):
    time.sleep(1)
    if value == 'clade':
        return dcc.Graph(figure=fig_clade_africa)
    elif value == 'who':
        return dcc.Graph(figure=fig_who_africa)
    elif value == 'number':
        return dcc.Graph(figure=fig_sum_africa)
    
@app.callback(Output('total-data-plot-Oceania', 'children'),
[Input('total-data-Oceania', 'value')])
def plot_totaldata(value):
    time.sleep(1)
    if value == 'clade':
        return dcc.Graph(figure=fig_clade_oceania)
    elif value == 'who':
        return dcc.Graph(figure=fig_who_oceania)
    elif value == 'number':
        return dcc.Graph(figure=fig_sum_oceania)


# # Singapore

# In[21]:


sin_layout = html.Div([
                        html.H3('Singapore trends (COVID-19)'),
                        html.P(['''Overall statistics of COVID-19 cases in Singapore. A time-series plot shows the percentage and number of cases for each clade or lineage (7-day average)''']),
                        html.Div(children=[
                                            dcc.Loading(
                                                id="loading-plot_sin",
                                                type="default",
                                                children=html.Div(id="total-data-plot-sin")
                                                       )
                                      ]),
                        html.Div(children=[
                                            dcc.Loading(
                                                id="who-plot_sin",
                                                type="default",
                                                children=html.Div(id="who-data-plot-sin")
                                                       )
                                      ]),
#                         html.Div(id='total-data-plot-thai') ,
                        html.Div(children=[dcc.RadioItems(id='total-data-sin',
                                    options=[{'label': 'WHO naming (%)', 'value': 'who'},
                                             {'label': 'GISAID clade (%)', 'value': 'clade'},
                                             {'label': 'Number of cases', 'value': 'number'},],
                                    value='who',
                                    labelStyle={'display': 'inline-block'}
                                    ) 
                                 ]),
                        
                        html.P(['''Top 5 PANGO lineages in last 30 days = ''',', '.join(country_lineage['Singapore'])]),
                        html.Hr(),
    
                        html.H6(['''Diversity at each nucleotide position''']),
                        
                        dcc.Graph(figure=fig_enp_sin), 
    
                        html.Hr(),
                        html.H6(['''Diversity of SARS-CoV-2 genomes''']),
                        html.P(['''Multiple SARS-CoV-2 genomes are clustered based on top 50 positions with highest diversity and mutations predefined in GISAID. Each dot represents a genome and colors represent different lineages (Left) and k-mean clustering results (Right). ''',
                                ]),

                        html.Div([
                            html.Div([
                                html.Img(src='data:image/png;base64,{}'.format(call_image('../download/fig_cluster_who_sin.png')), style={'width':'86%'})
#                                 dcc.Graph(id='g1', figure=fig_cluster_result_thai,)
                            ], className="six columns"), 
                            html.Div([
                                html.Img(src='data:image/png;base64,{}'.format(call_image('../download/fig_cluster_k_mean_sin.png')), style={'width':'86%'})
#                                 dcc.Graph(id='g2', figure=fig_cluster_k_mean_thai)
                            ], className="six columns"),
                        ], className="row"),
                        html.Hr(),
                         html.H6('Dendrogram of SARS-CoV-2 genomes'),

                        dcc.Loading(type="default",
                                     children= html.Div([ 
                                                 dcc.Graph(figure=fig_dendro_sin), 
#                                                  dcc.Graph(figure=fig_dendro_sin_kmean),
                                     ]
                                                       )),
                        html.Hr(),
                        
                        html.P(['''Last update : 30 June 2021''']),
                        html.Br(), 
                        html.Br(),
    
    ])



@app.callback(Output('total-data-plot-sin', 'children'),
[Input('total-data-sin', 'value')])
def plot_totaldata_thai(value):
    time.sleep(1)
    if value == 'clade':
        return dcc.Graph(figure=fig_clade_sin)
    elif value == 'who':
        return dcc.Graph(figure=fig_who_sin)
    elif value == 'number':
        return dcc.Graph(figure=fig_sum_sin)


# # Build Tree

# In[22]:


pd.options.display.float_format = '${:.2f}'.format


# In[23]:


build_layout = html.Div([
                    html.Br(),
                    html.H3(['''Upload your SARS-CoV-2 genomes''']),
                    html.P(['''Download example files''']),
                    html.A(children='S15_03-03-2020_C.fasta', href="../download/{}".format(urlquote('S15_03-03-2020_C.fasta'))),html.Br(),
                    html.A(children='S16_04-04-2020_D.fasta', href="../download/{}".format(urlquote('S16_04-04-2020_D.fasta'))),html.Br(),
                    html.A(children='S21_17-03-2020_H.fasta', href="../download/{}".format(urlquote('S21_17-03-2020_H.fasta'))),html.Br(),
                    html.A(children='S28_13-03-2020_H.fasta', href="../download/{}".format(urlquote('S28_13-03-2020_H.fasta'))),html.Br(),
                    html.A(children='S29_14-03-2020_C.fasta', href="../download/{}".format(urlquote('S29_14-03-2020_C.fasta'))),html.Br(),
    
                    dcc.Upload(
                        id='upload_data',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Files')
                        ]),
                        style={
                            'width': '90%',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '10px'
                        },
                        # Allow multiple files to be uploaded
                        multiple=True
                    ),
    
                    dcc.Loading(type="default",
                                children=html.Div(id="file-list")
                               ),
    
                    html.Hr(),
                    html.H3('Lineage statistics (last update : 16 July 2021)'),                        

                    html.P(['''Lineage  ''',
                            dcc.Input(id="input1-lin", type="text", value = 'B.1.1.7', size='8', style={"margin-left": "10px"}),
                            dcc.Checklist(id="check-sub-lin",
                                        options=[
                                            {'label': 'Include sublineage', 'value': 'yes'}, 
                                        ], value=[]) ]), 
                    html.Br(),
                    html.P(['''Month  ''',      
                            dcc.Input(id="input2-lin-month", type="text", value = '2021-02', size='8', style={"margin-left": "10px"})]),
                    html.P(['''Top  ''',      
                            dcc.Input(id="input3-lin-top", type="text", value = 10, size='8', style={"margin-left": "10px"})]),                    
                    html.Button('Submit', id='submit-val-lin', n_clicks=0, style={"margin-left": "20px"}),
                    html.Br(),html.Br(),

                    html.Div(id='stat-lin'),
                    html.Br(), html.Br(),
])


def save_file(name, content, path):
    """Decode and store a file uploaded with Plotly Dash."""
    t = time.localtime()
    timestamp = time.strftime('%b-%d-%Y_%H%M%S', t)
    new_NAME = timestamp + '_' + name
    
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(path, new_NAME), "wb") as fp:
        fp.write(base64.decodebytes(data))    
    return new_NAME


@app.callback(Output('file-list', 'children'),
              [Input('upload_data', 'filename'),
              Input('upload_data', 'contents')])
def update_output(uploaded_filenames, uploaded_file_contents):
    list_file = []
    file_min_length = []
    
    
    if uploaded_filenames is not None and uploaded_file_contents is not None:
        #copy file reference + align
        t = time.localtime()
        timestamp = time.strftime('%b-%d-%Y_%H%M%S', t)
        path = ('/home/ubuntu/COVID/upload/'+ timestamp +'/')
        os.makedirs(path)
        
        ref_name = 'Reference_30-12-2019_Reference.fasta'
        original = '/home/ubuntu/COVID/COVID_data/GCF_009858895.2_ASM985889v3_genomic.fasta'
        target = (path + ref_name)
        shutil.copyfile(original, target)
        
        
        lin_sample, clade_test, change_protein, align_array, snps_max, len_genome, snps_insertion, choose_pos, start_end = align_array_test(ref_name,path)

        X_id_list = [ref_name.split('_')[0]]
        clade_dendro = pd.DataFrame(X_id_list,columns=['id'])
        clade_dendro['result'] = clade_test
        clade_dendro['lineage'] = lin_sample
        clade_dendro['day'] = ref_name.split('_')[1].split('-')[0]
        clade_dendro['month'] = ref_name.split('_')[1].split('-')[1]
        clade_dendro['year'] = ref_name.split('_')[1].split('-')[2]
        clade_dendro['group'] = ref_name.split('_')[2].split('.')[0]
        clade_dendro['coverage'] = round(len_genome, 2)

        X_df = align_array 
        snps_max_df = snps_max            
        change_protein_df = change_protein
        clade_dendro_df = clade_dendro
        snps_insertion_df = snps_insertion
        choose_pos_df = choose_pos
        start_stop_df = start_end
        
        #align file upload
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            if '.fasta' in name:
                new_name = save_file(name, data, path)                
                lin_sample,clade_test, change_protein, align_array, snps_max, len_genome, snps_insertion, choose_pos, start_end = align_array_test(new_name,path)
#                 if len_genome>=95:
                list_file += [new_name]

                X_id_list = [new_name.split('_')[2]]

                clade_dendro = pd.DataFrame(X_id_list,columns=['id'])
                clade_dendro['result'] = clade_test
                clade_dendro['lineage'] = lin_sample
                clade_dendro['day'] = new_name.split('_')[3].split('-')[0]
                clade_dendro['month'] = new_name.split('_')[3].split('-')[1]
                clade_dendro['year'] = new_name.split('_')[3].split('-')[2]
                clade_dendro['group'] = new_name.split('_')[4].split('.')[0]
                clade_dendro['coverage'] = round(len_genome, 2)

                X_df = pd.concat([X_df, align_array], axis=0)
                snps_max_df = pd.concat([snps_max_df, snps_max], axis=0)
                change_protein_df = pd.concat([change_protein_df, change_protein], axis=0)
                clade_dendro_df = pd.concat([clade_dendro_df, clade_dendro], axis=0)
                snps_insertion_df = pd.concat([snps_insertion_df, snps_insertion], axis=0)
                choose_pos_df = pd.concat([choose_pos_df, choose_pos], axis=0)
                start_stop_df = pd.concat([start_stop_df, start_end], axis=0)
                
                if len_genome<95:
                    file_min_length += ['_'.join(new_name.split('_')[2:])]
                    print(file_min_length,len_genome)
        
        if len(file_min_length)==0:
            min_length = '-'
        else:
            min_length = ', '.join([str(elem) for elem in file_min_length])
                   
            
        if len(list_file) == 0:
            return  html.P(['''Genomes with coverage less than 95% of the reference : ''', min_length]) 
        else:
            insertion_sum_df = insertion_summary(snps_insertion_df,choose_pos_df)
            insertion_sum_df['n_insert_base'] = insertion_sum_df['n_insert_base'].astype(int)
#             insertion_sum_df.to_csv('../upload/insertion_sum_df.csv')
    #         insertion_sum_df = insertion_sum_df['position_ref'].astype(int)

            fig_dendro_built, fig_heatmap, result_upload_df, download_NAME = prepare_plot_dendro(X_df, clade_dendro_df, change_protein_df, snps_max_df, start_stop_df, insertion_sum_df)
#             fig_dendro_built.write_image("../pic/fig_dendro_built.png", scale=3)
            data = result_upload_df.to_dict('rows')
            columns=[
                    {"name": "Id", "id": "Id"},
                    {"name": "Lineage", "id": "Lineage"},
                    {"name": "Clade", "id": "Clade"},
                    {"name": "SNPs", "id": "SNPs"},
                    {"name": "Deletion", "id": "Deletion"},
                    {"name": "Insertion", "id": "Insertion"},
                    {"name": "Amino acid", "id": "Amino acid"},
                    ]
            return  html.Div([
                        html.P(['''Genomes with coverage less than 95% of the reference : ''', min_length]) ,
                    html.Br() ,
                    html.P(['''PANGO lineage and GISAID clade and mutation information (Pango Version : ''',pango_version,''')''']) ,
                    dcc.Graph(figure=fig_dendro_built),
                
                    html.P(['''Distances between genomes''']) ,
                    
                    dcc.Graph(figure=fig_heatmap),
                    html.Li(file_download_link(download_NAME)),
                     
            
                    dt.DataTable(data=data, columns=columns,page_size=5,
                            style_cell={'padding': '5px',
                                        'whiteSpace': 'pre-line',
                                        'height': 'auto',},
                            style_header={
                                        'textAlign': 'center'
                                    },
                            style_cell_conditional=[
                                {'if': {'column_id': 'Id'},
                                             'width': '75px','textAlign': 'center'},
                                {'if': {'column_id': 'Lineage'},
                                             'width': '75px'},
                                {'if': {'column_id': 'Clade'},
                                             'width': '75px','textAlign': 'center'},
                                {'if': {'column_id': 'SNPs'},
                                             'width': '180px','textAlign': 'left'},
                                {'if': {'column_id': 'Deletion'},
                                             'width': '70px','textAlign': 'left'},
                                {'if': {'column_id': 'Insertion'},
                                             'width': '75px','textAlign': 'left'},
                                {'if': {'column_id': 'Amino acid'},
                                             'width': '180px','textAlign': 'left'},
                            ],
                    
                            fixed_columns={'headers': True, 'data': 1},
                            style_table={'minWidth': '90%',})
                    ])


@app.callback(Output('stat-lin', 'children'),
[Input('submit-val-lin', 'n_clicks'),
 State('input1-lin', 'value'),
 State('input2-lin-month', 'value'),
 State('input3-lin-top', 'value'),
 State('check-sub-lin', 'value'),
])
def count_lin(n_clicks,input1,input2,input3,input4):

    if input1!=None and input2!=None:
        lin = str(input1)
        month = str(input2)
        top = int(input3)
        
        count_metadata_df=pd.read_csv('../COVID_data/count_metadata_df.csv')
#         print(check)
    
        month_count_df = count_metadata_df[count_metadata_df['month']==month]

        #suggestion lineage
        n = len(lin.split('.'))
        month_count_df['root-lin']=['.'.join(i.split('.')[:n]) for i in month_count_df['lin']]
        sug_lin=sorted(list(set(month_count_df[month_count_df['root-lin']==lin]['lin'])))

        if 'yes' in input4:
            month_count_df = month_count_df[month_count_df['root-lin']==lin]
            month_count_df = month_count_df.groupby(['country_1','root-lin','sum'])[['percentage','count']].sum()
            month_count_df = month_count_df.reset_index()
            month_count_df= month_count_df.rename(columns={'root-lin':'lin'})
        else:
            month_count_df = month_count_df[month_count_df['lin']==lin]
        
        month_count_df = month_count_df[['country_1', 'percentage','count','sum','lin']].sort_values('percentage',ascending=False).head(top)
        month_count_df.columns=['Country', 'Percentage','Count','Sum','Lineage']
        data = month_count_df.to_dict('rows')
        columns=[
                {"name": "Country", "id": "Country"},
                {"name": "Percentage", "id": "Percentage", "type": "numeric","format": {'specifier': '.2f'}},
                {"name": f"# {lin} genomes", "id": "Count"},
                {"name": "# total genomes", "id": "Sum"},
                {"name": "Lineage", "id": "Lineage"},
                ]

        return html.Div([
                        html.P(['''All sublineages within the selected month : ''', ',  '.join(sug_lin)]),
                        dt.DataTable(data=data, columns=columns,
                                    style_cell={'padding': '5px'},
                                    style_table={
        #                                     'maxHeight': '50ex',
                                            'width': '90%',
        #                                     'minWidth': '500px',
                                        },

                                        style_cell_conditional=[
                                            {'if': {'column_id': 'Country'},
                                             'width': '150px'},
                                            {'if': {'column_id': 'Percentage'},
                                             'width': '75px'},
                                            {'if': {'column_id': 'Count'},
                                             'width': '75px'},
                                            {'if': {'column_id': 'Sum'},
                                             'width': '75px'},
                                            {'if': {'column_id': 'Lineage'},
                                             'width': '75px'},
                                            {'if': {'column_id': 'Country'},
                                             'textAlign': 'left'
                                            },
                                        ])
                        ])
    
    
    else: 
        return html.Div([
                        html.H1('error'),
                    ])


# In[24]:


def insertion_summary(snps_insertion_df,choose_pos_df):
        id_insertion = list(set(snps_insertion_df['id']))
        len(id_insertion)
        #หาตำแหน่ง insertion ทั้งหมดพร้อมนับจำนวน
        insertion_list = built.find_insertion_position(id_insertion,snps_insertion_df)
        insertion_df = pd.DataFrame(insertion_list, columns=['id','position_ref','n_insert_base','start_end'])
        insertion_df = insertion_df.sort_values('id')
        
        #เรียงตำแหน่งใหม่่เทียบกับ ref
        insertion_sum_df = pd.merge(insertion_df, choose_pos_df, how="right", on=["id", "start_end"])
        insertion_sum_df = insertion_sum_df[insertion_sum_df['position_ref'].notna()]
        insertion_sum_df = insertion_sum_df[(insertion_sum_df['position_ref']>=insertion_sum_df['start_new']) & (insertion_sum_df['position_ref']<=insertion_sum_df['end_new'])]
        return insertion_sum_df


# In[25]:


def prepare_plot_dendro(X_df, clade_dendro_df, change_protein_df, snps_max_df,start_stop_df,insertion_sum_df):
    id_list = sorted(list(set(start_stop_df.index)))
    
    #no insertion
    if insertion_sum_df.shape[0] == 0:
#         print('no insertion')
        X_id_list = list(clade_dendro_df['id'])
        X_df.columns = X_df.columns.map(str)
        X_df.index = X_id_list
        X_df = X_df.sort_index()
        X_df = X_df.replace('','N')
        X = X_df.to_numpy()

    else : 
#         print('1')
        max_insertion = insertion_sum_df.groupby(['position_ref']).max()
        max_insertion = max_insertion.reset_index()
        max_insertion['sum'] = np.cumsum(max_insertion['n_insert_base'])
        max_insertion = max_insertion[['position_ref','sum']]
        n_insertion = list(max_insertion['sum'])[-1]

        #สร้าง array ขนาด sample x insertion จากนั้นใส่่ 1 ในตำแหน่่งที่มี
        insertion_array = np.zeros((len(id_list),n_insertion), dtype=str)                
        s = 0
        for p in list(max_insertion['position_ref']):
            temp_df = insertion_sum_df[insertion_sum_df['position_ref']==p]
            list_temp = list(temp_df['id'])

            for sample in list_temp:
                start = s
                num = int(temp_df[temp_df['id']==sample]['n_insert_base'])
                i = id_list.index(sample)
                insertion_array[i][start:start+num] = 1

            s = int(max_insertion[max_insertion['position_ref']==p]['sum'].values)

        insertion_array[insertion_array==''] = '0'

        X_id_list = list(clade_dendro_df['id'])
        X_df.columns = X_df.columns.map(str)
        X_df.index = X_id_list
        X_df = X_df.sort_index()
        X_df = X_df.replace('','N')
        X = X_df.to_numpy()
        X = np.concatenate((X, insertion_array),axis=1) #รวมกับ align


    X = X_df.to_numpy()
    X = np.where(X=='A', 1, X)
    X = np.where(X=='T', 2, X)
    X = np.where(X=='G', 3, X)
    X = np.where(X=='C', 4, X)
    X = np.where(X=='N', 0, X)
    X = np.where(X=='-', 5, X)
    X = np.where(X=='Y', 0, X)
    X = np.where(X=='R', 0, X)
    X = np.where(X=='S', 0, X)
    X = np.where(X=='K', 0, X)
    X = np.where(X=='M', 0, X)
    X = np.where(X=='W', 0, X)
    
    clade_dendro_df = clade_dendro_df.sort_values('id')
    clade_dendro_df = clade_dendro_df.reset_index()

    
    X_id_list = sorted(X_id_list)
    
    #for save data
    result_upload = []
    for sample in X_id_list:
        df = snps_max_df[snps_max_df['id']==sample].sort_values('position')
        deletion = ', '.join([str(elem) for elem in list(df[df['sbjct']=='-']['change'])])
        snps = ', '.join([str(elem) for elem in list(df[df['sbjct']!='-']['change'])])
        
        pro_df = change_protein_df[change_protein_df['id']==sample].sort_values('position_protein')
        gene = list(set(pro_df['gene']))
        protein_temp = []
        frame_shift_temp = []
        for g in sorted(gene):
#             list_gene = list(set(pro_df[pro_df['gene']==g]['change_protein']))
            list_gene = list(pro_df[pro_df['gene']==g]['change_protein'].unique())
            protein_temp += [str(g) +': '+ ', '.join(list_gene)]
            for s in list_gene:
                if 'X' in s:
                    frame_shift_temp += [g]
                    
        frame_shift_temp = sorted(list(set(frame_shift_temp)))
        frame_shift = ', '.join(str(fs) for fs in frame_shift_temp)
        
        protein = '\n'.join(str(elem) for elem in protein_temp)
        
        clade = clade_dendro_df[clade_dendro_df['id']==sample]['result'].values[0]
        lineage = clade_dendro_df[clade_dendro_df['id']==sample]['lineage'].values[0]
        coverage = clade_dendro_df[clade_dendro_df['id']==sample]['coverage'].values[0]

        if list(insertion_sum_df[insertion_sum_df['id']==sample]['position_ref'])==[]:
            insert_temp = []
            insertion = ', '.join([str(elem) for elem in insert_temp])
        else:
            in_df = insertion_sum_df[insertion_sum_df['id']==sample]
            insert_temp = tuple(zip(in_df['position_ref'].astype(int),in_df['n_insert_base'].astype(int)))
            insertion = ', '.join([str(elem) for elem in insert_temp])

        keys = ['id','lineage','clade','snps','deletion','insertion','protein','frame shift','coverage']
        values = [sample,lineage,clade,snps,deletion,insertion,protein,frame_shift,coverage]
        result_upload += [values]
    

    result_upload_df = pd.DataFrame(result_upload, columns=keys)
    result_upload_df.columns = ['Id', 'Lineage', 'Clade','SNPs','Deletion','Insertion','Amino acid','Frame shift','Coverage']
    result_upload_df.set_index('Id')
    

    t = time.localtime()
    timestamp = time.strftime('%m-%d-%Y_%H%M%S', t)
    download_NAME = ('../upload/RESULT_' + timestamp + '.xlsx')
    
    
 

    writer = pd.ExcelWriter(download_NAME, engine='xlsxwriter')
    result_upload_df.to_excel(writer, sheet_name='Sheet1', index=False)
    workbook=writer.book
    worksheet = writer.sheets['Sheet1']

    format = workbook.add_format({'text_wrap': True, 'valign': 'top'})
    for col_num, value in enumerate(result_upload_df.columns.values):
            worksheet.write(0, col_num, value, format)

    for i in range(len(result_upload_df.index)):
        for col_num, value in enumerate(result_upload_df.loc[result_upload_df.index[i],:].values):
            worksheet.write(i+1, col_num, value, format)

    
    fig,fig_heatmap,position_dendro_df,dict_color_label,heatmap_df = built.fig_dendro_built_tree(X, X_id_list, clade_dendro_df, snps_max_df, change_protein_df)
    
    heatmap_df.to_excel(writer, sheet_name='Sheet2')
    writer.save()
    
    return fig, fig_heatmap, result_upload_df, download_NAME


# In[ ]:





# In[ ]:





# In[ ]:





# # Alarm

# In[26]:


gene_df = pd.read_excel('../preprocessed/gene.xlsx')
gene_df['position'] = [str(s)+' - '+str(e) for s,e in zip(gene_df['start'],gene_df['end'])]
gene_df = gene_df[['gene','position']]
data = gene_df.to_dict('rows')
columns=[
        {"name": "Gene", "id": "gene"},
        {"name": "Position", "id": "position"},
                ]


alarm_layout = html.Div([                                                   
                        html.Div(className='row',  # Define the row element
                               children=[
                                   html.H3('Inspect region of interest in SARS-CoV-2 genomes'), 
                                   html.P(['''Statistics are based on 3,406 genomes in Singapore (downloaded from GISAID as of 30 June 2021).''']),
                                  html.Div(className='three columns div-user-controls',
                                          children = [
                                                 dt.DataTable(data=data[0:9], columns=columns,
                                                    style_cell={'padding': '5px'},
                                                    style_table={
                                                            'width': '100px','textAlign': 'center'},
                                                    style_cell_conditional=[
                                                            {'if': {'column_id': 'Gene'},
                                                             'width': '50px'},
                                                    ]) 
                                          ]),  # Define the left element
#                                   
                                   html.Div(className='three columns div-user-controls',
                                          children = [
                                                 dt.DataTable(data=data[9:18], columns=columns,
                                                    style_cell={'padding': '5px'},
                                                    style_table={
                                                            'width': '100px','textAlign': 'center'},
                                                    style_cell_conditional=[
                                                            {'if': {'column_id': 'Gene'},
                                                             'width': '50px'},
                                                    ]) 
                                          ]), # Define the middle element
                                   
                                   html.Div(className='three columns div-user-controls',
                                          children = [
                                                 dt.DataTable(data=data[18:], columns=columns,
                                                    style_cell={'padding': '5px'},
                                                    style_table={
                                                            'width': '100px','textAlign': 'center'},
                                                    style_cell_conditional=[
                                                            {'if': {'column_id': 'Gene'},
                                                             'width': '50px'},
                                                    ]) 
                                          ])
                                  ]), # Define the right element
                        
    

                        html.Br(),html.Br(),
                        html.P(['''Nucleotide position start ''',
                                dcc.Input(id="input1", type="text", value = 23400, size='8', style={"margin-left": "10px","margin-right": "10px"}),
                                ''' end ''',
                                dcc.Input(id="input2", type="text", value = 23450, size='8', style={"margin-left": "10px"}),
                                html.Button('Submit', id='submit-val', n_clicks=0, style={"margin-left": "20px"})]),
    
                        html.Div(id='alarm-plot-neu'),
                        
                        
                        html.Hr(),    
                        html.H3('Amino acid changes in a specific protein'),  
                        html.Div([
                                dcc.Dropdown(
                                    id='choose-protein',
                                    options=[
                                        {'label': ' S (1 - 1273)', 'value': 'S'},
                                        {'label': ' E (1 - 75)', 'value': 'E'},
                                        {'label': ' M (1 - 222)', 'value': 'M'},
                                        {'label': ' N (1 - 419)', 'value': 'N'},
                                        {'label': ' ORF3a (1 - 275)', 'value': 'ORF3a'},
                                        {'label': ' ORF6 (1 - 61)', 'value': 'ORF6'},
                                        {'label': ' ORF7a (1 - 121)', 'value': 'ORF7a'},
                                        {'label': ' ORF7b (1 - 43)', 'value': 'ORF7b'},
                                        {'label': ' ORF8 (1 - 121)', 'value': 'ORF8'},
                                        {'label': ' ORF10 (1 - 38)', 'value': 'ORF10'},
                                        {'label': ' NSP1 (1 - 180)', 'value': 'NSP1'},
                                        {'label': ' NSP2 (1 - 638)', 'value': 'NSP2'},
                                        {'label': ' NSP3 (1 - 1944)', 'value': 'NSP3'},
                                        {'label': ' NSP4 (1 - 499)', 'value': 'NSP4'},
                                        {'label': ' NSP5 (1 - 605)', 'value': 'NSP5'},
                                        {'label': ' NSP6 (1 - 289)', 'value': 'NSP6'},
                                        {'label': ' NSP7 (1 - 82)', 'value': 'NSP7'},
                                        {'label': ' NSP8 (1 - 197)', 'value': 'NSP8'},
                                        {'label': ' NSP9 (1 - 112)', 'value': 'NSP9'},
                                        {'label': ' NSP10 (1 - 138)', 'value': 'NSP10'},
                                        {'label': ' NSP12 (1 - 931)', 'value': 'NSP12'},
                                        {'label': ' NSP13 (1 - 600)', 'value': 'NSP13'},
                                        {'label': ' NSP14 (1 - 526)', 'value': 'NSP14'},
                                        {'label': ' NSP15 (1 - 345)', 'value': 'NSP15'},
                                        {'label': ' NSP16 (1 - 297)', 'value': 'NSP16'},                         
                                   ], style = dict(width = '60%'),
                                    value='S',
                                    clearable=False
                                ),
                            ]),
                        html.Br(),
                        html.P(['''Amino acid position start ''',
                                dcc.Input(id="input1-aa", type="text", value = 500, size='8', style={"margin-left": "10px","margin-right": "10px"}),
                                ''' end ''',
                                dcc.Input(id="input2-aa", type="text", value = 550, size='8',style={"margin-left": "10px"}),
                                html.Button('Submit', id='submit-val-aa', n_clicks=0, style={"margin-left": "20px"})]),
    
                        html.Div(id='alarm-plot-aa'),
                        html.Br(),html.Hr(),
    
                        html.H3('Inspect amino acid changes over time'),
    
                        html.P(['''Upload Excel file containing a list of amino acid changes ''', html.A(children='(example file)', href="../download/{}".format(urlquote('test_amino_acid_change.xlsx'))), ''', where the first columns is a protein name (i.e. S, E, N) and the second columns contains amino acid change (i.e. D614G, N501Y)''']),
                            dcc.Upload(
                                id='upload-alarm-excel',
                                children=html.Div([
                                    'Drag and Drop or ',
                                    html.A('Select File')
                                ]),
                                style={
                                    'width': '90%',
                                    'height': '60px',
                                    'lineHeight': '60px',
                                    'borderWidth': '1px',
                                    'borderStyle': 'dashed',
                                    'borderRadius': '5px',
                                    'textAlign': 'center',
                                    'margin': '10px'
                                },
                                # Allow multiple files to be uploaded
                                multiple=False
                            ),
                        html.Div(id='alarm-plot-excel'),
                        html.Br(),html.Hr(),

    
                        html.H3('Inspect co-occurrence amino acid changes over time'),
                        html.P(['''Note : ORF1ab = NSP1-NSP16 ''',
                              html.A(children='(more detail)', href='https://www.ncbi.nlm.nih.gov/sars-cov-2/', target='_blank', )]),
                            
                        html.Div(className='row',  # Define the row element
                               children=[
                                  html.Div(className='three columns div-user-controls',
                                          children = [
                                              html.P([''' S :  ''',
                                                    dcc.Input(id="input-S", type="text", value = 'N501Y, E484K', size='14')]),
                                              html.P([''' E :  ''',
                                                    dcc.Input(id="input-E", type="text", value = 'P71L', size='14')]),
                                              html.P([''' M :  ''',
                                                    dcc.Input(id="input-M", type="text", size='14')]),
                                              html.P([''' N :  ''',
                                                    dcc.Input(id="input-N", type="text", size='14')]),
                                              html.P(['''ORF3a : ''',
                                                    dcc.Input(id="input-ORF3a", type="text", size='14')]),
                                              html.P(['''ORF6 : ''',
                                                    dcc.Input(id="input-ORF6", type="text", size='14')]),
                                              html.P(['''ORF7a : ''',
                                                    dcc.Input(id="input-ORF7a", type="text", size='14')]),
                                              html.P(['''ORF7b : ''',
                                                    dcc.Input(id="input-ORF7b", type="text", size='14')]),
                                              html.P(['''ORF8 : ''',
                                                    dcc.Input(id="input-ORF8", type="text", size='14')]),
                                          
                                            ]),  # Define the left element
                                   
#                                   
                                   html.Div(className='three columns div-user-controls bg-grey',
                                          children=[
                                              
                                              html.P(['''ORF10 : ''',
                                                    dcc.Input(id="input-ORF10", type="text", size='14')]),
                                              html.P(['''NSP1 : ''',
                                                    dcc.Input(id="input-NSP1", type="text", size='16')]),
                                              html.P(['''NSP2 : ''',
                                                    dcc.Input(id="input-NSP2", type="text", size='16')]),
                                              html.P(['''NSP3 : ''',
                                                    dcc.Input(id="input-NSP3", type="text", value = 'K837N', size='16')]),
                                              html.P(['''NSP4 : ''',
                                                    dcc.Input(id="input-NSP4", type="text", size='16')]),
                                              html.P(['''NSP5 : ''',
                                                    dcc.Input(id="input-NSP5", type="text", size='16')]),
                                              html.P(['''NSP6 : ''',
                                                    dcc.Input(id="input-NSP6", type="text", size='16')]),
                                              html.P(['''NSP7 : ''',
                                                    dcc.Input(id="input-NSP7", type="text", size='16')]),
                                              html.P(['''NSP8 : ''',
                                                    dcc.Input(id="input-NSP8", type="text", size='16')]),
                                           
                                          ]),  # Define the middle element
#                                   
                                   html.Div(className='three columns div-user-controls bg-grey',
                                          children=[
                                              html.P(['''NSP9 : ''',
                                                    dcc.Input(id="input-NSP9", type="text", size='16')]),
                                              html.P(['''NSP10 : ''',
                                                    dcc.Input(id="input-NSP10", type="text", size='14')]),
                                              html.P(['''NSP12 : ''',
                                                    dcc.Input(id="input-NSP12", type="text", size='14')]),
                                              html.P(['''NSP13 : ''',
                                                    dcc.Input(id="input-NSP13", type="text", size='14')]),
                                              html.P(['''NSP14 : ''',
                                                    dcc.Input(id="input-NSP14", type="text", size='14')]),
                                              html.P(['''NSP15 : ''',
                                                    dcc.Input(id="input-NSP15", type="text", size='14')]),
                                              html.P(['''NSP16 : ''',
                                                    dcc.Input(id="input-NSP16", type="text", size='14')]),
                                              html.Button('Submit', id='submit-val-protein', n_clicks=0, style={"margin-left": "20px"}),
                                          ]) # Define the right element
                                  ]),
    
                        html.Div(id='alarm-plot-specific'),
                        html.Br(),html.Br(),html.Br(),
])


@app.callback(Output('alarm-plot-neu', 'children'),
[Input('submit-val', 'n_clicks'),State('input1', 'value'),State('input2', 'value')])
def plot_alarm_neu(n_clicks,input1,input2):

    if input1!=None and input2!=None:
        min_range = int(input1)
        max_range = int(input2)
        position = range(min_range,max_range+1)
        neucleotide_alarm_df = pd.read_csv('../preprocessed/{}_neucleotide_alarm_df.csv'.format('singapore'))
        del_alarm_df = pd.read_csv('../preprocessed/{}_del_alarm_df.csv'.format('singapore'))
        insert_alarm_df = pd.read_csv('../preprocessed/{}_insert_alarm_df.csv'.format('singapore'))
        return html.Div([
                html.H3('''SNPs'''),
                dcc.Graph(figure=plot_data.alarm_plot(neucleotide_alarm_df,position,'position')),
                html.Br(),html.H3('''Deletion'''),
                dcc.Graph(figure=plot_data.alarm_plot2(del_alarm_df,position,'position')),
                html.Br(),html.H3('''Insertion'''),
                dcc.Graph(figure=plot_data.alarm_plot2(insert_alarm_df,position,'position_ref')),
        ])
            
    else: 
        return html.Div([
    html.H1('error'),
])

                         
@app.callback(Output('alarm-plot-aa', 'children'),
[Input('submit-val-aa', 'n_clicks'),
 State('input1-aa', 'value'),
 State('input2-aa', 'value'),
 State('choose-protein', 'value')])
def plot_alarm_aa(n_clicks,input1,input2,protein):

    if input1!=None and input2!=None:
        min_range = int(input1)
        max_range = int(input2)
        position = range(min_range,max_range+1)
        protein_alarm = pickle.load(open('../preprocessed/{}_protein_alarm.pickle'.format('singapore'), "rb"))
        protein_plot = protein_alarm[str(protein)] #choose protein
        return dcc.Graph(figure=plot_data.alarm_plot(protein_plot,position,'position_protein'))
    else: 
        return html.Div([
    html.H1('error'),
])


# In[ ]:





# In[ ]:





# In[27]:


@app.callback(Output('alarm-plot-excel', 'children'),
              [Input('upload-alarm-excel', 'filename'),
              Input('upload-alarm-excel', 'contents')])
def update_output(uploaded_filenames, uploaded_file_contents):
    if uploaded_filenames is not None and uploaded_file_contents is not None:                
        #file upload
        path = '/home/ubuntu/COVID/upload/'

        if '.xls' in uploaded_filenames:
            new_name = save_file(uploaded_filenames, uploaded_file_contents, path)                
            gene_and_position = pd.read_excel(path + new_name)

            gene_list = list(gene_and_position[gene_and_position.columns[0]])
            pos_list = list(gene_and_position[gene_and_position.columns[1]])
            
            return dcc.Graph(figure=plot_data.alarm_excel_plot(gene_list,pos_list))
        else: 
            return html.Div([
                        html.H1('Please upload Excel file'),
                    ])


@app.callback(Output('alarm-plot-specific', 'children'),
[Input('submit-val-protein', 'n_clicks'),State('input-NSP1', 'value'),State('input-NSP2', 'value'),
 State('input-NSP3', 'value'),State('input-NSP4', 'value'),State('input-NSP5', 'value'),
 State('input-NSP6', 'value'),State('input-NSP7', 'value'),State('input-NSP8', 'value'),
 State('input-NSP9', 'value'),State('input-NSP10', 'value'),State('input-NSP12', 'value'),
 State('input-NSP13', 'value'),State('input-NSP14', 'value'),State('input-NSP15', 'value'),
 State('input-NSP16', 'value'),State('input-S', 'value'),State('input-ORF3a', 'value'),
 State('input-E', 'value'),State('input-M', 'value'),State('input-ORF6', 'value'),
 State('input-ORF7a', 'value'),State('input-ORF7b', 'value'),State('input-ORF8', 'value'),
 State('input-N', 'value'),State('input-ORF10', 'value')])
def plot_alarm_specific_protein(n_clicks,input1,input2,input3,input4,input5,input6,input7,input8,input9,input10,input11,input12,input13,input14,input15,input16,input17,input18,input19,input20,input21,input22,input23,input24,input25):
    gene_list = ['NSP1','NSP2','NSP3','NSP4','NSP5','NSP6','NSP7','NSP8','NSP9','NSP10','NSP12','NSP13','NSP14','NSP15','NSP16','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10']
    position_list = [input1,input2,input3,input4,input5,input6,input7,input8,input9,input10,input11,input12,input13,input14,input15,input16,input17,input18,input19,input20,input21,input22,input23,input24,input25]
    position_and_gene = []
    for i in range(len(gene_list)):
        if position_list[i]!=None:
            temp_list = [p.strip() for p in str(position_list[i]).split(',')]
            for k in temp_list:
                position_and_gene += [[gene_list[i],k]]

    position_and_gene_df = pd.DataFrame(position_and_gene,columns=['gene','position'])
    position_and_gene_df = position_and_gene_df[position_and_gene_df['position']!='']

    gene_list = list(position_and_gene_df['gene'])
    pos_list = list(position_and_gene_df['position'])
    if (gene_list!=[]) and (pos_list!=[]):
        return dcc.Graph(figure=plot_data.alarm_specific_protein(gene_list,pos_list))
    else:
        html.Div([html.H1(''),
                    ])


# # Define the app

# In[28]:


app.layout = html.Div(children=[
                      html.Div(className='row',  # Define the row element
                               children=[
                                  html.Div(className='two columns div-user-controls',
                                          children = [html.Img(src='data:image/png;base64,{}'.format(call_image('../download/logo_02.png')),style={'width':'100%'}),
                                                    html.Br(),html.Br(),
                                                    html.P([''' Analyse multiple SARS-CoV-2 genomes''',html.Br(),html.Br(),
                                                            ''' Study genomic patterns across different clades ''',html.Br(),html.Br(),
                                                            ''' Inspect possible variants of concern ''']), 

                                          ]),  # Define the left element
#                                   
                                   html.Div(className='ten columns div-user-controls bg-grey',
                                          children=[dcc.Location(id='url', refresh=False),

                                                html.Div(
                                                    id="tabs_holder",
                                                    children=[dcc.Tabs(id="tabs", value='/about')]  
                                                ),

                                                html.Div(id='page-content'),
                                            ]) # Define the right element
                                  ])
                        ])


# In[29]:


@app.callback([Output('page-content', 'children'),
               Output('tabs_holder', 'children')],
              [Input('url', 'pathname')])
def display_page(pathname):
    tabs = [
           dcc.Tabs(
                id="tabs",
                value=pathname,
                children=[
                    dcc.Tab(label='About', value='/about'),
                    dcc.Tab(label='Trends by Continent', value='/continents'),
                    dcc.Tab(label='Singapore Stats', value='/singapore'),
                    dcc.Tab(label='Analyse Genomes', value='/analyse'),
                    dcc.Tab(label='Inspect ROI (Singapore)', value='/inspect'),
                ],style={
                        'width': '100%',
                        'height': '50%'
                         },
                  colors={
                    "border": "white",
                    "primary": "#747879",
                    "background": "#F2F2F3"
                }
            )
        ]

    
    
    if pathname == '/about':
        return about_layout, tabs
    elif pathname == '/continents':
        return region_layout, tabs
    elif pathname == '/singapore':
        return sin_layout, tabs
    elif pathname == '/analyse':
        return build_layout, tabs
    elif pathname == '/inspect':
        return alarm_layout, tabs
    else:
        return html.Div([html.H1('Error 404 - Page not found')]), tabs


@app.callback(Output('url', 'pathname'),
              [Input('tabs', 'value')])
def tab_updates_url(value):
    return value


# #  Test data

# In[30]:


def align_array_test(file_name,upload_path):
    t = time.localtime()
    timestamp = time.strftime('%b-%d-%Y_%H%M%S', t)
    
    ref_path = '/home/ubuntu/COVID/COVID_data/GCF_009858895.2_ASM985889v3_genomic.fasta'
 
   
    BACKUP_NAME = ('alignment_' + timestamp + '.xml')
    Pango_outfile = (file_name.split('.fasta')[0]  + '.csv')
    
    
    command_pango_list = ['pangolin', '--outdir', f'{upload_path}', '--outfile ', f'{upload_path}{Pango_outfile}', f'{upload_path}{file_name}']
    
    pango_output = subprocess.check_output(' '.join(command_pango_list), shell=True)
    pango_df = pd.read_csv(f'{upload_path}{Pango_outfile}')
    lin_sample = pango_df['lineage'].values[0]


    command_list = ['blastn', '-query', ref_path, '-subject', f'{upload_path}{file_name}', '-outfmt', '5', '-out', f'{upload_path}{BACKUP_NAME}']
    proc = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.communicate()
    
    
    snps_result, start_end_result = utils.start_end_list_func_built_tree(f'{upload_path}{BACKUP_NAME}')

    start_end_df = pd.DataFrame(start_end_result, columns=['id', 'start_end', 'length'])    
    id_list = [file_name.split('_')[2]] #เปลี่ยน id ให้ตรงกับชื่อไฟล์
    sample_id = file_name.split('_')[2]
    
    snps_df = pd.DataFrame(snps_result, columns=['id','query','position','sbjct', 'start_end', 'length'])    
    snps_df['id'] = sample_id
    start_end_df['id'] = sample_id

    start_end_df = start_end_df.set_index('id') #for utils.length_max_snps

    
    snps_insertion_df = snps_df.copy()
    snps_insertion_df = snps_insertion_df[(snps_insertion_df['query']=='-')]  
    
    
    #หาตำแหน่งที่จะตัดออก (delete insertion)
    #ถ้ามี insertion ให้ตัดออก
    if snps_df[snps_df['query']=='-']['query'].sum()!=0:
        snps_df = utils.del_insertion(id_list,snps_df)
   
    snps_df['change']=snps_df['query']+snps_df['position'].astype(str)+snps_df['sbjct']

    
    #Load reference
    seq_ref_list = pickle.load(open("../preprocessed/seq_ref_list.pickle", "rb"))
    seq_ref_length = len(seq_ref_list)

    align_array = np.zeros((1,seq_ref_length), dtype=str)

          
    #for insertion ---------------------------------------------------
    snps_max_df = snps_df.copy() 
    sum_length,snps_max_df,align_array,choose_pos_list = utils.length_max_snps2(align_array,id_list,start_end_df,snps_max_df) #หาตำแหน่งของ insertion
    len_genome = int(sum_length[0][1])/seq_ref_length*100
    
    choose_pos_df = pd.DataFrame(choose_pos_list,columns=['id','start_new','end_new','start_end'])
 
    
    #---------------------------------------------------
    
   
    for j in range(seq_ref_length):
        align_array[0][j] = seq_ref_list[j] #สมมติว่าทุกตำแหน่่งมีค่าเหมือน ref


    sample_id_dict = dict(zip(id_list, range(len(id_list))))

    snps_max_df = snps_max_df[snps_max_df['sbjct']!='N'] #ตำแหน่งที่่มี N ให้เป็น ref
    
    # ใส่ subject ในตำแหน่ง snps
    for _,row in snps_max_df.iterrows():
        inx = sample_id_dict[row['id']]
        align_array[inx][row['position']-1] = row['sbjct']


    align_array_df = pd.DataFrame(align_array)
    align_array_df.index = id_list


    #predict clade
    clade_test_df = utils.predict_clade(id_list,snps_max_df,align_array)
    clade_test = clade_test_df['result']



    #check protein   
    position_gene_df = pd.read_csv('../preprocessed/position_gene_df.csv')

 
    snps_max_df = snps_max_df.sort_values('id')

    snps_max_df = snps_max_df.reset_index()

    change_df = snps_max_df.copy()
    change_df = change_df[['id','position']]
    change_df['position']-=1 #ลบเพราะจะเทียบกับ ref ซึ่่งเริ่มที่ 0
    change_df['check']='None'
    
    #สร้าง align array for protein
    align_array_pro,change_df = utils.find_align_array_pro(align_array_df,position_gene_df,change_df)
    change_df = change_df.reset_index()
    change_protein_df = utils.change_protein_table(change_df,position_gene_df,align_array_pro,seq_ref_list)
    change_protein_df = change_protein_df[change_protein_df['check']==False].drop_duplicates(keep='first')
    
    
    #change position & name gene NSP12
    change_protein_df.loc[change_protein_df['gene']=='NSP12_2','position_protein']+=9
    change_protein_df.loc[change_protein_df['gene']=='NSP12_2','gene']= 'NSP12'
    change_protein_df.loc[change_protein_df['gene']=='NSP12_1','gene']= 'NSP12'
    
    change_protein_df['change_protein'] = [str(q)+str(p)+str(s) for q,p,s in zip(change_protein_df['query_protein'],change_protein_df['position_protein'],change_protein_df['sbjct_protein'])]


    return lin_sample,clade_test, change_protein_df, align_array_df, snps_max_df, len_genome, snps_insertion_df, choose_pos_df, start_end_df


# # Run app

# In[31]:


if __name__ == '__main__':
    app.run_server(debug=False, port=8050)







