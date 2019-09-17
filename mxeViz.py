from __future__ import print_function


########################################################################
# File: mxeViz.py
#  executable: mxeViz.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 02/10/2019 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################

import sys, os
import json
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import numpy as np
from textwrap import dedent as d
from dash.dependencies import Input, Output
from plotly import tools


N = 1000
random_x = np.random.randn(N)
random_y = np.random.randn(N)

########################################################################
# Sample
########################################################################


class MXE(object):
    '''
    A simple class to handle sample related infoamtion
    '''
    def __init__(self, intronID=None):
        self.intronID = intronID
        self.chrom = intronID.split(":")[0] 
        self.c1, self.c2 = intronID.split(":")[-1].split("-")
        self.c1, self.c2 = int(self.c1), int(self.c2)
        self.mxes = set()

########################################################################
########################################################################

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

df = pd.read_csv(sys.argv[1], sep="\t", names=["sample","file","group1","group2"])
df2 = pd.read_csv(sys.argv[2], sep="\t", names=["event","iqr","norm median psi","outliers","delta psi"])


styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

# mxes = dict()
# with open(sys.argv[3]) as l:
#     for line in l:
#         c = line.split()
#         i = "%s:%s-%s" % (c[0],c[1],c[2])
#         mxes[i] = MXE(i)

# with open(sys.argv[4]) as l:
#     for line in l:
#         c = line.rstrip().split()
#         intron, mxe  = c
#         allMXE = set([i for i in mxe.split(",")])
#         mxes[intron].mxes = allMXE


fig2 = {
    'data': [{
        'x' : df2["outliers"],
        'y' : df2["delta psi"],
        'type' : 'scatter',
        'text' : df2["event"],
        'mode' : 'markers'#'name': event
        #'orientation':'v'
        }],
    'layout': {
                'clickmode': 'event+select'
            }
}

fig3 = {
    'data': [{
        'x':[0,4,8],
        'y':[1,1,1],
        'text':["exon1","exon2","exon3"],
        'mode':'markers+text',
        'textposition' : 'top center',
        'marker'  : {
                    'symbol': "triangle-down",
                    'size' : 15
                    },
        'type' : 'scatter',
        'name':"test"
        }],
    'layout': {
        'xaxis': {
            'range': [-5, 13],
            'showgrid': False,
            'fixedrange': True,
            'showticklabels':False,
            'zeroline':False,
            'hovermode' : 'closest'
        },
        'yaxis': {
            'range': [-2, 2],
            'showgrid': False,
            'fixedrange': True,
            'showticklabels':False,
            'zeroline':False
        },
        'shapes':[
        # unfilled Rectangle
        {
            'type': 'rect',
            'x0': -2,
            'y0': -1,
            'x1': 2,
            'y1': 0.90,
            'line': {
                'color': 'rgba(128, 0, 128, 1)',
            },
        },
        # filled Rectangle
        {
            'type': 'rect',
            'x0': 3,
            'y0': -1,
            'x1': 5,
            'y1': 0.90,
            'line': {
                'color': 'rgba(128, 0, 128, 1)',
                'width': 2,
            },
            'fillcolor': 'rgba(128, 0, 128, 0.7)',
        },
                            {
            'type': 'rect',
            'x0': 6,
            'y0': -1,
            'x1': 10,
            'y1': 0.90,
            'line': {
                'color': 'rgba(128, 0, 128, 1)',
                'width': 5,
            },
            'fillcolor': 'rgba(128, 0, 128, 0.7)',
        },
    ]
    }
    }




app.layout = html.Div([

   #  html.Div([
   #      html.Img(src="https://brookslabucsc.github.io/mesa/assets/img/logo.png", style={'vertical-align': 'middle','display': 'inline-block', 'width': '15%','height':'15%'}),
   #      html.Strong("\t\tResults Report", style={'vertical-align': 'middle', 'display': 'inline-block'}),
   # ]),

    html.Div([
           html.Div([
            html.Strong("Select sample group(s):"),
                ], style={'width':"33%", 'display': 'inline-block'}),
           html.Div([
            html.Strong("Select sample type(s):"),
                ], style={'width':"33%", 'display': 'inline-block'}),
           ]),

    html.Div([
        html.Div([
            dcc.Dropdown(
                id='group',
                options=[{'label': i, 'value': i} for i in df['group1'].unique()] + [{'label':"All groups","value":"All groups"}],
                value=['All groups'],
                multi=True),
        ],style={'width':"33%", 'display': 'inline-block'}),

        html.Div([
                dcc.Dropdown(
                id='sample_group_type',
                options=[{'label': i, 'value': i} for i in df['group2'].unique()] + [{'label':"All types","value":"All types"}],
                value=['All types'],
                multi=True)

        ], style={'width':"33%",'display': 'inline-block' })
        ]),
        

    html.Div([
        dcc.Graph(id='plot')
    ]),
    
    html.Div([
        html.Strong('Do other stuff here?', id='header-2'),
        dcc.Graph(
                id = "graph2",
                figure = fig2,
                config = {
                            'dubleClick': 'reset'
                }
        )]),
    html.Div([

            dcc.Markdown(d("""
                **Selection Data**

                Choose the lasso or rectangle tool in the graph's menu
                bar and then select points in the graph.

                Note that if `layout.clickmode = 'event+select'`, selection data also 
                accumulates (or un-accumulates) selected data if you hold down the shift
                button while clicking.
            """)),
            html.Div(id="plot2")], className="three columns")
            # html.Div(dcc.Graph(id='plot2', figure={'data': []}),
            #     style={'width': '49%', 'display': 'inline-block', 'vertical-align': 'middle'})
            # ])

    ])

@app.callback(
    dash.dependencies.Output('plot', 'figure'),
    [dash.dependencies.Input('group', 'value'),
    dash.dependencies.Input('sample_group_type', 'value')])

def update_graph(group, sample_group_type):
    

    if "All types" in sample_group_type:
        dff = df
    else:
        dff = df[df['group2'].isin(sample_group_type)]

    if "All groups" in group:
        pass
    else:
        dff = dff[dff['group1'].isin(group)]


    if 'All groups' not in group:

        return  {            
                'data': [go.Bar(
                x       = dff[dff.group2==i].group1.unique(),
                y       = [dff[(dff.group2==i) & (dff.group1==x)].group2.value_counts()[0] for x in dff[dff.group2==i].group1.unique()],
                opacity = 0.7,
                name    = i,
                orientation='v',
                      ) for i in dff.group2.unique()
                    ],
            'layout' : go.Layout(
                hovermode='x',
                margin={'t': 10},
                xaxis= {
                    'fixedrange': True,
                    'hoverformat' : ''
                    
                },
                yaxis= {
                    'fixedrange': True,
                    }
                    )
            
            }

    else:
        return  {            
                'data': [go.Bar(
                x       = dff[dff.group2==i].group1.unique(),
                y       = [dff[(dff.group2==i) & (dff.group1==x)].group2.value_counts()[0] for x in dff[dff.group2==i].group1.unique()],
                opacity = 0.7,
                name    = i,
                orientation='v',
                      ) for i in dff.group2.unique()
                    ],
            'layout' : go.Layout(
                barmode='stack',
                margin={'t': 5},
                hovermode='x',
                xaxis= {
                    'fixedrange': True,
                    'hoverformat' : ''
                },
                yaxis= {
                    'fixedrange': True,
                    }
                    )
            
            }
        

@app.callback(
    Output('plot2', 'children'),
    [Input('graph2', 'selectedData')])
def display_click_data(selectedData):
    plots = list()

    for i,j in selectedData.items():
        for num,k in enumerate(j):

            data = k['text']
            c1,c2 = map(int, data.split(":")[-1].split("-"))
            print(c1,c2)
            plots.append(html.Div([dcc.Graph(
            id='graph-{}'.format(i),
            figure={
                'data': [{
                    'x':[c1+((c2-c1)/2)],
                    'y':[1],
                    'mode':'markers+text',
                    'textposition' : 'top center',
                    'marker'  : {
                                'symbol': "triangle-down",
                                'size' : 15
                                },
                    'type' : 'scatter',
                    
                    }],
                'layout': {
                    'title': 'Graph {}'.format(i),
                    'xaxis': {
                        'range': [c1-100, c2+100],
                        'showgrid': False,
                        'fixedrange': True,
                        'showticklabels':False,
                        'zeroline':False,
                    },
                    'yaxis': {
                        'range': [-10, 10],
                        'showgrid': False,
                        'fixedrange': True,
                        'showticklabels':False,
                        'zeroline':False
                    },
                    'shapes': [{
                        'type': 'rect',
                        'x0': c1,
                        'y0': -3,
                        'x1': c2,
                        'y1': 3,
                        'line': {'color':
                            'rgba(128, 0, 128, 0.7)'
                        }},
                        {
                        'type': 'rect',
                        'x0': c1-100,
                        'y0': -8,
                        'x1': c1,
                        'y1': 8,
                        'line': {'color':
                            'rgba(128, 0, 128, 1)'
                        }},
                        {
                        'type': 'rect',
                        'x0': c2,
                        'y0': -8,
                        'x1': c2+100,
                        'y1': 8,
                        'line': {'color':
                            'rgba(128, 0, 128, 1)'
                        }}]
                        }}),],
            className="three columns"))


    return plots



if __name__ == '__main__':

    app.run_server(debug=True)




    # go.Layout( 
    #             xaxis= {
    #                     'range': [c1-100, c2+100],
    #                     'showgrid': False,
    #                     'fixedrange': True,
    #                     'showticklabels':False,
    #                     'zeroline':False,
    #             },
    #             yaxis= {
    #                     'range': [-10, 10],
    #                     'showgrid': False,
    #                     'fixedrange': True,
    #                     'showticklabels':False,
    #                     'zeroline':False
    #             },
    #             shapes= [{
    #                     'type': 'rect',
    #                     'x0': c1,
    #                     'y0': -3,
    #                     'x1': c2,
    #                     'y1': 3,
    #                     'line': {'color':
    #                         'rgba(128, 0, 128, 0.7)'
    #                     }},
    #                     {
    #                     'type': 'rect',
    #                     'x0': c1-100,
    #                     'y0': -8,
    #                     'x1': c1,
    #                     'y1': 8,
    #                     'line': {'color':
    #                         'rgba(128, 0, 128, 1)'
    #                     }},
    #                     {
    #                     'type': 'rect',
    #                     'x0': c2,
    #                     'y0': -8,
    #                     'x1': c2+100,
    #                     'y1': 8,
    #                     'line': {'color':
    #                         'rgba(128, 0, 128, 1)'
    #                     }}]))


            #  f = go.Scatter(
            #     x=[-100,c2-c1],
            #     y=[-10,10],
                
            #     #text=[c1]
            #     mode='markers',
            #     marker={
            #             'symbol': "triangle-down",
            #             'size' : 15,
            #             'opacity':0
            #     })

            # middle.append({'type': 'rect',
            #                 'x0': 0,
            #                 'y0': -3,
            #                 'x1': c2-c1,
            #                 'y1': 3,
            #                 'line': {'color':
            #                     'rgba(128, 0, 128, 0.7)'
            #                 }})
            # left.append({'type': 'rect',
            #                 'x0': -100,
            #                 'y0': -3,
            #                 'x1': 0,
            #                 'y1': 3,
            #                 'line': {'color':
            #                     'rgba(128, 100, 128, 1)'
            #                 }})
            # right.append({'type': 'rect',
            #                 'x0': c2-c1,
            #                 'y0': -3,
            #                 'x1': (c2-c1)+100,
            #                 'y1': 3,
            #                 'line': {'color':
            #                     'rgba(128, 10, 128, 1)'
            #                 }})