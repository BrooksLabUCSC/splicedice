from __future__ import print_function


########################################################################
# File: tbd.py
#  executable: tbd.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 03/10/2019 Created
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




fig2 = {
    'data': [{
        'x' : [1,2,3],
        'y' : [1,2,3],
        'type' : 'scatter',
       	'mode' : 'markers'
        }],
    'layout': {
                'clickmode': 'event+select'
            }



app.layout = html.Div([])



# @app.callback(
#     dash.dependencies.Output('plot', 'figure'),
#     [dash.dependencies.Input('group', 'value'),
#     dash.dependencies.Input('sample_group_type', 'value')])

# def update_graph(group, sample_group_type):