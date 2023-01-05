import sys
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from lenspy import DynamicPlot
import dash
from dash import Dash, dcc, html, Input, Output, ctx
import dash_bio as dashbio

#https://dash.plotly.com/dash-bio/ideogram

app = dash.Dash(__name__)

app.layout = html.Div([
    'Select which chromosomes to display on the ideogram below:',
    dcc.Dropdown(
        id='my-default-displayed-chromosomes',
        options=[{'label': str(i), 'value': str(i)} for i in range(1, 23)],
        multi=True,
        value=[str(i) for i in range(1, 23)]
    ),
    dashbio.Ideogram(
        id='my-default-dashbio-ideogram'
    ),
    html.Div(id='my-default-ideogram-rotated')
])

@app.callback(
    Output('my-default-dashbio-ideogram', 'chromosomes'),
    Input('my-default-displayed-chromosomes', 'value')
)
def update_ideogram(value):
    return value

@app.callback(
    Output('my-default-ideogram-rotated', 'children'),
    Input('my-default-dashbio-ideogram', 'rotated')
)
def update_ideogram_rotated(rot):
    return 'You have {} selected a chromosome.'.format(
        '' if rot else 'not')

if __name__ == '__main__':
    app.run_server(host='127.0.0.1', debug=True)


