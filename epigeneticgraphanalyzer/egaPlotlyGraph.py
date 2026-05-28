import os
import networkx as nx
import argparse
import plotly.graph_objs as go


def parse_args():
    parser = argparse.ArgumentParser(description='Generate a Plotly graph from a NetworkX graph.')
    parser.add_argument('graph', type=str, help='Path to the input graph file (GraphML).')
    parser.add_argument('output_path', type=str, help='Directory where the output HTML file will be saved.')
    parser.add_argument('file_name_prefix', type=str, help='Prefix for the output HTML file name.')
    parser.add_argument('--partitionId', type=str, default='', help='Optional partition id suffix used in the output filename.')
    return parser.parse_args()


def writePlotlyGraph(partition_id, output_path, graph, file_name_prefix):
    pos = nx.spring_layout(graph)
    edge_x = []
    edge_y = []
    for edge in graph.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')
    node_x = []
    node_y = []
    node_text = []
    for node in graph.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(str(node))
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(graph.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: '+str(len(adjacencies[1])))
    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text
    fig = go.Figure(data=[edge_trace, node_trace],
                layout=go.Layout(
                title='Network graph made with Plotly',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                annotations=[ dict(
                    text="Python code: <a href='https://plotly.com/python/'> https://plotly.com/python/</a>",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002 ) ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )
    fig.write_html(os.path.join(output_path, file_name_prefix + str(partition_id) + ".html"))


def main():
    args = parse_args()
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    graph = nx.read_graphml(args.graph)
    writePlotlyGraph(args.partitionId, args.output_path, graph, args.file_name_prefix)