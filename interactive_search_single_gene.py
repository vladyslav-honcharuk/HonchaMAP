#!/usr/bin/env python3
"""
Interactive Radial Pattern Search - SINGLE GENE MODE
Select a gene → search by that gene's radial pattern only!
"""

import dash
from dash import dcc, html, Input, Output, State, callback
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import zarr
from pathlib import Path
from scipy.spatial.distance import cosine
import pandas as pd

# ====================== CONFIG ======================
RADIUS = 5
N_SHELLS = 10
TOP_K = 5
PATCH_STEP = 10  # How densely patches were generated
CUSTOM_COLORSCALE = [
    [0.00, "rgb(0, 0, 160)"],
    [0.25, "rgb(0, 160, 160)"],
    [0.50, "rgb(0, 200, 0)"],
    [0.75, "rgb(240, 200, 0)"],
    [1.00, "rgb(240, 0, 0)"],
]

# ====================== GLOBALS (will be filled on startup) ======================
database = []
genes = []
data_dir = None
sample_paths = {}


# ====================== CORE FUNCTIONS ======================
def load_sample(sample_path):
    with open(sample_path / "genes.csv") as f:
        genes_local = [line.strip() for line in f if line.strip()]
    zarr_file = sample_path / "zarr" / "bins_size_160.zarr.zip"
    store = zarr.storage.ZipStore(str(zarr_file), mode='r')
    expression = np.asarray(zarr.open_array(store, mode='r'))
    store.close()
    return genes_local, expression


def extract_patch(expression, center_x, center_y, radius):
    n_genes, width, height = expression.shape
    bins = []
    for x in range(max(0, center_x - radius), min(width, center_x + radius + 1)):
        for y in range(max(0, center_y - radius), min(height, center_y + radius + 1)):
            if (x - center_x)**2 + (y - center_y)**2 <= radius**2:
                bins.append((x, y))
    if not bins:
        return None
    x_coords = np.array([b[0] for b in bins])
    y_coords = np.array([b[1] for b in bins])
    patch_data = expression[:, x_coords, y_coords]
    return x_coords, y_coords, patch_data


def encode_radial_stats(x_coords, y_coords, patch_data, n_shells=N_SHELLS):
    """
    Encode ALL genes with radial statistics.
    Returns: (n_genes * n_shells * 4) array
    """
    n_genes = patch_data.shape[0]
    cx = x_coords.mean()
    cy = y_coords.mean()
    dists = np.sqrt((x_coords - cx)**2 + (y_coords - cy)**2)
    max_dist = dists.max() if dists.size > 0 else 0
    if max_dist == 0:
        stats = patch_data.mean(axis=1)
        return np.tile(stats, n_shells * 4)

    ring_edges = np.linspace(0, max_dist, n_shells + 1)
    features = []
    for gene_idx in range(n_genes):
        for ring_idx in range(n_shells):
            in_ring = (dists >= ring_edges[ring_idx]) & (dists < ring_edges[ring_idx + 1])
            if in_ring.sum() > 0:
                vals = patch_data[gene_idx, in_ring]
                features.extend([vals.mean(), vals.std(), vals.max(), vals.min()])
            else:
                features.extend([0.0, 0.0, 0.0, 0.0])
    return np.array(features, dtype=np.float32)


def extract_gene_features(encoding, gene_idx, n_shells=N_SHELLS):
    """
    Extract features for a single gene from full encoding.
    Returns: (n_shells * 4) array with [mean, std, max, min] per ring
    """
    start_idx = gene_idx * n_shells * 4
    end_idx = start_idx + n_shells * 4
    return encoding[start_idx:end_idx]


def build_database(data_dir_path):
    global database, genes, sample_paths
    data_dir_path = Path(data_dir_path)
    samples = [d for d in data_dir_path.iterdir() if d.is_dir() and (d / "zarr").exists()]
    all_genes = set()
    for s in samples:
        with open(s / "genes.csv") as f:
            all_genes.update(line.strip() for line in f if line.strip())
    genes = sorted(all_genes)

    database = []
    for sample_dir in samples:
        print(f"Loading {sample_dir.name}...")
        local_genes, expression = load_sample(sample_dir)
        global_expr = np.zeros((len(genes), expression.shape[1], expression.shape[2]), dtype=expression.dtype)
        for i, g in enumerate(local_genes):
            if g in genes:
                global_expr[genes.index(g)] = expression[i]

        w, h = expression.shape[1], expression.shape[2]
        for cx in range(RADIUS, w - RADIUS, PATCH_STEP):
            for cy in range(RADIUS, h - RADIUS, PATCH_STEP):
                result = extract_patch(global_expr, cx, cy, RADIUS)
                if result is not None:
                    x_coords, y_coords, patch_data = result
                    encoding = encode_radial_stats(x_coords, y_coords, patch_data)
                    database.append({
                        'sample': sample_dir.name,
                        'x': cx,
                        'y': cy,
                        'encoding': encoding  # Full encoding for all genes
                    })
        sample_paths[sample_dir.name] = sample_dir

    print(f"Database built: {len(database)} patches, {len(genes)} genes")
    return database, genes


def search_single_gene(query_encoding, gene_name, exclude_sample=None, top_k=TOP_K):
    """
    Search using ONLY the specified gene's radial pattern.

    Args:
        query_encoding: Full encoding with all genes
        gene_name: Gene to use for comparison
        exclude_sample: Sample name to exclude from results
        top_k: Number of top matches to return
    """
    if gene_name not in genes:
        print(f"Warning: {gene_name} not in gene list")
        return []

    gene_idx = genes.index(gene_name)

    # Extract features for just this gene from query
    query_gene_features = extract_gene_features(query_encoding, gene_idx, N_SHELLS)

    similarities = []
    valid_db = []

    for entry in database:
        if exclude_sample and entry['sample'] == exclude_sample:
            continue

        # Extract features for just this gene from database entry
        db_gene_features = extract_gene_features(entry['encoding'], gene_idx, N_SHELLS)

        # Compare ONLY this gene's features
        sim = 1 - cosine(query_gene_features, db_gene_features)
        similarities.append(sim)
        valid_db.append(entry)

    if not similarities:
        return []

    indices = np.argsort(similarities)[-top_k:][::-1]
    results = []
    for i in indices:
        entry = valid_db[i]
        results.append({
            'sample': entry['sample'],
            'x': entry['x'],
            'y': entry['y'],
            'similarity': similarities[i]
        })
    return results


def create_spatial_plot(sample_name, center_x=None, center_y=None, title="", gene_name=None, show_full_view=False):
    if gene_name is None or gene_name not in genes:
        gene_name = genes[0]

    sample_path = sample_paths[sample_name]
    local_genes, expression = load_sample(sample_path)
    gene_idx = genes.index(gene_name) if gene_name in genes else 0

    global_expr = np.zeros((len(genes), expression.shape[1], expression.shape[2]))
    for i, g in enumerate(local_genes):
        if g in genes:
            global_expr[genes.index(g)] = expression[i]

    w, h = expression.shape[1:]

    # Show full view if requested, otherwise show local region
    if show_full_view:
        x_min, x_max = 0, w
        y_min, y_max = 0, h
    else:
        margin = 30
        view_center_x = center_x if center_x is not None else w // 2
        view_center_y = center_y if center_y is not None else h // 2

        x_min = max(0, view_center_x - margin)
        x_max = min(w, view_center_x + margin)
        y_min = max(0, view_center_y - margin)
        y_max = min(h, view_center_y + margin)

    xs, ys, vals = [], [], []
    for x in range(int(x_min), int(x_max)):
        for y in range(int(y_min), int(y_max)):
            xs.append(x)
            ys.append(y)
            vals.append(float(global_expr[gene_idx, x, y]))

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode='markers',
        marker=dict(color=vals, colorscale=CUSTOM_COLORSCALE, size=8,
                    colorbar=dict(title=gene_name), showscale=True),
        text=[f"({x},{y})<br>{gene_name}: {v:.2f}" for x, y, v in zip(xs, ys, vals)],
        hoverinfo="text"
    ))

    if center_x is not None and center_y is not None:
        fig.add_trace(go.Scatter(x=[center_x], y=[center_y], mode="markers",
                                 marker=dict(color="red", size=14, symbol="x"),
                                 name="Center"))
        fig.add_shape(type="circle",
                      x0=center_x - RADIUS, y0=center_y - RADIUS,
                      x1=center_x + RADIUS, y1=center_y + RADIUS,
                      line=dict(color="red", width=2, dash="dash"))

    fig.update_layout(
        title=title,
        xaxis_title="X (bins)",
        yaxis_title="Y (bins)",
        width=500,
        height=500,
        showlegend=False,
        hovermode="closest",
        dragmode='lasso',
        yaxis=dict(
            scaleanchor="x",
            scaleratio=1,
            autorange="reversed"
        ),
        xaxis=dict(
            constrain="domain"
        )
    )

    return fig


# ====================== DASH APP ======================
app = dash.Dash(__name__, title="Single Gene Radial Search")

# Build database on startup
DATA_DIR = "/home/vlad/xenium_mundus/data_full/data/GPL33762/GSE243168"  # <<< CHANGE THIS PATH !!!
build_database(DATA_DIR)

# Pick a default sample and gene
default_sample = list(sample_paths.keys())[0] if sample_paths else ""
default_gene = genes[0] if genes else "Unknown"

app.layout = html.Div([
    html.H1("Interactive Radial Pattern Search - Single Gene Mode", style={'textAlign': 'center'}),
    html.P("Select a gene → lasso/box-select a point → see patches with similar radial patterns for THAT GENE ONLY!",
           style={'textAlign': 'center', 'fontSize': '16px', 'color': '#2c3e50'}),

    dcc.Store(id='selected-point', data={'x': None, 'y': None, 'sample': default_sample}),

    html.Div([
        html.Div([
            html.H3("Query Sample - Select Point"),
            html.Div([
                html.Label("Sample:"),
                dcc.Dropdown(id='query-sample-dropdown',
                           options=[{'label': s, 'value': s} for s in sample_paths.keys()],
                           value=default_sample,
                           style={'marginBottom': '10px'}),

                html.Label("Gene (searches by this gene's radial pattern ONLY):"),
                dcc.Dropdown(id='gene-dropdown',
                           options=[{'label': g, 'value': g} for g in genes],
                           value=default_gene,
                           style={'marginBottom': '10px'}),
            ]),
            dcc.Graph(id='query-plot', config={'displayModeBar': True, 'toImageButtonOptions': {'format': 'png'}}),
            html.Div(id='selected-info', style={'margin': '10px', 'fontWeight': 'bold', 'color': '#2980b9'})
        ], style={'width': '48%', 'display': 'inline-block', 'verticalAlign': 'top'}),

        html.Div([
            html.H3(f"Top {TOP_K} Matches (by selected gene's radial pattern)"),
            html.Div(id='matches-container', style={'display': 'flex', 'flexWrap': 'wrap', 'justifyContent': 'space-around'})
        ], style={'width': '48%', 'display': 'inline-block', 'verticalAlign': 'top', 'marginLeft': '4%'})
    ])
])


@callback(
    Output('query-plot', 'figure'),
    Input('query-sample-dropdown', 'value'),
    Input('gene-dropdown', 'value'),
    Input('selected-point', 'data')
)
def update_query_plot(sample_name, gene_name, selected_data):
    if not sample_name:
        return go.Figure().update_layout(title="No sample loaded")
    x = selected_data['x'] if selected_data else None
    y = selected_data['y'] if selected_data else None
    title = f"Query: {sample_name} | Gene: {gene_name}<br><sup>Click or lasso select a point</sup>"
    return create_spatial_plot(sample_name, center_x=x, center_y=y, title=title, gene_name=gene_name)


@callback(
    Output('selected-info', 'children'),
    Output('selected-point', 'data'),
    Output('matches-container', 'children'),
    Input('query-plot', 'selectedData'),
    Input('query-plot', 'clickData'),
    State('query-sample-dropdown', 'value'),
    State('selected-point', 'data'),
    State('gene-dropdown', 'value')
)
def on_selection(selectedData, clickData, current_sample, stored_point, gene_name):
    ctx = dash.callback_context
    triggered = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None

    new_x = new_y = None
    info = "No point selected yet."

    if triggered == 'query-plot':
        if selectedData and selectedData.get('points'):
            points = selectedData['points']
            if len(points) > 1:
                x_coords = [pt['x'] for pt in points]
                y_coords = [pt['y'] for pt in points]
                new_x = int(round(np.mean(x_coords)))
                new_y = int(round(np.mean(y_coords)))
                info = f"Selected region center: ({new_x}, {new_y}) | Searching by {gene_name} radial pattern"
            elif len(points) == 1:
                pt = points[0]
                new_x = int(round(pt['x']))
                new_y = int(round(pt['y']))
                info = f"Selected: ({new_x}, {new_y}) | Searching by {gene_name} radial pattern"
        elif clickData and clickData.get('points'):
            pt = clickData['points'][0]
            new_x = int(round(pt['x']))
            new_y = int(round(pt['y']))
            info = f"Clicked: ({new_x}, {new_y}) | Searching by {gene_name} radial pattern"

    if new_x is None and stored_point and stored_point.get('x') is not None:
        new_x = stored_point['x']
        new_y = stored_point['y']
        info = f"Current query: ({new_x}, {new_y}) | Gene: {gene_name}"
        if triggered != 'query-plot':
            return info, stored_point, dash.no_update

    new_stored = {'x': new_x, 'y': new_y, 'sample': current_sample}

    match_plots = [html.P("Select a point to see matches...", style={'color': 'gray'})]

    if new_x is not None and current_sample and triggered == 'query-plot':
        try:
            sample_path = sample_paths[current_sample]
            local_genes, expr = load_sample(sample_path)
            global_expr = np.zeros((len(genes), expr.shape[1], expr.shape[2]))
            for i, g in enumerate(local_genes):
                if g in genes:
                    global_expr[genes.index(g)] = expr[i]

            patch_result = extract_patch(global_expr, new_x, new_y, RADIUS)
            if patch_result is not None:
                x_coords, y_coords, patch_data = patch_result
                query_enc = encode_radial_stats(x_coords, y_coords, patch_data, n_shells=N_SHELLS)

                # SINGLE GENE SEARCH!
                results = search_single_gene(query_enc, gene_name, exclude_sample=current_sample, top_k=TOP_K)

                match_plots = []

                for i, res in enumerate(results):
                    fig = create_spatial_plot(
                        res['sample'],
                        res['x'],
                        res['y'],
                        title=f"#{i+1} | {res['sample']}<br>({res['x']},{res['y']})<br>Sim: {res['similarity']:.3f}",
                        gene_name=gene_name,
                        show_full_view=True
                    )
                    match_plots.append(
                        html.Div([
                            dcc.Graph(figure=fig, config={'displayModeBar': False},
                                    style={'width': '400px', 'height': '400px'})
                        ], style={'margin': '10px'})
                    )
            else:
                match_plots = [html.P("No valid patch at this location", style={'color': 'orange'})]

        except Exception as e:
            match_plots = [html.P(f"Error: {str(e)}", style={'color': 'red'})]
            print(f"Callback error: {e}")
            import traceback
            traceback.print_exc()
    elif stored_point and stored_point.get('x') is not None:
        return info, new_stored, dash.no_update

    return info, new_stored, match_plots


if __name__ == '__main__':
    try:
        print(f"Database: {len(database)} patches")
        print(f"Genes: {len(genes)} total")
        print(f"Samples: {list(sample_paths.keys())}")
        print("=" * 60)
        print("SINGLE GENE MODE ENABLED")
        print("Select a gene → search by ONLY that gene's radial pattern")
        print(f"Each gene has {N_SHELLS * 4} features (mean/std/max/min per ring)")
        print("=" * 60)
        print("Starting Dash app on http://127.0.0.1:8050")
        app.run(debug=True, port=8050)
    except Exception as e:
        print(f"Failed to start: {e}")
        import traceback
        traceback.print_exc()
