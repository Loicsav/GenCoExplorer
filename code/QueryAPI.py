from flask import Flask, flash, jsonify, render_template_string, request, make_response
import numpy as np
import pandas as pd
import os
import re
from io import BytesIO

app = Flask(__name__)

# Distintos ficheros
annotations_file = "./final_subgraphs_preprocessed_IC.csv" 
modules_dir = "./scCoExpNets-master/inst/data/networks/"

bulk_annotations_file = "./APP_ROSMAP_TGCNs_annotations.csv"
bulk_modules_file = "./APP_ROSMAP_TGCNs_modules.csv"

minimally_expressed_df = pd.read_csv("./Minimally_Expressed_Statistics.csv")
relevant_at_t0_df = pd.read_csv("./Relevant_At_T0_Statistics.csv")
relevant_in_all_iterations_df = pd.read_csv("./Relevant_In_All_Iterations_Statistics.csv")

def load_csv(file_path):
    # Cargar un archivo CSV desde la ruta especificada.
    if os.path.exists(file_path):
        return pd.read_csv(file_path)
    return None

def extract_iteration(file_name):
    match = re.search(r"_T(\d+)_", file_name)
    return f"T{match.group(1)}" if match else ""

def extract_cell_type(file_name):
    # "DA_like_neurons_11_T0_modules.csv" -> "DA_like_neurons"
    match = re.match(r"(.+?)_\d+_", file_name)
    return match.group(1).replace("_", " ") if match else ""

def format_cell_type(cell_type):
    # Formatea el nombre del Cell type eliminando el guion bajo y separando el número.
    if "_" in cell_type:
        parts = cell_type.split("_")
        return f"{parts[0]} {parts[1]}"  
    return cell_type



def format_p_value(p_value):
    # Pocos decimales de la notación científica
    if isinstance(p_value, (int, float)):
        if p_value <= 1e-3:
            return f"{p_value:.2E}"
        return f"{p_value:.4f}"
    return str(p_value)

def format_intersection(intersection_str):
    if not intersection_str:
        return ""
    genes = sorted([gene.strip() for gene in intersection_str.split(",") if gene.strip()]) 
    if len(genes) > 5:
        displayed_genes = genes[:5]
        return ", ".join(displayed_genes) + ", ..." 
    else:
        return ", ".join(genes) 

def format_ic(ic_value):
    # Formatea el valor de IC a 3 decimales y lo devuelve como float.
    try:
        return round(float(ic_value), 3)
    except (ValueError, TypeError):
        return 0.0

#Funciones para el autocompletado
def get_available_genes():
    genes = set()
    for file_name in os.listdir(modules_dir):
        if file_name.endswith(".csv"):
            dataset = pd.read_csv(os.path.join(modules_dir, file_name))
            if 'gene' in dataset.columns:
                genes.update(dataset['gene'].dropna().unique())
    return sorted(genes)

# Endpoint para obtener genes que coincidan con un patrón
@app.route('/api/genes')
def api_genes():
    search_term = request.args.get('term', '').upper()
    genes = get_available_genes()
    matching_genes = [gene for gene in genes if gene.startswith(search_term)][:100] 
    return jsonify(matching_genes)

def get_available_terms():
    annotations = load_csv(annotations_file)
    if annotations is not None and 'term_name' in annotations.columns:
        return sorted(annotations['term_name'].dropna().unique())
    return []

@app.route('/api/terms')
def api_terms():
    search_term = request.args.get('term', '').lower()
    terms = get_available_terms()
    matching_terms = [term for term in terms if term.lower().startswith(search_term)][:100]
    return jsonify(matching_terms)

#Funciones para filtros dinámicos
def get_available_clusters(cell_type):
    clusters = set()
    for file_name in os.listdir(modules_dir):
        if file_name.endswith(".csv"):
            current_cell_type = extract_cell_type(file_name)
            if current_cell_type == cell_type:
                dataset = pd.read_csv(os.path.join(modules_dir, file_name))
                if 'subcluster' in dataset.columns:
                    clusters.update(dataset['subcluster'].unique())
    return sorted(clusters)

def get_available_iterations(cell_type, cluster=None):
    iterations = set()
    for file_name in os.listdir(modules_dir):
        if file_name.endswith(".csv"):
            current_cell_type = extract_cell_type(file_name)
            if current_cell_type == cell_type:
                dataset = pd.read_csv(os.path.join(modules_dir, file_name))
                if cluster is None or ('subcluster' in dataset.columns and cluster in dataset['subcluster'].unique()):
                    iteration = extract_iteration(file_name)
                    if iteration:
                        iterations.add(iteration)
    return sorted(iterations)

@app.route('/api/clusters')
def api_clusters():
    cell_type = request.args.get('cell_type', '').strip()
    if not cell_type:
        return jsonify([])
    # Obtener clusters únicos para el cell_type especificado
    clusters = set()
    for file_name in os.listdir(modules_dir):
        if file_name.endswith(".csv"):
            current_cell_type = extract_cell_type(file_name)
            if current_cell_type == cell_type:
                try:
                    df = pd.read_csv(os.path.join(modules_dir, file_name))
                    if 'subcluster' in df.columns:
                        clusters.update(int(cluster) for cluster in df['subcluster'].unique())
                except Exception as e:
                    print(f"Error reading {file_name}: {e}")
    
    return jsonify(sorted(clusters))

@app.route('/api/iterations')
def api_iterations():
    cell_type = request.args.get('cell_type', '').strip()
    cluster = request.args.get('cluster', '').strip()
    if not cell_type:
        return jsonify([])
    iterations = set()
    for file_name in os.listdir(modules_dir):
        if file_name.endswith(".csv"):
            current_cell_type = extract_cell_type(file_name)
            if current_cell_type == cell_type:
                try:
                    df = pd.read_csv(os.path.join(modules_dir, file_name))
                    if cluster:
                        if 'subcluster' in df.columns:
                            if str(cluster) in map(str, df['subcluster'].unique()):
                                iteration = extract_iteration(file_name)
                                if iteration:
                                    iterations.add(iteration)
                    else:
                        iteration = extract_iteration(file_name)
                        if iteration:
                            iterations.add(iteration)
                except Exception as e:
                    print(f"Error reading {file_name}: {e}")
    
    return jsonify(sorted(iterations))

# Solo para modules
def query_dataset(file_type, search_term, cell_type_filter=None, iteration_filter=None, cluster_filter=None, module_filter=None, percentile_filter=None):
    data_dir = modules_dir
    search_column = ["gene"]
    results = []
    for file_name in os.listdir(data_dir):
        if file_name.endswith(".csv"):
            dataset = load_csv(os.path.join(data_dir, file_name))
            if dataset is None:
                continue
            for col in search_column:
                if col in dataset.columns:
                    filtered_data = dataset[dataset[col] == search_term]
                    if not filtered_data.empty:
                        for _, row in filtered_data.iterrows():
                            if percentile_filter is None or row["percentile"] >= percentile_filter:
                                results.append({
                                    "Iteration": extract_iteration(file_name),
                                    "Cell type": format_cell_type(extract_cell_type(file_name)),
                                    "Cluster": row["subcluster"],
                                    "Module": row["module"],
                                    "Module size": row["module_size"],
                                    "Gene": row["gene"],
                                    "Module membership": row["module_membership"],
                                    "Percentile (%)": row["percentile"],
                                })
    if cell_type_filter:
        filters = [f.strip().lower() for f in cell_type_filter.split(',')]
        results = [row for row in results if row["Cell type"].lower() in filters]
    if iteration_filter:
        filters = [f.strip() for f in iteration_filter.split(',')]
        results = [row for row in results if row["Iteration"] in filters]
    if cluster_filter:
        filters = [f.strip() for f in cluster_filter.split(',')]
        results = [row for row in results if str(row["Cluster"]) in filters]
    if module_filter:
        filters = [f.strip() for f in module_filter.split(',')]
        results = [row for row in results if str(row["Module"]) in filters]
    
    results.sort(key=lambda x: -x["Percentile (%)"])
    return results

# Solo para anotaciones
def query_annotations(search_term, cell_type_filter=None, iteration_filter=None, cluster_filter=None, module_filter=None):
    dataset = load_csv(annotations_file)
    if dataset is None:
        return []
    results = []
    search_columns = ["term_name", "term_id"]
    for col in search_columns:
        if col in dataset.columns:
            filtered_data = dataset[dataset[col] == search_term]
            if not filtered_data.empty:
                for _, row in filtered_data.iterrows():
                    cell_type = row["cell_type"].replace("_", " ")
                    p_value_formatted = format_p_value(row["p_value"])
                    results.append({
                        "Iteration": row["iteration"],
                        "Cell type": cell_type,
                        "Cluster": row["cluster"],
                        "Module": row["module"],
                        "Term id": row["term_id"],
                        "Term name": row["term_name"],
                        "P-value": p_value_formatted,
                        "Intersection": row.get("intersection", ""),
                        "Length of Intersection": row.get("length_intersection", ""),
                        "Source": row["source"],
                        "Subgraph ID": row["subgraph_id"],
                        "Subgraph size": row.get("subgraph_size", ""),
                        "IC": format_ic(row.get("IC", "")) 
                    })
    if cell_type_filter:
        filters = [f.strip().lower().replace("_", " ") for f in cell_type_filter.split(',')]
        results = [row for row in results if row["Cell type"].lower() in filters]
    if iteration_filter:
        filters = [f.strip() for f in iteration_filter.split(',')]
        results = [row for row in results if row["Iteration"] in filters]
    if cluster_filter:
        filters = [f.strip() for f in cluster_filter.split(',')]
        results = [row for row in results if str(row["Cluster"]) in filters]
    if module_filter:
        filters = [f.strip() for f in module_filter.split(',')]
        results = [row for row in results if str(row["Module"]) in filters]
    
    results.sort(key=lambda x: float(x["P-value"]) if x["P-value"] else float('inf'))
    return results

# Página principal de la API
@app.route('/', methods=['GET', 'POST'])
@app.route('/home')
def home():
    return render_template_string("""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GeneCoExplorer</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {
            font-family: Arial, sans-serif;
            background-color: #1a2a4f;
            padding: 40px;
        }
        .header {
            text-align: center;
            margin-bottom: 60px;
        }
        .header h1 {
            font-size: 3.5rem;
            color: #ffffff; 
            font-weight: bold;
        }
        .content {
            max-width: 1000px;
            margin: 0 auto;
            background-color: #a0c4ff; 
            padding: 40px;
            border-radius: 12px;
            box-shadow: 0 0 20px rgba(0, 0, 0, 0.15);
        }
        .content h2 {
            font-size: 2.2rem;
            color: #1a2a4f;
            margin-bottom: 30px;
            font-weight: bold;
        }
        .content p {
            font-size: 1.4rem;
            line-height: 1.8;
            color: #1a2a4f; 
            margin-bottom: 30px;
        }
        .button-container {
            display: grid;
            grid-template-columns: repeat(3, 1fr); 
            gap: 20px; 
        }
        .button-container a {
            display: flex;
            justify-content: center;
            align-items: center;
            text-decoration: none;
            color: white;
            font-weight: bold;
            font-size: 1.5rem;
            background-color: #1a2a4f; 
            padding: 20px;
            border-radius: 10px; 
            transition: background-color 0.3s ease; 
        }
        .button-container a:hover {
            background-color: #2a3a5f; 
        }
        .footer {
            text-align: center;
            margin-top: 60px;
            font-size: 1.2rem;
            color: #ffffff; 
        }
        .footer p {
            margin-bottom: 10px;
        }
        .stats {

            margin-bottom: 40px;
            color: #ffffff;
            font-size: 1.8rem;
        }
        .stats span {
            margin: 0 10px;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>GeneCoExplorer</h1>
    </div>
    <div class="content">
        <h2>A new resource to identify potential cell type markers and predict genes functions at cell type level in the context of Parkinson's Disease</h2>
        <p>Select the category of interest:</p>
        <div class="button-container">
            <a href="/gene_symbol">Gene symbol</a>
            <a href="/gene_ontology_terms">Gene ontology terms</a>
            <a href="/cell_type">Cell type</a>
        </div>
    </div>
    <div class="footer">
        <div class="stats">
            <span><strong>Networks:</strong> 112</span> |
            <span><strong>Cell types:</strong> 24</span> |
            <span><strong>Annotations:</strong> 2,423</span> |
            <span><strong>Predicted functions:</strong> 601,656</span>
        </div>
        <p>GeneCoExplorer was developed by Manuel Salas Díaz, BSc student in Computer Science</p>
        <p>GeneCoExplorer has been developed under the supervision of Professor José T. Pepe Palma and postdoctoral researcher Alicia Gómez Pascual from the University of Murcia</p>
        <p>We also thank the collaboration of Professor Juan A. Botía Blaya, also from the University of Murcia</p>
    </div>
</body>
</html>
""")

# Tres principales tipos de queries
@app.route('/gene_symbol')
def gene_symbol():
    return render_template_string("""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Symbol</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: grid;
            grid-template-rows: auto auto 1fr;
            height: 100vh;
        }
        .main-header {
            grid-column: 1 / -1;
            background-color: white;
            color: white;
            padding: 0;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            text-align: center;
            margin-bottom: 25px;
        }
        .main-header a {
            background-color: #1a2a3f;
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 20px;
            display: block;
            font-size: 1.4rem;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.8);
        }
        .main-header a:hover {
            background-color: #6c84b4;
        }
        .main-header a.active {
            background-color: #6c84b4;
        }
        .main-header a:nth-child(1) {
            margin-right: 30px;
            border-bottom-right-radius: 8px;
        }
        .main-header a:nth-child(2) {
            margin: 0 30px;
            border-bottom-left-radius: 8px;
            border-bottom-right-radius: 8px;
        }
        .main-header a:nth-child(3) {
            margin-left: 30px;
            border-bottom-left-radius: 8px;
        }
        .sub-header {
            grid-column: 1 / -1;
            background-color: white;  
            display: grid;
            grid-template-columns: 1fr minmax(150px, 500px) 100px minmax(150px, 500px) 100px minmax(150px, 500px) 1fr;
            text-align: center;
            padding: 0;
            align-items: stretch;
            margin-bottom: 10px;
        }
        .sub-header a {
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 15px 5px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.25rem;
            background-color: #354e7c;
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 1);
        }
        .sub-header a:nth-child(1) { grid-column: 2; }
        .sub-header a:nth-child(2) { grid-column: 4; }
        .sub-header a:nth-child(3) { grid-column: 6; }
        .sub-header a:hover {
            background-color: #6c84b4;
        }
        .sub-header a.active {
            background-color: #6c84b4;
        }
        .content {
            grid-row: 3;
            padding: 20px;
            overflow-y: auto;
        }
    </style>
</head>
<body>
    <div class="main-header">
        <a href="/gene_symbol" class="active">Gene symbol</a>
        <a href="/gene_ontology_terms">Gene ontology terms</a>
        <a href="/cell_type">Cell type</a>
    </div>
    <div class="sub-header">
        <a href="/gene_relevance">Gene relevance</a>
        <a href="/gene_functions">Gene functions</a>
        <a href="/new_gene_functions">New gene functions</a>
    </div>
    <div class="content">
        <h2>Gene Symbol Analysis</h2>
        <p>Select one of the subcategories above to explore gene-related data.</p>
    </div>
</body>
</html>
""")

@app.route('/gene_ontology_terms')
def gene_ontology_terms():
    return render_template_string("""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Ontology Terms</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: grid;
            grid-template-rows: auto auto 1fr;
            height: 100vh;
        }
        .main-header {
            grid-column: 1 / -1;
            background-color: white;
            color: white;
            padding: 0;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            text-align: center;
            margin-bottom: 25px;
        }
        .main-header a {
            background-color: #1a2a3f;
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 20px;
            display: block;
            font-size: 1.4rem;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.8);
        }
        .main-header a:hover {
            background-color: #6c84b4;
        }
        .main-header a.active {
            background-color: #6c84b4;
        }
        .main-header a:nth-child(1) {
            margin-right: 30px;
            border-bottom-right-radius: 8px;
        }
        .main-header a:nth-child(2) {
            margin: 0 30px;
            border-bottom-left-radius: 8px;
            border-bottom-right-radius: 8px;
        }
        .main-header a:nth-child(3) {
            margin-left: 30px;
            border-bottom-left-radius: 8px;
        }
        .sub-header {
            grid-column: 1 / -1;
            background-color: white;  
            display: grid;
            grid-template-columns: 1fr minmax(250px, 500px) 1fr;
            text-align: center;
            padding: 0;
            align-items: stretch;
            margin-bottom: 10px;
        }
        .sub-header a {
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 15px 5px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.25rem;
            background-color: #354e7c;
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 1);
        }
        .sub-header a:nth-child(1) { grid-column: 2; }
        .sub-header a:nth-child(2) { grid-column: 4; }
        .sub-header a:nth-child(3) { grid-column: 6; }
        .sub-header a:hover {
            background-color: #6c84b4;
        }
        .sub-header a.active {
            background-color: #6c84b4;
        }
        .content {
            grid-row: 3;
            padding: 20px;
            overflow-y: auto;
        }
    </style>
</head>
<body>
    <div class="main-header">
        <a href="/gene_symbol">Gene symbol</a>
        <a href="/gene_ontology_terms" class="active">Gene ontology terms</a>
        <a href="/cell_type">Cell type</a>
    </div>
    <div class="sub-header">
        <a href="/go_term_relevance">GO term relevance</a>
    </div>
    <div class="content">
        <h2>Gene Ontology Terms Analysis</h2>
        <p>Select the subcategory above to explore GO term-related data.</p>
    </div>
</body>
</html>
""")

@app.route('/cell_type')
def cell_type():
    return render_template_string("""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Cell Type</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: grid;
            grid-template-rows: auto auto 1fr;
            height: 100vh;
        }
        .main-header {
            grid-column: 1 / -1;
            background-color: white;
            color: white;
            padding: 0;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            text-align: center;
            margin-bottom: 25px;
        }
        .main-header a {
            background-color: #1a2a3f;
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 20px;
            display: block;
            font-size: 1.4rem;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.8);
        }
        .main-header a:hover {
            background-color: #6c84b4;
        }
        .main-header a.active {
            background-color: #6c84b4;
        }
        .main-header a:nth-child(1) {
            margin-right: 30px;
            border-bottom-right-radius: 8px;
        }
        .main-header a:nth-child(2) {
            margin: 0 30px;
            border-bottom-left-radius: 8px;
            border-bottom-right-radius: 8px;
        }
        .main-header a:nth-child(3) {
            margin-left: 30px;
            border-bottom-left-radius: 8px;
        }
        .sub-header {
            grid-column: 1 / -1;
            background-color: white;  
            display: grid;
            grid-template-columns: 1fr minmax(200px, 500px) 100px minmax(200px, 500px) 1fr;
            text-align: center;
            padding: 0;
            align-items: stretch;
            margin-bottom: 10px;
        }
        .sub-header a {
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 15px 5px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.25rem;
            background-color: #354e7c;
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 1);
        }
        .sub-header a:nth-child(1) { grid-column: 2; }
        .sub-header a:nth-child(2) { grid-column: 4; }
        .sub-header a:nth-child(3) { grid-column: 6; }
        .sub-header a:hover {
            background-color: #6c84b4;
        }
        .sub-header a.active {
            background-color: #6c84b4;
        }
        .content {
            grid-row: 3;
            padding: 20px;
            overflow-y: auto;
        }
    </style>
</head>
<body>
    <div class="main-header">
        <a href="/gene_symbol">Gene symbol</a>
        <a href="/gene_ontology_terms">Gene ontology terms</a>
        <a href="/cell_type" class="active">Cell type</a>
    </div>
    <div class="sub-header">
        <a href="/exclusive_relevant_genes">Exclusive relevant genes</a>
        <a href="/exclusive_go_terms">Exclusive GO terms</a>
    </div>
    <div class="content">
        <h2>Cell Type Analysis</h2>
        <p>Select one of the subcategories above to explore cell type-specific data.</p>
    </div>
</body>
</html>
""")

# Primer subquery de gene symbol
@app.route('/', methods=['GET', 'POST'])
@app.route('/gene_relevance', methods=['GET', 'POST'])
def gene_relevance():
    data_source = request.form.get('data_source', 'scRNA').strip()
    search_term = request.form.get('gene_name', '').strip()
    
    if data_source == 'scRNA':
        cell_type_filter = request.form.get('cell_type_filter', '').strip()
        iteration_filter = request.form.get('iteration_filter', '').strip()
        cluster_filter = request.form.get('cluster_filter', '').strip()
        module_filter = request.form.get('module_filter', '').strip()
        percentile_filter = request.form.get('percentile_filter', '').strip()
        percentile_filter = float(percentile_filter) if percentile_filter else None

        file_type = "modules"
        results = query_dataset(file_type, search_term, cell_type_filter, iteration_filter, 
                              cluster_filter, module_filter, percentile_filter)

        minimally_expressed_stats = minimally_expressed_df[minimally_expressed_df["Gene"] == search_term]
        relevant_at_t0_stats = relevant_at_t0_df[relevant_at_t0_df["Gene"] == search_term]
        relevant_in_all_iterations_stats = relevant_in_all_iterations_df[relevant_in_all_iterations_df["Gene"] == search_term]

        minimally_expressed = minimally_expressed_stats["Statistic"].values[0] if not minimally_expressed_stats.empty else "N/A"
        minimally_expressed_percentage = minimally_expressed_stats["Percentage"].values[0] if not minimally_expressed_stats.empty else "N/A"
        minimally_expressed_cell_types = minimally_expressed_stats["Cell Types"].values[0] if not minimally_expressed_stats.empty else "N/A"

        relevant_at_t0 = relevant_at_t0_stats["Statistic"].values[0] if not relevant_at_t0_stats.empty else "N/A"
        relevant_at_t0_percentage = relevant_at_t0_stats["Percentage"].values[0] if not relevant_at_t0_stats.empty else "N/A"
        relevant_at_t0_cell_types = relevant_at_t0_stats["Cell Types"].values[0] if not relevant_at_t0_stats.empty else "N/A"

        relevant_in_all_iterations = relevant_in_all_iterations_stats["Statistic"].values[0] if not relevant_in_all_iterations_stats.empty else "N/A"
        relevant_in_all_iterations_percentage = relevant_in_all_iterations_stats["Percentage"].values[0] if not relevant_in_all_iterations_stats.empty else "N/A"
        relevant_in_all_iterations_cell_types = relevant_in_all_iterations_stats["Cell Types"].values[0] if not relevant_in_all_iterations_stats.empty else "N/A"

        headers = ["Iteration", "Cell type", "Cluster", "Module", "Module size", "Gene", "Module membership", "Percentile (%)"]
        
        cell_types = set()
        for file_name in os.listdir(modules_dir):
            if file_name.endswith(".csv"):
                cell_type = extract_cell_type(file_name)
                if cell_type:
                    cell_types.add(cell_type)
        
        formatted_cell_types = sorted(cell_types, key=lambda x: (x.split()[0], int(x.split()[1]) if len(x.split()) > 1 and x.split()[1].isdigit() else 0))
        
        filter_dropdown = f"""
            <div class="mb-3">
                <label for="cell_type_filter" class="form-label">Cell type:</label>
                <select class="form-select" id="cell_type_filter" name="cell_type_filter" onchange="updateClusters()">
                    <option value="">All cell types</option>
                    {"".join([f'<option value="{cell_type}" {"selected" if cell_type == cell_type_filter else ""}>{cell_type}</option>' 
                    for cell_type in formatted_cell_types])}
                </select>
            </div>
            <div class="mb-3">
                <label for="cluster_filter" class="form-label">Cluster:</label>
                <select class="form-select" id="cluster_filter" name="cluster_filter" onchange="updateIterations()">
                    <option value="">All clusters</option>
                    {"".join([f'<option value="{cluster}" {"selected" if str(cluster) == cluster_filter else ""}>{cluster}</option>' 
                    for cluster in get_available_clusters(cell_type_filter) if cell_type_filter])}
                </select>
            </div>
            <div class="mb-3">
                <label for="iteration_filter" class="form-label">Iteration:</label>
                <select class="form-select" id="iteration_filter" name="iteration_filter">
                    <option value="">All iterations</option>
                    {"".join([f'<option value="{iteration}" {"selected" if iteration == iteration_filter else ""}>{iteration}</option>' 
                    for iteration in get_available_iterations(cell_type_filter, cluster_filter if cluster_filter else None) if cell_type_filter])}
                </select>
            </div>
            <div class="mb-3">
                <label for="module_filter" class="form-label">Module:</label>
                <input type="text" class="form-control" id="module_filter" name="module_filter" value="{module_filter}" placeholder="e.g., red, blue...">
            </div>
            <div class="mb-3">
                <label for="percentile_filter" class="form-label">Minimum Percentile:</label>
                <input type="number" class="form-control" id="percentile_filter" name="percentile_filter" value="{percentile_filter if percentile_filter is not None else ''}" placeholder="e.g., 90" min="0" max="100">
            </div>
        """
        
        # Función auxiliar dente de gene relevance para conseguir quitar los _ de los tipos celulares
        def format_cell_types(cell_types_str):
            if pd.isna(cell_types_str) or cell_types_str == "N/A":
                return "No data available"
            
            cell_types_list = [
                cell_type.replace("_", " ").strip() 
                for cell_type in cell_types_str.split("; ")
            ]
            cell_types_list.sort()
            return "<br>".join(cell_types_list)
            
        metric_text = f"""
            <div class="metric-container">
                <h3>{search_term} Module relevance statistics for all cell types</h3>
                <details>
                    <summary>Minimally expressed in {minimally_expressed} ({minimally_expressed_percentage})</summary>
                    <div class="cell-types-content">{format_cell_types(minimally_expressed_cell_types)}</div>
                </details>
                <details>
                    <summary>Relevant gene at T0 in {relevant_at_t0} ({relevant_at_t0_percentage})</summary>
                    <div class="cell-types-content">{format_cell_types(relevant_at_t0_cell_types)}</div>
                </details>
                <details>
                    <summary>Relevant gene in all iterations for {relevant_in_all_iterations} ({relevant_in_all_iterations_percentage})</summary>
                    <div class="cell-types-content">{format_cell_types(relevant_in_all_iterations_cell_types)}</div>
                </details>
            </div>
        """ if search_term else ""
        
    else:
        target_filter = request.form.get('target_filter', '').strip()
        tissue_filter = request.form.get('tissue_filter', '').strip()
        cutoff_filter = request.form.get('cutoff_filter', '').strip()
        module_filter = request.form.get('module_filter', '').strip()
        min_correlation = request.form.get('min_correlation', '').strip()
        min_correlation = float(min_correlation) if min_correlation else None

        bulk_modules = pd.read_csv(bulk_modules_file)
        
        results = []
        if search_term:
            filtered_data = bulk_modules[bulk_modules['gene'] == search_term]
            
            if target_filter:
                filtered_data = filtered_data[filtered_data['target'] == target_filter]
            if tissue_filter:
                filtered_data = filtered_data[filtered_data['tissue'] == tissue_filter]
            if cutoff_filter:
                filtered_data = filtered_data[filtered_data['cutoff'].astype(str) == cutoff_filter]
            if module_filter:
                filtered_data = filtered_data[filtered_data['module'] == module_filter]
            if min_correlation is not None:
                filtered_data = filtered_data[filtered_data['correlation'] >= min_correlation]
            
            filtered_data = filtered_data.sort_values(by='correlation', ascending=False)
            results = filtered_data.to_dict('records')

        headers = ["Cutoff", "Target", "Tissue", "Phenotype", "Module", "Module size", "Gene", "Correlation"]
        
        available_targets = sorted(bulk_modules['target'].unique())
        available_tissues = sorted(bulk_modules['tissue'].unique())
        available_cutoffs = sorted(bulk_modules['cutoff'].unique())
        
        filter_dropdown = f"""
            <div class="mb-3">
                <label for="target_filter" class="form-label">Target:</label>
                <select class="form-select" id="target_filter" name="target_filter">
                    <option value="">All targets</option>
                    {"".join([f'<option value="{target}" {"selected" if target == target_filter else ""}>{target}</option>' 
                    for target in available_targets])}
                </select>
            </div>
            <div class="mb-3">
                <label for="tissue_filter" class="form-label">Tissue:</label>
                <select class="form-select" id="tissue_filter" name="tissue_filter">
                    <option value="">All tissues</option>
                    {"".join([f'<option value="{tissue}" {"selected" if tissue == tissue_filter else ""}>{tissue}</option>' 
                    for tissue in available_tissues])}
                </select>
            </div>
            <div class="mb-3">
                <label for="cutoff_filter" class="form-label">Cutoff:</label>
                <select class="form-select" id="cutoff_filter" name="cutoff_filter">
                    <option value="">All cutoffs</option>
                    {"".join([f'<option value="{cutoff}" {"selected" if str(cutoff) == cutoff_filter else ""}>{cutoff}</option>' 
                    for cutoff in available_cutoffs])}
                </select>
            </div>
            <div class="mb-3">
                <label for="module_filter" class="form-label">Module:</label>
                <input type="text" class="form-control" id="module_filter" name="module_filter" 
                    value="{module_filter}" placeholder="e.g., PPP2CA">
            </div>
            <div class="mb-3">
                <label for="min_correlation" class="form-label">Minimum Correlation:</label>
                <input type="number" class="form-control" id="min_correlation" name="min_correlation" 
                    value="{min_correlation if min_correlation is not None else ''}" 
                    placeholder="e.g., 0.5" min="0" max="1" step="0.01">
            </div>
        """
        
        metric_text = ""

    table_rows = ""
    for result in results:
        if data_source == 'scRNA':
            result["Cell type"] = result["Cell type"].split("_")[0]
            
            row = "".join([
                f"<td class='text-center'>{result.get(col, '')}</td>" 
                for col in headers
            ])
        else:
            row = "".join([
                f"<td class='text-center'>{result.get(col.lower().replace(' ', '_'), '')}</td>"
                for col in headers
            ])
        
        table_rows += f"<tr>{row}</tr>"

    column_descriptions = {
        "scRNA": {
            "Iteration": "Number of times the pseudo-cell algorithm has been executed (e.g., T0, T5, etc.)",
            "Cell type": "Cell type where the gene is expressed",
            "Cluster": "Subgroup within the cell type",
            "Module": "Name of the module",
            "Module size": "Number of genes that make up the module",
            "Gene": "Gene symbol",
            "Module membership": "Connectivity of a gene in a module and its relevance (higher values indicate stronger association)",
            "Percentile (%)": "Relevance of a gene in a module represented from 1 to 100",
        },
        "bulk": {
            "Cutoff": "Indicates the stringency level, ranging from 1 to 10 where higher values correspond to more reliable modules",
            "Target": "Main variable of interest, typically a gene or a covariate associated with each individual",
            "Tissue": "Tissue sample origin",
            "Phenotype": "Disease or condition under study",
            "Module": "Name of the co-expression module",
            "Module size": "Number of genes in the module",
            "Gene": "Gene symbol",
            "Correlation": "A measure of gene importance within the network",
        }
    }

    current_descriptions = column_descriptions['scRNA'] if data_source == 'scRNA' else column_descriptions['bulk']
    headers_with_tooltips = []
    for header in headers:
        description = current_descriptions.get(header, "No description available")
        headers_with_tooltips.append(
            f'''
            <th class="text-center">
                {header}
                <div class="tooltip-icon">?</div>
                <div class="tooltip-text">{description}</div>
            </th>
            '''
        )

    return render_template_string(f"""   
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Search</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: grid;
            grid-template-columns: 425px 1fr 475px;
            grid-template-rows: auto auto 1fr;
            height: 100vh;
        }}
        .main-header {{
            grid-column: 1 / -1;
            background-color: white;
            color: white;
            padding: 0;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            text-align: center;
            margin-bottom: 25px;
        }}
        .main-header a {{
            background-color: #1a2a3f;
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 20px;
            display: block;
            font-size: 1.4rem;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.8);
        }}
        .main-header a:hover {{
            background-color: #6c84b4;
        }}
        .main-header a.active {{
            background-color: #6c84b4;
        }}
        .main-header a:nth-child(1) {{
            margin-right: 30px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(2) {{
            margin: 0 30px;
            border-bottom-left-radius: 8px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(3) {{
            margin-left: 30px;
            border-bottom-left-radius: 8px;
        }}
        .sub-header {{
            grid-column: 1 / -1;
            background-color: white;  
            display: grid;
            grid-template-columns: 1fr minmax(150px, 500px) 100px minmax(150px, 500px) 100px minmax(150px, 500px) 1fr;
            text-align: center;
            padding: 0;
            align-items: stretch;
            margin-bottom: 10px;
        }}
        .sub-header a {{
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 15px 5px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.25rem;
            background-color: #354e7c;
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 1);
        }}
        .sub-header a:nth-child(1) {{ grid-column: 2; }}
        .sub-header a:nth-child(2) {{ grid-column: 4; }}
        .sub-header a:nth-child(3) {{ grid-column: 6; }}
        .sub-header a:hover {{
            background-color: #6c84b4;
        }}
        .sub-header a.active {{
            background-color: #6c84b4;
        }}
        .left-panel {{
            grid-column: 1;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .center-panel {{
            grid-column: 2;
            padding: 0px;
            overflow-y: auto;
        }}
        .right-panel {{
            grid-column: 3;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .table {{
            width: 100%;
            border-collapse: collapse;
        }}
        .table th, .table td {{
            padding: 10px;
            border: 1px solid #ddd;
            text-align: center;
        }}
        .table th {{
            background-color: #a0c4ff;
            color: black;
            font-weight: bold;
        }}
        .table td {{
            white-space: normal !important;
            word-wrap: break-word;
        }}
        .results-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0;
            padding: 10px;
            background-color: white;
            position: sticky;
            top: 0;
            z-index: 1000;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        .table-container {{
            margin-top: 0;
        }}
        .metric-container {{
            margin-top: 20px;
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric-container h3 {{
            color: #1a2a4f;
            border-bottom: 2px solid #a0c4ff;
            padding-bottom: 10px;
        }}
        .metric-container p {{
            margin-bottom: 15px;
            font-size: 1.1rem;
            line-height: 1.6;
        }}
        .cell-types-content {{
            margin-left: 20px;  
            padding-left: 10px;
            border-left: 2px solid #a0c4ff;  
            margin-top: 5px;
            margin-bottom: 10px;
        }}
        .filter-dropdown {{
            margin-top: 10px;
        }}
        .tooltip-icon {{
            display: inline-block;
            width: 18px;
            height: 18px;
            background-color: #1a2a4f;
            color: white;
            border-radius: 50%;
            text-align: center;
            font-size: 12px;
            line-height: 18px;
            cursor: help;
            margin-left: 5px;
        }}
        .tooltip-text {{
            visibility: hidden;
            width: 200px;
            background-color: white;
            color: black;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 12px;
            position: absolute;
            z-index: 1000;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
            opacity: 0;
            transition: opacity 0.3s;
            font-size: 14px;
            line-height: 1.5;
            text-align: left;
        }}
        th:hover .tooltip-text {{
            visibility: visible;
            opacity: 1;
        }}
        .ui-autocomplete {{
            max-height: 230px;
            overflow-y: auto;
            overflow-x: hidden;
            background-color: white;
            border: 1px solid #ddd;
            padding: 5px;
        }}
        .ui-menu-item {{
            padding: 5px;
        }}
        .ui-menu-item:hover {{
            background-color: #f0f0f0;
            cursor: pointer;
        }}
        .run-example-btn {{
            background-color: #e3a42b;  
            color: #ffffff;
            margin-left: auto;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .run-example-btn:hover {{
            background-color: #f9aa33;
            color: white;
        }}
    </style>
</head>
<body>
    <div class="main-header">
        <a href="/gene_symbol" class="active">Gene symbol</a>
        <a href="/gene_ontology_terms">Gene ontology terms</a>
        <a href="/cell_type">Cell type</a>
    </div>
    <div class="sub-header">
        <a href="/gene_relevance" class="active">Gene relevance</a>
        <a href="/gene_functions">Gene functions</a>
        <a href="/new_gene_functions">New gene functions</a>
    </div>
    <div class="left-panel">
        <form method="POST">
            <div class="mb-3">
                <label for="gene_name" class="form-label">Search for gene symbol:</label>
                <input type="text" class="form-control" id="gene_name" name="gene_name" 
                    value="{search_term}" placeholder="e.g., SNCA, AATF" 
                    autocomplete="off" required>
            </div>
            <div class="mb-3">
                <label for="data_source" class="form-label">Database:</label>
                <select class="form-select" id="data_source" name="data_source">
                    <option value="scRNA" {"selected" if data_source == 'scRNA' else ""}>scRNA-seq PD scCoExpNets</option>
                    <option value="bulk" {"selected" if data_source == 'bulk' else ""}>bulk RNA-seq AD TGCNs</option>
                </select>
            </div>
            <div class="mb-3">
                <button type="button" class="btn btn-success" onclick="toggleFilters()">Filters</button>
                <div id="filter-dropdown" class="filter-dropdown">
                    {filter_dropdown}
                </div>
            </div>
            <div class="d-flex justify-content-between align-items-center">
                <button type="submit" class="btn btn-primary">Search</button>
                <button type="button" class="btn run-example-btn" onclick="runExample()">Run example</button>
            </div>
        </form>
    </div>
    <div class="center-panel">
        <div class="results-header">
            <h2>{search_term} {"Module relevance per cell type" if data_source == 'scRNA' else "Modules relevance per cutoff"} </h2>
            <form method="POST" action="/download">
                <input type="hidden" name="data_source" value="{data_source}">
                <input type="hidden" name="file_type" value="modules">
                <input type="hidden" name="gene_name" value="{search_term}">
                {"".join([f'<input type="hidden" name="{name}" value="{value}">' 
                for name, value in request.form.items() if name not in ['data_source', 'file_type', 'gene_name']])}
                <div class="mb-3">
                    <label for="download_format" class="form-label">Choose format:</label>
                    <select class="form-select d-inline-block w-auto" id="download_format" name="download_format">
                        <option value="csv">CSV</option>
                        <option value="xlsx">XLSX</option>
                        <option value="html">HTML</option>
                    </select>
                    <button type="submit" class="btn btn-success">Download</button>
                </div>
            </form>
        </div>
        <div class="table-container">
            <table class="table table-bordered table-hover">
                <thead class="table-primary">
                    <tr>
                        {''.join(headers_with_tooltips)}
                    </tr>
                </thead>
                <tbody>
                    {table_rows if search_term else f"<tr><td colspan='{len(headers)}' class='text-center'>Enter a gene symbol to search</td></tr>"}
                </tbody>
            </table>
        </div>
    </div>
    <div class="right-panel">
        {metric_text}
    </div>
    <script>
        function toggleFilters() {{
            const filterDropdown = document.getElementById('filter-dropdown');
            filterDropdown.style.display = filterDropdown.style.display === 'none' ? 'block' : 'none';
        }}
        
        document.getElementById('data_source').addEventListener('change', function() {{
            this.form.submit();
        }});
    </script>
    <script src="/static/sort_table.js"></script>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://code.jquery.com/ui/1.13.1/jquery-ui.min.js"></script>
    <script>
    $(function() {{
        $("#gene_name").autocomplete({{
            source: function(request, response) {{
                $.getJSON("/api/genes", {{
                    term: request.term
                }}, response);
            }},
            minLength: 2,
            select: function(event, ui) {{
                $("#gene_name").val(ui.item.value);
                $(this).closest("form").submit();
            }},
            focus: function(event, ui) {{
                $("#gene_name").val(ui.item.value);
                return false;
            }}
        }});
    }});
    </script>
    <script src="/static/update_filters.js"></script>
    <script>
    function runExample() {{
        const form = document.querySelector('.left-panel form');
        form.gene_name.value = "SNCA";
        form.data_source.value = "scRNA";
        document.getElementById('cell_type_filter').value = "Oligodendrocytes";
        document.getElementById('module_filter').value = "turquoise";
        document.getElementById('percentile_filter').value = "80";
        form.submit();
    }}
    </script>
</body>
</html>
""")

# Query para predecir nuevas funciones
@app.route('/gene_functions', methods=['GET', 'POST'])
def gene_functions():
    data_source = request.form.get('data_source', 'scRNA').strip()
    search_term = request.form.get('gene_name', '').strip().upper()
    
    if data_source == 'scRNA':
        annotations_file_path = "./final_subgraphs_preprocessed_IC.csv"
        modules_dir_path = "./scCoExpNets-master/inst/data/networks/"
        minimally_expressed_file = "./Minimally_Expressed_Statistics.csv"
    else:
        annotations_file_path = "./APP_ROSMAP_TGCNs_annotations.csv"
        modules_dir_path = None  
    
    stats = {
        "minimally_expressed": "N/A",
        "minimally_expressed_percentage": "N/A",
        "mean_mm": "N/A",
        "mean_percentile": "N/A",
        "participation_percentage": "N/A",
        "mean_annotations": "N/A",
        "mean_ic": "N/A"
    }

    if data_source == 'scRNA':
        cell_type_filter = request.form.get('cell_type_filter', '').strip()
        iteration_filter = request.form.get('iteration_filter', '').strip()
        cluster_filter = request.form.get('cluster_filter', '').strip()
        module_filter = request.form.get('module_filter', '').strip()
    else:
        target_filter = request.form.get('target_filter', '').strip()
        tissue_filter = request.form.get('tissue_filter', '').strip()
        cutoff_filter = request.form.get('cutoff_filter', '').strip()
        module_filter = request.form.get('module_filter', '').strip()

    results = []
    if search_term:
        annotations_data = load_csv(annotations_file_path)
        
        if annotations_data is not None:
            mask = annotations_data['intersection'].str.contains(search_term, case=False, na=False)
            filtered_data = annotations_data[mask].copy()
            
            if data_source == 'scRNA':
                if cell_type_filter:
                    cell_filters = [f.strip().lower().replace("_", " ") for f in cell_type_filter.split(',')]
                    filtered_data = filtered_data[
                        filtered_data['cell_type'].str.replace("_", " ").str.lower().isin(cell_filters)
                    ]
                
                if iteration_filter:
                    iter_filters = [f.strip() for f in iteration_filter.split(',')]
                    filtered_data = filtered_data[filtered_data['iteration'].isin(iter_filters)]
                
                if cluster_filter:
                    cluster_filters = [f.strip() for f in cluster_filter.split(',')]
                    filtered_data = filtered_data[filtered_data['cluster'].astype(str).isin(cluster_filters)]
                
                if module_filter:
                    module_filters = [f.strip() for f in module_filter.split(',')]
                    filtered_data = filtered_data[filtered_data['module'].astype(str).isin(module_filters)]
            else:
                if target_filter:
                    filtered_data = filtered_data[filtered_data['target'] == target_filter]
                if tissue_filter:
                    filtered_data = filtered_data[filtered_data['tissue'] == tissue_filter]
                if cutoff_filter:
                    filtered_data = filtered_data[filtered_data['cutoff'].astype(str) == cutoff_filter]
                if module_filter:
                    filtered_data = filtered_data[filtered_data['module'] == module_filter]
            
            filtered_data = filtered_data.sort_values(by='p_value')
            
            results = filtered_data.to_dict('records')

            if data_source == 'scRNA' and len(results) > 0:
                if os.path.exists(minimally_expressed_file):
                    minimally_expressed_stats = pd.read_csv(minimally_expressed_file)
                    gene_stats = minimally_expressed_stats[minimally_expressed_stats['Gene'] == search_term]
                    if not gene_stats.empty:
                        stats["minimally_expressed"] = gene_stats.iloc[0]['Statistic']
                        stats["minimally_expressed_percentage"] = gene_stats.iloc[0]['Percentage']

                mm_values = []
                percentile_values = []
                for file_name in os.listdir(modules_dir_path):
                    if file_name.endswith(".csv"):
                        df_module = pd.read_csv(os.path.join(modules_dir_path, file_name))
                        gene_data = df_module[df_module['gene'] == search_term]
                        if not gene_data.empty:
                            mm_values.append(gene_data.iloc[0]['module_membership'])
                            percentile_values.append(gene_data.iloc[0]['percentile'])
                
                if mm_values:
                    stats["mean_mm"] = f"{sum(mm_values)/len(mm_values):.3f}"
                if percentile_values:
                    stats["mean_percentile"] = f"{sum(percentile_values)/len(percentile_values):.1f}"
                
                total_annotations = len(results)
                unique_cell_types = filtered_data.groupby(['cell_type', 'cluster']).size().reset_index().shape[0]
                stats["participation_percentage"] = f"{(total_annotations/unique_cell_types):.1f}" if unique_cell_types > 0 else "0"
                
                ic_values = []
                for r in results:
                    try:
                        ic = float(r.get('IC', 0))
                        if ic == float('inf'):
                            max_ic = annotations_data['IC'].replace([np.inf, -np.inf], np.nan).max()
                            ic = max_ic + 10 if not pd.isna(max_ic) else 0
                        ic_values.append(ic)
                    except (ValueError, TypeError):
                        pass
                
                if ic_values:
                    stats["mean_ic"] = f"{sum(ic_values)/len(ic_values):.2f}"
                    stats["mean_annotations"] = f"{(total_annotations/unique_cell_types):.1f}" if unique_cell_types > 0 else "0"

    if data_source == 'scRNA':
        headers = ["Iteration", "Cell type", "Cluster", "Module", "Term id", "Term name", 
                   "P-value", "Intersection", "Length of Intersection", "Source", 
                   "Subgraph ID", "Subgraph size", "IC"]
    else:
        headers = ["Cutoff", "Target", "Tissue", "Phenotype", "Module", "Term id", "Term name", 
                   "P-value", "Intersection", "Length of Intersection", "Source", "IC"]
    
    table_rows = ""
    for result in results:
        p_value_formatted = format_p_value(result.get("p_value", ""))
        intersection_formatted = format_intersection(result.get("intersection", ""))
        ic_formatted = format_ic(result.get("IC", ""))
        
        if data_source == 'scRNA':
            row = f"""
            <tr>
                <td class='text-center'>{result.get('iteration', '')}</td>
                <td class='text-center'>{result.get('cell_type', '').replace('_', ' ')}</td>
                <td class='text-center'>{result.get('cluster', '')}</td>
                <td class='text-center'>{result.get('module', '')}</td>
                <td class='text-center'>{result.get('term_id', '')}</td>
                <td class='text-center'>{result.get('term_name', '')}</td>
                <td class='text-center'>{p_value_formatted}</td>
                <td class='text-center'>{intersection_formatted}</td>
                <td class='text-center'>{result.get('length_intersection', '')}</td>
                <td class='text-center'>{result.get('source', '')}</td>
                <td class='text-center'>{result.get('subgraph_id', '')}</td>
                <td class='text-center'>{result.get('subgraph_size', '')}</td>
                <td class='text-center'>{ic_formatted}</td>
            </tr>
            """
        else:
            row = f"""
            <tr>
                <td class='text-center'>{result.get('cutoff', '')}</td>
                <td class='text-center'>{result.get('target', '')}</td>
                <td class='text-center'>{result.get('tissue', '')}</td>
                <td class='text-center'>{result.get('phenotype', '')}</td>
                <td class='text-center'>{result.get('module', '')}</td>
                <td class='text-center'>{result.get('term_id', '')}</td>
                <td class='text-center'>{result.get('term_name', '')}</td>
                <td class='text-center'>{p_value_formatted}</td>
                <td class='text-center'>{intersection_formatted}</td>
                <td class='text-center'>{result.get('length_intersection', '')}</td>
                <td class='text-center'>{result.get('source', '')}</td>
                <td class='text-center'>{ic_formatted}</td>
            </tr>
            """
        table_rows += row

    cell_types = set()
    if data_source == 'scRNA' and os.path.exists(modules_dir_path):
        for file_name in os.listdir(modules_dir_path):
            if file_name.endswith(".csv"):
                cell_type = extract_cell_type(file_name)
                if cell_type:
                    cell_types.add(cell_type)

    if data_source == 'scRNA':
        filter_dropdown = f"""
            <div class="mb-3">
                <label for="cell_type_filter" class="form-label">Cell type:</label>
                <select class="form-select" id="cell_type_filter" name="cell_type_filter" onchange="updateClusters()">
                    <option value="">All cell types</option>
                    {"".join([f'<option value="{cell_type}" {"selected" if cell_type == cell_type_filter else ""}>{cell_type}</option>' 
                    for cell_type in sorted(cell_types) if cell_type])}
                </select>
            </div>
            <div class="mb-3">
                <label for="cluster_filter" class="form-label">Cluster:</label>
                <select class="form-select" id="cluster_filter" name="cluster_filter" onchange="updateIterations()">
                    <option value="">All clusters</option>
                    {"".join([f'<option value="{cluster}" {"selected" if str(cluster) == cluster_filter else ""}>{cluster}</option>' 
                    for cluster in get_available_clusters(cell_type_filter) if cell_type_filter])}
                </select>
            </div>
            <div class="mb-3">
                <label for="iteration_filter" class="form-label">Iteration:</label>
                <select class="form-select" id="iteration_filter" name="iteration_filter">
                    <option value="">All iterations</option>
                    {"".join([f'<option value="{iteration}" {"selected" if iteration == iteration_filter else ""}>{iteration}</option>' 
                    for iteration in get_available_iterations(cell_type_filter, cluster_filter if cluster_filter else None) if cell_type_filter])}
                </select>
            </div>
            <div class="mb-3">
                <label for="module_filter" class="form-label">Module:</label>
                <input type="text" class="form-control" id="module_filter" name="module_filter" value="{module_filter}" placeholder="e.g., red, blue...">
            </div>
        """
    else:
        filter_dropdown = """
            <div class="mb-3">
                <label for="target_filter" class="form-label">Target:</label>
                <input type="text" class="form-control" id="target_filter" name="target_filter" placeholder="e.g., APP">
            </div>
            <div class="mb-3">
                <label for="tissue_filter" class="form-label">Tissue:</label>
                <input type="text" class="form-control" id="tissue_filter" name="tissue_filter" placeholder="e.g., Temporal Cortex">
            </div>
            <div class="mb-3">
                <label for="cutoff_filter" class="form-label">Cutoff:</label>
                <input type="text" class="form-control" id="cutoff_filter" name="cutoff_filter" placeholder="e.g., 5">
            </div>
            <div class="mb-3">
                <label for="module_filter" class="form-label">Module:</label>
                <input type="text" class="form-control" id="module_filter" name="module_filter" placeholder="e.g., PPP2CA">
            </div>
        """

    column_descriptions = {
        "scRNA": {
            "Iteration": "Number of times the pseudo-cell algorithm has been executed (e.g., T0, T5, etc.)",
            "Cell type": "Cell type where the gene is expressed",
            "Cluster": "Subgroup within the cell type",
            "Module": "Name of the module",
            "Term id": "Identifier of an annotation in the GO database",
            "Term name": "Full name of a GO annotation",
            "P-value": "Significance of an annotation given a set of genes in a co-expression module (lower = more significant)",
            "Intersection": "Name of the genes in a module that intersect with a GO database annotation",
            "Length of Intersection": "Number of genes in a module that intersect with a GO database annotation",
            "Source": "Database used to obtain the annotations (biological_process, cellular_component and molecular_function)",
            "Subgraph ID": "GO annotation subgraph identifier",
            "Subgraph size": "Number of GO annotations that make up a subgraph",
            "IC": "Index that tells us how informative an annotation is. The higher the CI, the more specific and informative the annotation is."
        },
        "bulk": {
            "Cutoff": "Indicates the stringency level, ranging from 1 to 10 where higher values correspond to more reliable modules",
            "Target": "Main variable of interest, typically a gene or a covariate associated with each individual",
            "Tissue": "Tissue sample origin",
            "Phenotype": "Disease or condition under study",
            "Module": "Name of the co-expression module",
            "Term id": "Identifier of an annotation in the GO database",
            "Term name": "Full name of a GO annotation",
            "P-value": "Significance of an annotation given a set of genes in a co-expression module (lower = more significant)",
            "Intersection": "Name of the genes in a module that intersect with a GO database annotation",
            "Length of Intersection": "Number of genes in a module that intersect with a GO database annotation",
            "Source": "Database used to obtain the annotations (biological_process, cellular_component and molecular_function)",
            "IC": "Index that tells us how informative an annotation is. The higher the CI, the more specific and informative the annotation is."
        }
    }

    return render_template_string(f"""   
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Function Search</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: grid;
            grid-template-columns: 425px 1fr 475px;
            grid-template-rows: auto auto 1fr;
            height: 100vh;
        }}
        .main-header {{
            grid-column: 1 / -1;
            background-color: white;
            color: white;
            padding: 0;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            text-align: center;
            margin-bottom: 25px;
        }}
        .main-header a {{
            background-color: #1a2a3f;
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 20px;
            display: block;
            font-size: 1.4rem;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.8);
        }}
        .main-header a:hover {{
            background-color: #6c84b4;
        }}
        .main-header a.active {{
            background-color: #6c84b4;
        }}
        .main-header a:nth-child(1) {{
            margin-right: 30px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(2) {{
            margin: 0 30px;
            border-bottom-left-radius: 8px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(3) {{
            margin-left: 30px;
            border-bottom-left-radius: 8px;
        }}
        .sub-header {{
            grid-column: 1 / -1;
            background-color: white;  
            display: grid;
            grid-template-columns: 1fr minmax(150px, 500px) 100px minmax(150px, 500px) 100px minmax(150px, 500px) 1fr;
            text-align: center;
            padding: 0;
            align-items: stretch;
            margin-bottom: 10px;
        }}
        .sub-header a {{
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 15px 5px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.25rem;
            background-color: #354e7c;
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 1);
        }}
        .sub-header a:nth-child(1) {{ grid-column: 2; }}
        .sub-header a:nth-child(2) {{ grid-column: 4; }}
        .sub-header a:nth-child(3) {{ grid-column: 6; }}
        .sub-header a:hover {{
            background-color: #6c84b4;
        }}
        .sub-header a.active {{
            background-color: #6c84b4;
        }}
        .left-panel {{
            grid-column: 1;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .center-panel {{
            grid-column: 2;
            padding: 0px;
            overflow-y: auto;
        }}
        .right-panel {{
            grid-column: 3;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .table {{
            width: 100%;
            border-collapse: collapse;
        }}
        .table th, .table td {{
            padding: 10px;
            border: 1px solid #ddd;
            text-align: center;
        }}
        .table th {{
            background-color: #a0c4ff;
            color: black;
            font-weight: bold;
        }}
        .table td {{
            white-space: normal !important;
            word-wrap: break-word;
        }}
        .results-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0;
            padding: 10px;
            background-color: white;
            position: sticky;
            top: 0;
            z-index: 1000;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        .table-container {{
            margin-top: 0;
        }}
        .metric-container {{
            margin-top: 20px;
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric-container h3 {{
            color: #1a2a4f;
            border-bottom: 2px solid #a0c4ff;
            padding-bottom: 10px;
        }}
        .metric-container p {{
            margin-bottom: 15px;
            font-size: 1.1rem;
            line-height: 1.6;
        }}
        .filter-dropdown {{
            margin-top: 10px;
        }}
        .tooltip-icon {{
            display: inline-block;
            width: 18px;
            height: 18px;
            background-color: #1a2a4f;
            color: white;
            border-radius: 50%;
            text-align: center;
            font-size: 12px;
            line-height: 18px;
            cursor: help;
            margin-left: 5px;
        }}
        .tooltip-text {{
            visibility: hidden;
            width: 220px;
            background-color: white;
            color: #333;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 12px;
            position: absolute;
            z-index: 1000;
            box-shadow: 0 3px 12px rgba(0, 0, 0, 0.15);
            opacity: 0;
            transition: opacity 0.2s;
            font-size: 14px;
            text-align: left;
        }}
        th:hover .tooltip-text {{
            visibility: visible;
            opacity: 1;
        }}
        .ui-autocomplete {{
            max-height: 200px;
            overflow-y: auto;
            overflow-x: hidden;
            background-color: white;
            border: 1px solid #ddd;
            padding: 5px;
        }}
        .ui-menu-item {{
            padding: 5px;
        }}
        .ui-menu-item:hover {{
            background-color: #f0f0f0;
            cursor: pointer;
        }}
        .run-example-btn {{
            background-color: #e3a42b;  
            color: #ffffff;
            margin-left: auto;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .run-example-btn:hover {{
            background-color: #f9aa33;
            color: white;
        }}
    </style>
</head>
<body>
    <div class="main-header">
        <a href="/gene_symbol" class="active">Gene symbol</a>
        <a href="/gene_ontology_terms">Gene ontology terms</a>
        <a href="/cell_type">Cell type</a>
    </div>
    <div class="sub-header">
        <a href="/gene_relevance">Gene relevance</a>
        <a href="/gene_functions" class="active">Gene functions</a>
        <a href="/new_gene_functions">New gene functions</a>
    </div>
    <div class="left-panel">
        <form method="POST">
            <div class="mb-3">
                <label for="gene_name" class="form-label">Search for gene symbol:</label>
                <input type="text" class="form-control" id="gene_name" name="gene_name" 
                    value="{search_term}" placeholder="e.g., SNCA, AATF" 
                    autocomplete="off" required>
            </div>
            <div class="mb-3">
                <label for="data_source" class="form-label">Database:</label>
                <select class="form-select" id="data_source" name="data_source">
                    <option value="scRNA" {"selected" if data_source == 'scRNA' else ""}>scRNA-seq PD scCoExpNets</option>
                    <option value="bulk" {"selected" if data_source == 'bulk' else ""}>bulk RNA-seq AD TGCNs</option>
                </select>
            </div>
            <div class="mb-3">
                <button type="button" class="btn btn-success" onclick="toggleFilters()">Filters</button>
                <div id="filter-dropdown" class="filter-dropdown">
                    {filter_dropdown}
                </div>
            </div>
            <div class="d-flex justify-content-between align-items-center">
                <button type="submit" class="btn btn-primary">Search</button>
                <button type="button" class="btn run-example-btn" onclick="runExample()">Run example</button>
            </div>
        </form>
    </div>
    <div class="center-panel">
        <div class="results-header">
            <h2>Functional annotations containing {search_term}</h2>
            <form method="POST" action="/download_gene_functions">
                <input type="hidden" name="data_source" value="{data_source}">
                <input type="hidden" name="gene_name" value="{search_term}">
                <input type="hidden" name="cell_type_filter" value="{cell_type_filter if data_source == 'scRNA' else ''}">
                <input type="hidden" name="iteration_filter" value="{iteration_filter if data_source == 'scRNA' else ''}">
                <input type="hidden" name="cluster_filter" value="{cluster_filter if data_source == 'scRNA' else ''}">
                <input type="hidden" name="module_filter" value="{module_filter}">
                <input type="hidden" name="target_filter" value="{target_filter if data_source == 'bulk' else ''}">
                <input type="hidden" name="tissue_filter" value="{tissue_filter if data_source == 'bulk' else ''}">
                <input type="hidden" name="cutoff_filter" value="{cutoff_filter if data_source == 'bulk' else ''}">
                <div class="mb-3">
                    <label for="download_format" class="form-label">Choose format:</label>
                    <select class="form-select d-inline-block w-auto" id="download_format" name="download_format">
                        <option value="csv">CSV</option>
                        <option value="xlsx">XLSX</option>
                        <option value="html">HTML</option>
                    </select>
                    <button type="submit" class="btn btn-success">Download</button>
                </div>
            </form>
        </div>
        <div class="table-container">
            <table class="table table-bordered table-hover">
                <thead class="table-primary">
                    <tr>
                        {''.join([
                            f'<th class="text-center">{header} <div class="tooltip-icon">?</div><div class="tooltip-text">{column_descriptions[data_source].get(header, "")}</div></th>'
                            for header in headers
                        ])}
                    </tr>
                </thead>
                <tbody>
                    {table_rows if search_term else f"<tr><td colspan='{len(headers)}' class='text-center'>Enter a gene symbol to search</td></tr>"}
                </tbody>
            </table>
        </div>
    </div>
    <div class="right-panel">
        <div class="metric-container">
            <h3>Functional annotation statistics for {search_term}</h3>
            {"<p>• " + search_term + " is minimally expressed in " + stats["minimally_expressed"] + " out of 24 cell types (" + stats["minimally_expressed_percentage"] + ")</p>" if data_source == 'scRNA' and search_term else ""}
            {"<p>• " + search_term + " has a mean module membership of " + stats["mean_mm"] + " (percentile: " + stats["mean_percentile"] + ") per cell type</p>" if data_source == 'scRNA' and search_term else ""}
            {"<p>• Participates in " + stats["participation_percentage"] + "% of its module's known annotations</p>" if data_source == 'scRNA' and search_term else ""}
            {"<p>• " + search_term + " is involved in " + stats["mean_annotations"] + " annotations on average per cell type, with a mean IC of " + stats["mean_ic"] + "</p>" if data_source == 'scRNA' and search_term else ""}
        </div>
    </div>
    <script>
        function toggleFilters() {{
            const filterDropdown = document.getElementById('filter-dropdown');
            filterDropdown.style.display = filterDropdown.style.display === 'none' ? 'block' : 'none';
        }}
        
        document.getElementById('data_source').addEventListener('change', function() {{
            this.form.submit();
        }});
    </script>
    <script src="/static/sort_table.js"></script>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://code.jquery.com/ui/1.13.1/jquery-ui.min.js"></script>
    <script>
    $(function() {{
        $("#gene_name").autocomplete({{
            source: function(request, response) {{
                $.getJSON("/api/genes", {{
                    term: request.term
                }}, response);
            }},
            minLength: 2,
            select: function(event, ui) {{
                $("#gene_name").val(ui.item.value);
                $(this).closest("form").submit();
            }},
            focus: function(event, ui) {{
                $("#gene_name").val(ui.item.value);
                return false;
            }}
        }});
    }});
    </script>
    {f'<script src="/static/update_filters.js"></script>' if data_source == 'scRNA' else ''}
    <script>
    function runExample() {{
        const form = document.querySelector('.left-panel form');
        form.gene_name.value = "SNCA";
        form.data_source.value = "scRNA";
        document.getElementById('cell_type_filter').value = "Oligodendrocytes";
        document.getElementById('module_filter').value = "turquoise";
        form.submit();
    }}
    </script>
</body>
</html>
""")


@app.route('/exclusive_relevant_genes', methods=['GET', 'POST'])
def exclusive_relevant_genes():
    if request.method == 'POST':
        cell_type_filter = request.form.get('cell_type_filter', '').strip()
    else:
        cell_type_filter = request.args.get('cell_type_filter', '').strip()

    unique_cell_types = set()
    for df in [minimally_expressed_df, relevant_at_t0_df, relevant_in_all_iterations_df]:
        cell_types = df["Cell Types"].dropna().astype(str).str.split("; ").explode()
        cell_types = cell_types.str.replace("_", " ")
        unique_cell_types.update(cell_types)

    minimally_expressed_genes_all = set()
    minimally_expressed_genes_only = set()
    minimally_expressed_percentage = 0.0

    relevant_at_t0_genes_all = set()
    relevant_at_t0_genes_only = set()
    relevant_at_t0_percentage = 0.0

    relevant_in_all_iterations_genes_all = set()
    relevant_in_all_iterations_genes_only = set()
    relevant_in_all_iterations_percentage = 0.0

    genes_table = []

    if cell_type_filter:
        normalized_filter = cell_type_filter.replace("_", " ")
        
        minimally_expressed_genes_all = set(minimally_expressed_df[minimally_expressed_df["Cell Types"].str.replace("_", " ").str.contains(normalized_filter, na=False)]["Gene"])
        relevant_at_t0_genes_all = set(relevant_at_t0_df[relevant_at_t0_df["Cell Types"].str.replace("_", " ").str.contains(normalized_filter, na=False)]["Gene"])
        relevant_in_all_iterations_genes_all = set(relevant_in_all_iterations_df[relevant_in_all_iterations_df["Cell Types"].str.replace("_", " ").str.contains(normalized_filter, na=False)]["Gene"])

        minimally_expressed_genes_only = set()
        for gene in minimally_expressed_genes_all:
            cell_types_for_gene = minimally_expressed_df[minimally_expressed_df["Gene"] == gene]["Cell Types"].str.replace("_", " ").str.split("; ").explode()
            if all(cell_type == normalized_filter for cell_type in cell_types_for_gene):
                minimally_expressed_genes_only.add(gene)

        relevant_at_t0_genes_only = set()
        for gene in relevant_at_t0_genes_all:
            cell_types_for_gene = relevant_at_t0_df[relevant_at_t0_df["Gene"] == gene]["Cell Types"].str.replace("_", " ").str.split("; ").explode()
            if all(cell_type == normalized_filter for cell_type in cell_types_for_gene):
                relevant_at_t0_genes_only.add(gene)

        relevant_in_all_iterations_genes_only = set()
        for gene in relevant_in_all_iterations_genes_all:
            cell_types_for_gene = relevant_in_all_iterations_df[relevant_in_all_iterations_df["Gene"] == gene]["Cell Types"].str.replace("_", " ").str.split("; ").explode()
            if all(cell_type == normalized_filter for cell_type in cell_types_for_gene):
                relevant_in_all_iterations_genes_only.add(gene)

        minimally_expressed_percentage = (len(minimally_expressed_genes_only) / len(minimally_expressed_genes_all)) * 100 if len(minimally_expressed_genes_all) > 0 else 0
        relevant_at_t0_percentage = (len(relevant_at_t0_genes_only) / len(relevant_at_t0_genes_all)) * 100 if len(relevant_at_t0_genes_all) > 0 else 0
        relevant_in_all_iterations_percentage = (len(relevant_in_all_iterations_genes_only) / len(relevant_in_all_iterations_genes_all)) * 100 if len(relevant_in_all_iterations_genes_all) > 0 else 0

    cell_type_name = None
    cluster_number = None
    if cell_type_filter:
        match = re.match(r"(.*?)\s*(\d+)$", cell_type_filter)
        if match:
            cell_type_name = match.group(1).strip()  
            cluster_number = int(match.group(2))   

    statistics_list = [
        f"• {len(minimally_expressed_genes_all)} minimally expressed genes, {len(minimally_expressed_genes_only)} of them only expressed in {cell_type_filter if cell_type_filter else 'this cell type'} ({minimally_expressed_percentage:.2f}%)",
        f"• {len(relevant_at_t0_genes_all)} relevant genes at T0, {len(relevant_at_t0_genes_only)} of them only relevant in {cell_type_filter if cell_type_filter else 'this cell type'} ({relevant_at_t0_percentage:.2f}%)",
        f"• {len(relevant_in_all_iterations_genes_all)} relevant genes in all iterations, {len(relevant_in_all_iterations_genes_only)} of them only relevant in {cell_type_filter if cell_type_filter else 'this cell type'} ({relevant_in_all_iterations_percentage:.2f}%)"
    ]

    column_descriptions = {
        "Cell type": "Cell type and cluster where the gene is a potential marker",
        "Cluster": "Subgroup identifier within the cell type",
        "Criteria": "Type of relevance: minimally expressed, relevant at T0, or in all iterations",
        "Gene Name": "Gene symbol identified as a potential marker"
    }
    headers = ["Cell type", "Cluster", "Criteria", "Gene Name"]
    headers_with_tooltips = []
    for header in headers:
        description = column_descriptions.get(header, "No description available")
        headers_with_tooltips.append(
            f'''<th class="text-center">
                {header}
                <span class="tooltip-icon">?
                    <span class="tooltip-text">{description}</span>
                </span>
            </th>'''
        )

    return render_template_string(f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Cell Type Statistics</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: grid;
            grid-template-columns: 425px 1fr 475px;
            grid-template-rows: auto auto 1fr;
            height: 100vh;
        }}
        .main-header {{
            grid-column: 1 / -1;
            background-color: white;
            color: white;
            padding: 0;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            text-align: center;
            margin-bottom: 25px;
        }}
        .main-header a {{
            background-color: #1a2a3f;
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 20px;
            display: block;
            font-size: 1.4rem;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.8);
        }}
        .main-header a:hover {{
            background-color: #6c84b4;
        }}
        .main-header a.active {{
            background-color: #6c84b4;
        }}
        .main-header a:nth-child(1) {{
            margin-right: 30px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(2) {{
            margin: 0 30px;
            border-bottom-left-radius: 8px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(3) {{
            margin-left: 30px;
            border-bottom-left-radius: 8px;
        }}
        .sub-header {{
            grid-column: 1 / -1;
            background-color: white;  
            display: grid;
            grid-template-columns: 1fr minmax(200px, 500px)  100px minmax(200px, 500px) 1fr;
            text-align: center;
            padding: 0;
            align-items: stretch;
            margin-bottom: 10px;
        }}
        .sub-header a {{
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 15px 5px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.25rem;
            background-color: #354e7c;
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 1);
        }}
        .sub-header a:nth-child(1) {{ grid-column: 2; }}
        .sub-header a:nth-child(2) {{ grid-column: 4; }}
        .sub-header a:nth-child(3) {{ grid-column: 6; }}
        .sub-header a:hover {{
            background-color: #6c84b4;
        }}
        .sub-header a.active {{
            background-color: #6c84b4;
        }}
        .left-panel {{
            grid-column: 1;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .center-panel {{
            grid-column: 2;
            overflow-y: auto;
        }}
        .right-panel {{
            grid-column: 3;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .metric-container {{
            margin-top: 20px;
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric-container h3 {{
            color: #1a2a4f;
            border-bottom: 2px solid #a0c4ff;
            padding-bottom: 10px;
        }}
        .metric-container p {{
            margin-bottom: 15px;
            font-size: 1.1rem;
            line-height: 1.6;
            cursor: pointer;  
        }}
        .metric-container p:hover {{
            text-decoration: underline;
            color: #1a2a4f; 
        }}
        .table {{
            width: 100%;
            border-collapse: collapse;
        }}
        .table th, .table td {{
            padding: 10px;
            border: 1px solid #ddd;
            text-align: center;
        }}
        .table th {{
            background-color: #a0c4ff;
            color: black;
            font-weight: bold;
            cursor: pointer;
        }}
        .table th:hover {{
            background-color: #89b4ff;
        }}
        .results-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-top: 0;
            padding: 10px;
            background-color: white;
            position: sticky;
            top: 0;
            z-index: 1000;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        .tooltip-icon {{
            display: inline-block;
            width: 18px;
            height: 18px;
            background-color: #1a2a4f;
            color: white;
            border-radius: 50%;
            text-align: center;
            font-size: 12px;
            line-height: 18px;
            cursor: help;
            margin-left: 5px;
        }}
        .tooltip-text {{
            visibility: hidden;
            width: 220px;
            background-color: white;
            color: #333;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 12px;
            position: absolute;
            z-index: 1000;
            box-shadow: 0 3px 12px rgba(0, 0, 0, 0.15);
            opacity: 0;
            transition: opacity 0.2s;
            font-size: 14px;
            text-align: left;
            line-height: 1.5;
        }}
        th:hover .tooltip-text {{
            visibility: visible;
            opacity: 1;
        }}
        .run-example-btn {{
            background-color: #fcbf49;
            color: #ffffff;
            font-weight: bold;
            margin-left: auto;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .run-example-btn:hover {{
            background-color: #f9aa33;
            color: white;
        }}
    </style>
</head>
<body>
    <div class="main-header">
        <a href="/gene_symbol">Gene symbol</a>
        <a href="/gene_ontology_terms">Gene ontology terms</a>
        <a href="/cell_type" class="active">Cell type</a>
    </div>
    <div class="sub-header">
        <a href="/exclusive_relevant_genes" class="active">Exclusive relevant genes</a>
        <a href="/exclusive_go_terms">Exclusive GO terms</a>
    </div>
    <div class="left-panel">
        <form method="POST">
            <div class="mb-3">
                <label for="cell_type_filter" class="form-label">Select a cell type:</label>
                <select class="form-select" id="cell_type_filter" name="cell_type_filter">
                    <option value="">Select a cell type</option>
                    {"".join([f'<option value="{cell_type}" {"selected" if cell_type == cell_type_filter else ""}>{cell_type}</option>' for cell_type in sorted(unique_cell_types)])}
                </select>
            </div>
            <div class="mb-3">
                <label for="data_source" class="form-label">Database:</label>
                <select class="form-select" id="data_source" name="data_source" disabled>
                    <option value="scRNA" selected>scRNA-seq PD scCoExpNets</option>
                </select>
            </div>
            <div class="d-flex justify-content-between">
                <button type="submit" class="btn btn-primary">Search</button>
                <button type="button" class="btn run-example-btn" onclick="runExample()">Run example</button>
            </div>
        </form>
    </div>
    <div class="center-panel">
        <div class="results-header">
            <h2>Potential biomarkers exclusively found in {cell_type_filter}</h2>
            <form method="POST" action="/download_exclusive_genes">
                <input type="hidden" name="cell_type_filter" value="{cell_type_filter}">
                <div class="mb-3">
                    <label for="download_format" class="form-label">Choose format:</label>
                    <select class="form-select d-inline-block w-auto" id="download_format" name="download_format">
                        <option value="csv">CSV</option>
                        <option value="xlsx">XLSX</option>
                        <option value="html">HTML</option>
                    </select>
                    <button type="submit" class="btn btn-success">Download</button>
                </div>
            </form>
        </div>
        <table class="table table-bordered table-hover" id="sortable-table">
            <thead class="table-primary">
                <tr>
                    <tr>
                        {''.join(headers_with_tooltips)}
                    </tr>
                </tr>
            </thead>
            <tbody id="genes-table-body">
            </tbody>
        </table>
    </div>
    <div class="right-panel">
        <div class="metric-container">
            <h3>Gene relevance distribution for {cell_type_filter}</h3>
            {f"""
            <p onclick="loadGenes('minimally_expressed')">• {len(minimally_expressed_genes_all)} minimally expressed genes, {len(minimally_expressed_genes_only)} of them only expressed in {cell_type_filter} ({minimally_expressed_percentage:.2f}%)</p>
            <p onclick="loadGenes('relevant_at_t0')">• {len(relevant_at_t0_genes_all)} relevant genes at T0, {len(relevant_at_t0_genes_only)} of them only relevant in {cell_type_filter} ({relevant_at_t0_percentage:.2f}%)</p>
            <p onclick="loadGenes('relevant_in_all_iterations')">• {len(relevant_in_all_iterations_genes_all)} relevant genes in all iterations, {len(relevant_in_all_iterations_genes_only)} of them only relevant in {cell_type_filter} ({relevant_in_all_iterations_percentage:.2f}%)</p>
            """ if cell_type_filter else "<p>Please select a cell type to view statistics</p>"}
        </div>
    </div>
    <script>
        function loadGenes(criteria) {{
            const cellTypeFilter = "{cell_type_filter}";
            let genes = [];
            
            // Simular la carga de datos según el criterio seleccionado
            if (criteria === 'minimally_expressed') {{
                genes = {list(minimally_expressed_genes_only)};
            }} else if (criteria === 'relevant_at_t0') {{
                genes = {list(relevant_at_t0_genes_only)};
            }} else if (criteria === 'relevant_in_all_iterations') {{
                genes = {list(relevant_in_all_iterations_genes_only)};
            }}

            // Construir la tabla con los genes
            const tableBody = document.getElementById('genes-table-body');
            tableBody.innerHTML = ''; // Limpiar la tabla

            genes.forEach(gene => {{
                const row = document.createElement('tr');
                row.innerHTML = `
                    <td class="text-center">{cell_type_name if cell_type_name else "N/A"}</td>
                    <td class="text-center">{cluster_number if cluster_number is not None else "N/A"}</td>
                    <td class="text-center">${{criteria.replace(/_/g, ' ')}}</td>
                    <td class="text-center">${{gene}}</td>
                `;
                tableBody.appendChild(row);
            }});
        }}
    </script>
    <script src="/static/sort_table.js"></script>
    <script>
    function runExample() {{
        const form = document.querySelector('.left-panel form');
        const selector = document.getElementById('cell_type_filter');
        selector.value = "T cells 18";
        form.submit();
    }}

    // Al volver del POST, si se ha seleccionado un cell type, cargar tabla automáticamente
    document.addEventListener("DOMContentLoaded", function () {{
        if ("{{ cell_type_filter }}") {{
            loadGenes("minimally_expressed");
        }}
    }});
    </script>

</body>
</html>
""")

@app.route('/go_term_relevance', methods=['GET', 'POST'])
def go_term_relevance():
    data_source = request.form.get('data_source', 'scRNA').strip()
    search_term = request.form.get('search_term', '').strip()
    cell_type_filter = request.form.get('cell_type_filter', '').strip()
    iteration_filter = request.form.get('iteration_filter', '').strip()
    cluster_filter = request.form.get('cluster_filter', '').strip()
    module_filter = request.form.get('module_filter', '').strip()
    target_filter = request.form.get('target_filter', '').strip()
    tissue_filter = request.form.get('tissue_filter', '').strip()
    cutoff_filter = request.form.get('cutoff_filter', '').strip()

    if data_source == 'scRNA':
        results = query_annotations(search_term, cell_type_filter, iteration_filter, cluster_filter, module_filter)

        if results:
            annotations_data = load_csv(annotations_file)
            if annotations_data is not None:
                ic_values = pd.to_numeric(annotations_data['IC'], errors='coerce')
                max_ic = ic_values[ic_values != float('inf')].max()
                replacement_ic = max_ic + 10 if not pd.isna(max_ic) else 0

            term_ic = format_ic(results[0].get('IC', 0)) if results else "N/A"

            unique_cell_types = set(f"{result['Cell type']} {result['Cluster']}" for result in results)
            cell_types_count = len(unique_cell_types)
            cell_types_percentage = (cell_types_count / 24) * 100

            subgraph_sizes = [result["Subgraph size"] for result in results if "Subgraph size" in result]
            mean_subgraph_size = int(round(sum(subgraph_sizes) / len(subgraph_sizes))) if subgraph_sizes else "N/A"

            p_values = [float(result["P-value"]) for result in results if "P-value" in result and result["P-value"]]
            mean_neg_log_pvalue = sum(-np.log10(p_values)) / len(p_values) if p_values else "N/A"

            ic_values = []
            for result in results:
                ic = result.get('IC', 0)
                try:
                    ic = float(ic)
                    if ic == float('inf'):
                        ic = replacement_ic
                    ic_values.append(ic)
                except (ValueError, TypeError):
                    continue
            
            mean_ic = round(sum(ic_values) / len(ic_values), 2) if ic_values else "N/A"

            stats = {
                "term_ic": term_ic,
                "cell_types_count": cell_types_count,
                "cell_types_percentage": round(cell_types_percentage, 2),
                "mean_subgraph_size": mean_subgraph_size,
                "mean_neg_log_pvalue": round(mean_neg_log_pvalue, 2) if mean_neg_log_pvalue != "N/A" else "N/A",
                "mean_ic": mean_ic
            }
        else:
            stats = {
                "term_ic": "N/A",
                "cell_types_count": 0,
                "cell_types_percentage": 0.0,
                "mean_subgraph_size": "N/A",
                "mean_neg_log_pvalue": "N/A",
                "mean_ic": "N/A"
            }

        headers = ["Iteration", "Cell type", "Cluster", "Module", "Term id", "Term name", 
                   "P-value", "Intersection", "Length of Intersection", "Source", 
                   "Subgraph ID", "Subgraph size", "IC"] 
        table_rows = ""
        for result in results:
            p_value_formatted = format_p_value(result.get("P-value", ""))
            intersection_formatted = format_intersection(result.get("Intersection", ""))
            
            row = f"""
            <tr>
                <td class='text-center'>{result.get('Iteration', '')}</td>
                <td class='text-center'>{result.get('Cell type', '')}</td>
                <td class='text-center'>{result.get('Cluster', '')}</td>
                <td class='text-center'>{result.get('Module', '')}</td>
                <td class='text-center'>{result.get('Term id', '')}</td>
                <td class='text-center'>{result.get('Term name', '')}</td>
                <td class='text-center'>{p_value_formatted}</td>
                <td class='text-center'>{intersection_formatted}</td>
                <td class='text-center'>{result.get('Length of Intersection', '')}</td>
                <td class='text-center'>{result.get('Source', '')}</td>
                <td class='text-center'>{result.get('Subgraph ID', '')}</td>
                <td class='text-center'>{result.get('Subgraph size', '')}</td>
                <td class='text-center'>{result.get('IC', '')}</td>
                <td class='text-center'>
                    <button class='btn btn-sm btn-outline-primary toggle-subgraph' data-subgraph='{result.get("Subgraph ID", "")}'>+</button>
                </td>
            </tr>
            """
            table_rows += row

        full_data = load_csv(annotations_file)
        full_data_json = full_data.to_json(orient='records') if full_data is not None else "[]"

        cell_types = set()
        for file_name in os.listdir(modules_dir):
            if file_name.endswith(".csv"):
                cell_type = extract_cell_type(file_name)
                if cell_type:
                    cell_types.add(cell_type)

        filter_dropdown = f"""
            <div class="mb-3">
                <label for="cell_type_filter" class="form-label">Cell type:</label>
                <select class="form-select" id="cell_type_filter" name="cell_type_filter" onchange="updateClusters()">
                    <option value="">All cell types</option>
                    {"".join([f'<option value="{cell_type}" {"selected" if cell_type == cell_type_filter else ""}>{cell_type}</option>' 
                    for cell_type in sorted(cell_types, key=lambda x: (x.split()[0], int(x.split()[1]) if len(x.split()) > 1 and x.split()[1].isdigit() else (x, 0)))])}
                </select>
            </div>
            <div class="mb-3">
                <label for="cluster_filter" class="form-label">Cluster:</label>
                <select class="form-select" id="cluster_filter" name="cluster_filter" onchange="updateIterations()">
                    <option value="">All clusters</option>
                    {"".join([f'<option value="{cluster}" {"selected" if str(cluster) == cluster_filter else ""}>{cluster}</option>' 
                    for cluster in get_available_clusters(cell_type_filter) if cell_type_filter])}
                </select>
            </div>
            <div class="mb-3">
                <label for="iteration_filter" class="form-label">Iteration:</label>
                <select class="form-select" id="iteration_filter" name="iteration_filter">
                    <option value="">All iterations</option>
                    {"".join([f'<option value="{iteration}" {"selected" if iteration == iteration_filter else ""}>{iteration}</option>' 
                    for iteration in get_available_iterations(cell_type_filter, cluster_filter if cluster_filter else None) if cell_type_filter])}
                </select>
            </div>
            <div class="mb-3">
                <label for="module_filter" class="form-label">Module:</label>
                <input type="text" class="form-control" id="module_filter" name="module_filter" value="{module_filter}" placeholder="e.g., red, blue...">
            </div>
        """
        
    else:
        bulk_annotations = pd.read_csv(bulk_annotations_file)
        
        results = []
        if search_term:
            filtered_data = bulk_annotations[
                (bulk_annotations['term_id'].str.lower() == search_term.lower()) | 
                (bulk_annotations['term_name'].str.lower() == search_term.lower())
            ]
            
            if target_filter:
                filtered_data = filtered_data[filtered_data['target'] == target_filter]
            if tissue_filter:
                filtered_data = filtered_data[filtered_data['tissue'] == tissue_filter]
            if cutoff_filter:
                filtered_data = filtered_data[filtered_data['cutoff'].astype(str) == cutoff_filter]
            if module_filter:
                filtered_data = filtered_data[filtered_data['module'] == module_filter]
            
            filtered_data = filtered_data.sort_values(by='p_value')
            
            results = filtered_data.to_dict('records')

        headers = ["Cutoff", "Target", "Tissue", "Phenotype", "Module", "Term id", 
                   "Term name", "P-value", "Intersection", "Length of Intersection", 
                   "Source", "IC"]
        
        table_rows = ""
        for result in results:
            p_value_formatted = format_p_value(result.get("p_value", ""))
            intersection_formatted = format_intersection(result.get("intersection", ""))
            ic_formatted = format_ic(result.get("IC", ""))
            
            row = f"""
            <tr>
                <td class='text-center'>{result.get('cutoff', '')}</td>
                <td class='text-center'>{result.get('target', '')}</td>
                <td class='text-center'>{result.get('tissue', '')}</td>
                <td class='text-center'>{result.get('phenotype', '')}</td>
                <td class='text-center'>{result.get('module', '')}</td>
                <td class='text-center'>{result.get('term_id', '')}</td>
                <td class='text-center'>{result.get('term_name', '')}</td>
                <td class='text-center'>{p_value_formatted}</td>
                <td class='text-center'>{intersection_formatted}</td>
                <td class='text-center'>{result.get('length_intersection', '')}</td>
                <td class='text-center'>{result.get('source', '')}</td>
                <td class='text-center'>{ic_formatted}</td>
            </tr>
            """
            table_rows += row

        available_targets = sorted(bulk_annotations['target'].unique())
        available_tissues = sorted(bulk_annotations['tissue'].unique())
        available_cutoffs = sorted(bulk_annotations['cutoff'].unique())
        
        filter_dropdown = f"""
            <div class="mb-3">
                <label for="target_filter" class="form-label">Target:</label>
                <select class="form-select" id="target_filter" name="target_filter">
                    <option value="">All targets</option>
                    {"".join([f'<option value="{target}" {"selected" if target == target_filter else ""}>{target}</option>' 
                    for target in available_targets])}
                </select>
            </div>
            <div class="mb-3">
                <label for="tissue_filter" class="form-label">Tissue:</label>
                <select class="form-select" id="tissue_filter" name="tissue_filter">
                    <option value="">All tissues</option>
                    {"".join([f'<option value="{tissue}" {"selected" if tissue == tissue_filter else ""}>{tissue}</option>' 
                    for tissue in available_tissues])}
                </select>
            </div>
            <div class="mb-3">
                <label for="cutoff_filter" class="form-label">Cutoff:</label>
                <select class="form-select" id="cutoff_filter" name="cutoff_filter">
                    <option value="">All cutoffs</option>
                    {"".join([f'<option value="{cutoff}" {"selected" if str(cutoff) == cutoff_filter else ""}>{cutoff}</option>' 
                    for cutoff in available_cutoffs])}
                </select>
            </div>
            <div class="mb-3">
                <label for="module_filter" class="form-label">Module:</label>
                <input type="text" class="form-control" id="module_filter" name="module_filter" 
                    value="{module_filter}" placeholder="e.g., PPP2CA">
            </div>
        """
        
        stats = {}
        full_data_json = "[]"

    column_descriptions = {
        "scRNA": {
            "Iteration": "Number of times the pseudo-cell algorithm has been executed (e.g., T0, T5, etc.)",
            "Cell type": "Cell type where the gene is expressed",
            "Cluster": "Subgroup within the cell type",
            "Module": "Name of the module",
            "Term id": "Identifier of an annotation in the GO database",
            "Term name": "Full name of a GO annotation",
            "P-value": "Significance of an annotation given a set of genes in a co-expression module (lower = more significant)",
            "Intersection": "Name of the genes in a module that intersect with a GO database annotation",
            "Length of Intersection": "Number of genes in a module that intersect with a GO database annotation",
            "Source": "Database used to obtain the annotations (biological_process, cellular_component and molecular_function)",
            "Subgraph ID": "GO annotation subgraph identifier",
            "Subgraph size": "Number of GO annotations that make up a subgraph",
            "IC": "Index that tells us how informative an annotation is. The higher the CI, the more specific and informative the annotation is."
        },
        "bulk": {
            "Cutoff": "Indicates the stringency level, ranging from 1 to 10 where higher values correspond to more reliable modules",
            "Target": "Main variable of interest, typically a gene or a covariate associated with each individual",
            "Tissue": "Tissue sample origin",
            "Phenotype": "Disease or condition under study",
            "Module": "Name of the co-expression module",
            "Term id": "Identifier of an annotation in the GO database",
            "Term name": "Full name of a GO annotation",
            "P-value": "Significance of an annotation given a set of genes in a co-expression module (lower = more significant)",
            "Intersection": "Name of the genes in a module that intersect with a GO database annotation",
            "Length of Intersection": "Number of genes in a module that intersect with a GO database annotation",
            "Source": "Database used to obtain the annotations (biological_process, cellular_component and molecular_function)",
            "IC": "Index that tells us how informative an annotation is. The higher the CI, the more specific and informative the annotation is."
        }
    }

    current_descriptions = column_descriptions['scRNA'] if data_source == 'scRNA' else column_descriptions['bulk']
    headers_with_tooltips = []
    for header in headers:
        description = current_descriptions.get(header, "No description available")
        headers_with_tooltips.append(
            f'''
            <th class="text-center">
                {header}
                <div class="tooltip-icon">?</div>
                <div class="tooltip-text">{description}</div>
            </th>
            '''
        )

    return render_template_string(f"""   
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GO Term Search</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: grid;
            grid-template-columns: 425px 1fr 475px;
            grid-template-rows: auto auto 1fr;
            height: 100vh;
        }}
        .main-header {{
            grid-column: 1 / -1;
            background-color: white;
            color: white;
            padding: 0;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            text-align: center;
            margin-bottom: 25px;
        }}
        .main-header a {{
            background-color: #1a2a3f;
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 20px;
            display: block;
            font-size: 1.4rem;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.8);
        }}
        .main-header a:hover {{
            background-color: #6c84b4;
        }}
        .main-header a.active {{
            background-color: #6c84b4;
        }}
        .main-header a:nth-child(1) {{
            margin-right: 30px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(2) {{
            margin: 0 30px;
            border-bottom-left-radius: 8px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(3) {{
            margin-left: 30px;
            border-bottom-left-radius: 8px;
        }}
        .sub-header {{
            grid-column: 1 / -1;
            background-color: white;  
            display: grid;
            grid-template-columns: 1fr minmax(250px, 500px)  1fr;
            text-align: center;
            padding: 0;
            align-items: stretch;
            margin-bottom: 10px;
        }}
        .sub-header a {{
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 15px 5px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.25rem;
            background-color: #354e7c;
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 1);
        }}
        .sub-header a:nth-child(1) {{ grid-column: 2; }}
        .sub-header a:nth-child(2) {{ grid-column: 4; }}
        .sub-header a:nth-child(3) {{ grid-column: 6; }}
        .sub-header a:hover {{
            background-color: #6c84b4;
        }}
        .sub-header a.active {{
            background-color: #6c84b4;
        }}
        .left-panel {{
            grid-column: 1;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .center-panel {{
            grid-column: 2;
            padding: 0px;
            overflow-y: auto;
        }}
        .right-panel {{
            grid-column: 3;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .table {{
            width: 100%;
            border-collapse: collapse;
        }}
        .table th, .table td {{
            padding: 10px;
            border: 1px solid #ddd;
            text-align: center;
        }}
        .table th {{
            background-color: #a0c4ff;
            color: black;
            font-weight: bold;
        }}
        .table td {{
            white-space: normal !important;
            word-wrap: break-word;
        }}
        .results-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0;
            padding: 10px;
            background-color: white;
            position: sticky;
            top: 0;
            z-index: 1000;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        .table-container {{
            margin-top: 0;
        }}
        .metric-container {{
            margin-top: 20px;
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric-container h3 {{
            color: #1a2a4f;
            border-bottom: 2px solid #a0c4ff;
            padding-bottom: 10px;
        }}
        .metric-container p {{
            margin-bottom: 15px;
            font-size: 1.1rem;
            line-height: 1.6;
        }}
        .filter-dropdown {{
            margin-top: 10px;
        }}
        .additional-table th {{
            background-color: #c3ebc3;
            color: black;
            font-weight: bold;
        }}
        .tooltip-icon {{
            display: inline-block;
            width: 18px;
            height: 18px;
            background-color: #1a2a4f;
            color: white;
            border-radius: 50%;
            text-align: center;
            font-size: 12px;
            line-height: 18px;
            cursor: help;
            margin-left: 5px;
        }}
        .tooltip-text {{
            visibility: hidden;
            width: 220px;
            background-color: white;
            color: #333;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 12px;
            position: absolute;
            z-index: 1000;
            box-shadow: 0 3px 12px rgba(0, 0, 0, 0.15);
            opacity: 0;
            transition: opacity 0.2s;
            font-size: 14px;
            text-align: left;
        }}
        th:hover .tooltip-text {{
            visibility: visible;
            opacity: 1;
        }}
        .ui-autocomplete {{
            max-height: 200px;
            overflow-y: auto;
            overflow-x: hidden;
            background-color: white;
            border: 1px solid #ddd;
            padding: 5px;
        }}
        .ui-menu-item {{
            padding: 5px;
        }}
        .ui-menu-item:hover {{
            background-color: #f0f0f0;
            cursor: pointer;
        }}
        .run-example-btn {{
            background-color: #e3a42b;  
            color: #ffffff;
            margin-left: auto;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .run-example-btn:hover {{
            background-color: #f9aa33;
            color: white;
        }}
    </style>
</head>
<body>
    <div class="main-header">
        <a href="/gene_symbol">Gene symbol</a>
        <a href="/gene_ontology_terms" class="active">Gene ontology terms</a>
        <a href="/cell_type">Cell type</a>
    </div>
    <div class="sub-header">
        <a href="/go_term_relevance" class="active">GO term relevance</a>
    </div>
    <div class="left-panel">
        <form method="POST">
            <div class="mb-3">
                <label for="search_term" class="form-label">Search for GO ID or term name:</label>
                <input type="text" class="form-control ui-autocomplete-input" 
                       id="search_term" name="search_term" 
                       value="{search_term}" 
                       placeholder="e.g., GO:0043229 or intracellular organelle"
                       autocomplete="off">
            </div>
            <div class="mb-3">
                <label for="data_source" class="form-label">Database:</label>
                <select class="form-select" id="data_source" name="data_source">
                    <option value="scRNA" {"selected" if data_source == 'scRNA' else ""}>scRNA-seq PD scCoExpNets</option>
                    <option value="bulk" {"selected" if data_source == 'bulk' else ""}>bulk RNA-seq AD TGCNs</option>
                </select>
            </div>
            <div class="mb-3">
                <button type="button" class="btn btn-success" onclick="toggleFilters()">Filters</button>
                <div id="filter-dropdown" class="filter-dropdown">
                    {filter_dropdown}
                </div>
            </div>
            <div class="d-flex justify-content-between align-items-center">
                <button type="submit" class="btn btn-primary">Search</button>
                <button type="button" class="btn run-example-btn" onclick="runExample()">Run example</button>
            </div>
        </form>
    </div>
    <div class="center-panel">
        <div class="results-header">
            <h2>Modules associated with {search_term} biological function</h2>
            <form method="POST" action="/download_go_terms">
                <input type="hidden" name="data_source" value="{data_source}">
                <input type="hidden" name="search_term" value="{search_term}">
                <input type="hidden" name="cell_type_filter" value="{cell_type_filter if data_source == 'scRNA' else ''}">
                <input type="hidden" name="iteration_filter" value="{iteration_filter if data_source == 'scRNA' else ''}">
                <input type="hidden" name="cluster_filter" value="{cluster_filter if data_source == 'scRNA' else ''}">
                <input type="hidden" name="module_filter" value="{module_filter}">
                <input type="hidden" name="target_filter" value="{target_filter if data_source == 'bulk' else ''}">
                <input type="hidden" name="tissue_filter" value="{tissue_filter if data_source == 'bulk' else ''}">
                <input type="hidden" name="cutoff_filter" value="{cutoff_filter if data_source == 'bulk' else ''}">
                <div class="mb-3">
                    <label for="download_format" class="form-label">Choose format:</label>
                    <select class="form-select d-inline-block w-auto" id="download_format" name="download_format">
                        <option value="csv">CSV</option>
                        <option value="xlsx">XLSX</option>
                        <option value="html">HTML</option>
                    </select>
                    <button type="submit" class="btn btn-success">Download</button>
                </div>
            </form>
        </div>
        <div class="table-container">
            <table class="table table-bordered table-hover">
                <thead class="table-primary">
                    <tr>
                        {''.join(headers_with_tooltips)}
                    </tr>
                </thead>
                <tbody>
                    {table_rows if search_term else "<tr><td colspan='" + str(len(headers)) + "' class='text-center'>Enter a GO term or ID to search</td></tr>"}
                </tbody>
            </table>
        </div>
    </div>
    <div class="right-panel">
        {f"""
        <div class="metric-container">
            <h3>Functional enrichment statistics for {search_term}</h3>
            <p>- IC of the annotation: {stats.get("term_ic", "N/A")}</p>
            <p>- Found in {stats.get("cell_types_count", 0)} of 24 cell types ({stats.get("cell_types_percentage", 0):.2f}%)</p>
            <p>- Mean subgraph size: {stats.get("mean_subgraph_size", "N/A")}</p>
            <p>- Mean –log10(p-value): {stats.get("mean_neg_log_pvalue", "N/A")}</p>
            <p>- Mean IC of subgraphs: {stats.get("mean_ic", "N/A")}</p>
        </div>
        """ if search_term and data_source == 'scRNA' else ""}
    </div>
    <script>
        const fullData = {full_data_json};
        function toggleFilters() {{
            const filterDropdown = document.getElementById('filter-dropdown');
            filterDropdown.style.display = filterDropdown.style.display === 'none' ? 'block' : 'none';
        }}
        
        document.getElementById('data_source').addEventListener('change', function() {{
            this.form.submit();
        }});
    </script>
    <script src="/static/sort_table.js"></script>
    {f'<script src="/static/toggle_subgraph.js"></script>' if data_source == 'scRNA' else ''}
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://code.jquery.com/ui/1.13.1/jquery-ui.min.js"></script>
    <script>
    $(function() {{
        $("#search_term").autocomplete({{
            source: function(request, response) {{
                // Only activate for term names (not IDs starting with GO:)
                if (!request.term.toUpperCase().startsWith("GO:")) {{
                    $.getJSON("/api/terms", {{
                        term: request.term
                    }}, response);
                }}
            }},
            minLength: 3,
            select: function(event, ui) {{
                $("#search_term").val(ui.item.value);
            }}
        }});
    }});
    </script>
    {f'<script src="/static/update_filters.js"></script>' if data_source == 'scRNA' else ''}
    <script>
    function runExample() {{
        const form = document.querySelector('.left-panel form');
        form.search_term.value = "GO:0048583";
        form.data_source.value = "scRNA";
        document.getElementById('cell_type_filter').value = "Microglia";
        document.getElementById('module_filter').value = "turquoise";
        form.submit();
    }}
    </script>
</body>
</html>
""")

@app.route('/exclusive_go_terms', methods=['GET', 'POST'])
def exclusive_go_terms():
    if request.method == 'POST':
        cell_type_filter = request.form.get('cell_type_filter', '').strip()
    else:
        cell_type_filter = request.args.get('cell_type_filter', '').strip()

    annotations_data = load_csv(annotations_file)
    
    unique_cell_types = set()
    if annotations_data is not None:
        cell_types = annotations_data["cell_type"].str.replace("_", " ") + " " + annotations_data["cluster"].astype(str)
        unique_cell_types.update(cell_types.unique())

    max_ic = 0
    if annotations_data is not None:
        ic_values = pd.to_numeric(annotations_data['IC'], errors='coerce')
        max_ic = ic_values[ic_values != float('inf')].max()
    replacement_ic = max_ic + 10 if not pd.isna(max_ic) else 0

    exclusive_annotations = []
    stats = {
        "exclusive_count": 0,
        "exclusive_percentage": 0.0,
        "mean_neg_log_pvalue": 0.0,
        "mean_intersection_size": 0.0,
        "mean_ic": 0.0  
    }

    if cell_type_filter and annotations_data is not None:
        cell_type_mask = (annotations_data["cell_type"].str.replace("_", " ") + " " + 
                         annotations_data["cluster"].astype(str)) == cell_type_filter
        cell_type_data = annotations_data[cell_type_mask]
        
        if not cell_type_data.empty:
            unique_terms = cell_type_data[["term_id", "term_name"]].drop_duplicates()
            total_unique_terms = len(unique_terms)
            
            exclusive_terms = []
            for _, term_row in unique_terms.iterrows():
                term_id = term_row["term_id"]
                other_occurrences = annotations_data[
                    (annotations_data["term_id"] == term_id) & ~cell_type_mask
                ]
                if other_occurrences.empty:
                    exclusive_terms.append(term_id)
            
            stats["exclusive_count"] = len(exclusive_terms)
            if total_unique_terms > 0:
                stats["exclusive_percentage"] = (stats["exclusive_count"] / total_unique_terms) * 100
            
            if exclusive_terms:
                exclusive_annotations = cell_type_data[
                    cell_type_data["term_id"].isin(exclusive_terms)
                ].to_dict('records')
                
                p_values = []
                intersection_sizes = []
                ic_values = [] 
                
                for ann in exclusive_annotations:
                    try:
                        p_value = float(ann.get("p_value", 0))
                        if p_value > 0:
                            p_values.append(-np.log10(p_value))
                    except (ValueError, TypeError):
                        pass
                    
                    try:
                        intersection_sizes.append(int(ann.get("length_intersection", 0)))
                    except (ValueError, TypeError):
                        pass
                    
                    try:
                        ic = float(ann.get("IC", 0))
                        if ic == float('inf'):
                            ic = replacement_ic
                        ic_values.append(ic)
                    except (ValueError, TypeError):
                        pass
                
                if p_values:
                    stats["mean_neg_log_pvalue"] = sum(p_values) / len(p_values)
                
                if intersection_sizes:
                    stats["mean_intersection_size"] = sum(intersection_sizes) / len(intersection_sizes)
                
                if ic_values:
                    stats["mean_ic"] = sum(ic_values) / len(ic_values)
                

    headers = ["Iteration", "Cell type", "Cluster", "Module", "Term id", "Term name", 
               "P-value", "Intersection", "Length of Intersection", "Source", 
               "Subgraph ID", "Subgraph size", "IC"]  

    exclusive_annotations.sort(key=lambda x: float(x.get("p_value", float('inf'))))

    table_rows = ""
    for result in exclusive_annotations:
        p_value_formatted = format_p_value(result.get("p_value", ""))
        intersection_formatted = format_intersection(result.get("intersection", ""))
        IC_formatted = format_ic(result.get("IC", ""))
        
        row = f"""
        <tr>
            <td class='text-center'>{result.get('iteration', '')}</td>
            <td class='text-center'>{result.get('cell_type', '').replace('_', ' ')}</td>
            <td class='text-center'>{result.get('cluster', '')}</td>
            <td class='text-center'>{result.get('module', '')}</td>
            <td class='text-center'>{result.get('term_id', '')}</td>
            <td class='text-center'>{result.get('term_name', '')}</td>
            <td class='text-center'>{p_value_formatted}</td>
            <td class='text-center'>{intersection_formatted}</td>
            <td class='text-center'>{result.get('length_intersection', '')}</td>
            <td class='text-center'>{result.get('source', '')}</td>
            <td class='text-center'>{result.get('subgraph_id', '')}</td>
            <td class='text-center'>{result.get('subgraph_size', '')}</td>
            <td class='text-center'>{IC_formatted}</td>
        </tr>
        """
        table_rows += row

    full_data_json = annotations_data.to_json(orient='records') if annotations_data is not None else "[]"

    column_descriptions = {
        "Iteration": "Number of times the pseudo-cell algorithm has been executed (e.g., T0, T5, etc.)",
        "Cell type": "Cell type where the gene is expressed",
        "Cluster": "Subgroup within the cell type",
        "Module": "Name of the module",
        "Term id": "Identifier of an annotation in the GO database",
        "Term name": "Full name of a GO annotation",
        "P-value": "Significance of an annotation given a set of genes in a co-expression module (lower = more significant)",
        "Intersection": "Name of the genes in a module that intersect with a GO database annotation",
        "Length of Intersection": "Number of genes in a module that intersect with a GO database annotation",
        "Source": "Database used to obtain the annotations (biological_process, cellular_component and molecular_function) ",
        "Subgraph ID": "GO annotation subgraph identifier",
        "Subgraph size": "Number of GO annotations that make up a subgraph",
        "IC": "Index that tells us how informative an annotation is. The higher the CI, the more specific and informative the annotation is."
    }

    return render_template_string(f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Cell Type Specific Annotations</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: grid;
            grid-template-columns: 425px 1fr 475px;
            grid-template-rows: auto auto 1fr;
            height: 100vh;
        }}
        .main-header {{
            grid-column: 1 / -1;
            background-color: white;
            color: white;
            padding: 0;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            text-align: center;
            margin-bottom: 25px;
        }}
        .main-header a {{
            background-color: #1a2a3f;
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 20px;
            display: block;
            font-size: 1.4rem;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.8);
        }}
        .main-header a:hover {{
            background-color: #6c84b4;
        }}
        .main-header a.active {{
            background-color: #6c84b4;
        }}
        .main-header a:nth-child(1) {{
            margin-right: 30px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(2) {{
            margin: 0 30px;
            border-bottom-left-radius: 8px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(3) {{
            margin-left: 30px;
            border-bottom-left-radius: 8px;
        }}
        .sub-header {{
            grid-column: 1 / -1;
            background-color: white;
            display: grid;
            grid-template-columns: 1fr minmax(200px, 500px) 100px minmax(200px, 500px) 1fr;
            text-align: center;
            padding: 0;
            align-items: stretch;
            margin-bottom: 10px;
        }}
        .sub-header a {{
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 15px 5px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.25rem;
            background-color: #354e7c;
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 1);
        }}
        .sub-header a:nth-child(1) {{ grid-column: 2; }}
        .sub-header a:nth-child(2) {{ grid-column: 4; }}
        .sub-header a:nth-child(3) {{ grid-column: 6; }}
        .sub-header a:hover {{
            background-color: #6c84b4;
        }}
        .sub-header a.active {{
            background-color: #6c84b4;
        }}
        .left-panel {{
            grid-column: 1;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .center-panel {{
            grid-column: 2;
            padding: 0px;
            overflow-y: auto;
        }}
        .right-panel {{
            grid-column: 3;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .table {{
            width: 100%;
            border-collapse: collapse;
        }}
        .table th, .table td {{
            padding: 10px;
            border: 1px solid #ddd;
            text-align: center;
        }}
        .table th {{
            background-color: #a0c4ff;
            color: black;
            font-weight: bold;
            position: sticky;
            top: 0;
        }}
        .table td {{
            white-space: normal !important;
            word-wrap: break-word;
        }}
        .results-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0;
            padding: 10px;
            background-color: white;
            position: sticky;
            top: 0;
            z-index: 1000;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        .table-container {{
            margin-top: 0;
        }}
        .metric-container {{
            margin-top: 20px;
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric-container h3 {{
            color: #1a2a4f;
            border-bottom: 2px solid #a0c4ff;
            padding-bottom: 10px;
        }}
        .metric-container p {{
            margin-bottom: 15px;
            font-size: 1.1rem;
            line-height: 1.6;
        }}
        .tooltip-icon {{
            display: inline-block;
            width: 18px;
            height: 18px;
            background-color: #1a2a4f;
            color: white;
            border-radius: 50%;
            text-align: center;
            font-size: 12px;
            line-height: 18px;
            cursor: help;
            margin-left: 5px;
        }}
        .tooltip-text {{
            visibility: hidden;
            width: 180px;
            background-color: white;
            color: #333;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 12px;
            position: fixed;
            z-index: 1000;
            box-shadow: 0 3px 12px rgba(0, 0, 0, 0.15);
            opacity: 0;
            transition: opacity 0.2s;
            font-size: 14px;
            text-align: left;
            top: auto;
            left: auto;
            transform: translateX(-50%);
            margin-top: 5px;
        }}
        th:hover .tooltip-text {{
            visibility: visible;
            opacity: 1;
        }}
        .run-example-btn {{
            background-color: #fcbf49;
            color: #ffffff;
            font-weight: bold;
            margin-left: auto;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .run-example-btn:hover {{
            background-color: #f9aa33;
            color: white;
        }}
    </style>
</head>
<body>
    <div class="main-header">
        <a href="/gene_symbol">Gene symbol</a>
        <a href="/gene_ontology_terms">Gene ontology terms</a>
        <a href="/cell_type" class="active">Cell type</a>
    </div>
    <div class="sub-header">
        <a href="/exclusive_relevant_genes">Exclusive relevant genes</a>
        <a href="/exclusive_go_terms" class="active">Exclusive GO terms</a>
    </div>
    <div class="left-panel">
        <form method="POST">
            <div class="mb-3">
                <label for="cell_type_filter" class="form-label">Select a cell type:</label>
                <select class="form-select" id="cell_type_filter" name="cell_type_filter">
                    <option value="">Select a cell type</option>
                    {"".join([f'<option value="{cell_type}" {"selected" if cell_type == cell_type_filter else ""}>{cell_type}</option>' for cell_type in sorted(unique_cell_types)])}
                </select>
            </div>
            <div class="mb-3">
                <label for="data_source" class="form-label">Database:</label>
                <select class="form-select" id="data_source" name="data_source" disabled>
                    <option value="scRNA" selected>scRNA-seq PD scCoExpNets</option>
                </select>
            </div>
            <div class="d-flex justify-content-between align-items-center">
                <button type="submit" class="btn btn-primary">Search</button>
                <button type="button" class="btn run-example-btn" onclick="runExample()">Run example</button>
            </div>
        </form>
    </div>
    <div class="center-panel">
        <div class="results-header">
            <h2>Exclusive biological functions for {cell_type_filter}</h2>
            <form method="POST" action="/download_exclusive_go_terms">
                <input type="hidden" name="cell_type_filter" value="{cell_type_filter}">
                <div class="mb-3">
                    <label for="download_format" class="form-label">Choose format:</label>
                    <select class="form-select d-inline-block w-auto" id="download_format" name="download_format">
                        <option value="csv">CSV</option>
                        <option value="xlsx">XLSX</option>
                        <option value="html">HTML</option>
                    </select>
                    <button type="submit" class="btn btn-success">Download</button>
                </div>
            </form>
        </div>
        <div class="table-container">
            <table class="table table-bordered table-hover">
                <thead class="table-primary">
                    <tr>
                        {''.join([
                            f'<th class="text-center">{header} <div class="tooltip-icon">?</div><div class="tooltip-text">{column_descriptions.get(header, "")}</div></th>'
                            for header in headers
                        ])}
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                </tbody>
            </table>
        </div>
    </div>
    <div class="right-panel">
        <div class="metric-container">
            <h3>Functional specificity profile for {cell_type_filter}</h3>
            {f"""
            <p>• {stats['exclusive_count']} annotations are exclusive to {cell_type_filter} ({stats['exclusive_percentage']:.2f}%)</p>
            <p>• Exclusive annotations show a mean –log10(p-value) of {stats['mean_neg_log_pvalue']:.2f}</p>
            <p>• Exclusive annotations show a mean intersection size of {stats['mean_intersection_size']:.2f}</p>
            <p>• Mean IC of exclusive annotations: {stats['mean_ic']:.2f}</p>
            """ if cell_type_filter else "<p>Please select a cell type to view statistics</p>"}
        </div>
    </div>
    <script>
        const fullData = {full_data_json};
    </script>
    <script src="/static/sort_table.js"></script>
    <script src="/static/toggle_subgraph.js"></script>
    <script>
    function runExample() {{
        const form = document.querySelector('.left-panel form');
        const selector = document.getElementById('cell_type_filter');
        selector.value = "microglia 19";
        form.submit();
    }}

    </script>
</body>
</html>
""")

@app.route('/new_gene_functions', methods=['GET', 'POST'])
def new_gene_functions():
    search_term = request.form.get('gene_name', '').strip().upper()
    cell_type_filter = request.form.get('cell_type_filter', '').strip()
    data_source = request.form.get('data_source', 'scRNA').strip()
    
    predicts_dir = "./DistintosPredicts/Predicts/" if data_source == 'scRNA' else "./DistintosPredicts/PredictsBulk/"
    minimal_expr_file = "./Minimally_Expressed_Statistics.csv"
    
    results = []
    new_annotations = []
    show_annotations = False  # Nuevo estado para controlar la visualización
    stats = {
        "minimal_expression": "N/A",
        "minimal_percentage": "N/A",
        "mean_mm": "N/A",
        "mean_percentile": "N/A",
        "mean_participation": "N/A",
        "mean_participation_percentage": "N/A",
        "mean_new_annotations": "N/A",
        "mean_new_percentage": "N/A",
        "mean_known_ic": "N/A",
        "mean_new_ic": "N/A",

        "known_functions": "N/A",
        "known_percentage": "N/A",
        "new_functions": "N/A",
        "new_percentage": "N/A",
        "total_annotations": "N/A",
        "mean_known_ic_bulk": "N/A",
        "mean_new_ic_bulk": "N/A"
    }
    
    if search_term:
        predict_file = f"predict_{search_term}.csv"
        predict_path = os.path.join(predicts_dir, predict_file)
        
        if os.path.exists(predict_path):
            df = pd.read_csv(predict_path)
            
            if cell_type_filter and data_source == 'scRNA':
                df = df[df['tipo_celular'].str.contains(cell_type_filter, case=False, na=False)]
            
            if data_source == 'scRNA':
                df = df.sort_values(by='new_percentage', ascending=False)
            else: 
                df = df.sort_values(by='new_functions', ascending=False)
            
            results = df.to_dict('records')
            
            if len(results) > 0:
                if data_source == 'scRNA':
                    if os.path.exists(minimal_expr_file):
                        df_minimal = pd.read_csv(minimal_expr_file)
                        gene_minimal = df_minimal[df_minimal['Gene'] == search_term]
                        if not gene_minimal.empty:
                            stats["minimal_expression"] = gene_minimal.iloc[0]['Statistic']
                            stats["minimal_percentage"] = gene_minimal.iloc[0]['Percentage']
                    
                    df_gene = df[df['module_size'].notna()]
                    
                    if len(df_gene) > 0:
                        mm_col = f"{search_term}_MM"
                        percentile_col = f"{search_term}_MM_percentile"
                        
                        if mm_col in df_gene.columns:
                            stats["mean_mm"] = f"{df_gene[mm_col].mean():.3f}"
                        if percentile_col in df_gene.columns:
                            stats["mean_percentile"] = f"{df_gene[percentile_col].mean():.1f}"
                        
                        if 'known_annotations' in df_gene.columns:
                            stats["mean_participation"] = f"{df_gene['known_annotations'].mean():.1f}"
                        
                        if 'new_annotations' in df_gene.columns:
                            stats["mean_new_annotations"] = f"{df_gene['new_annotations'].mean():.1f}"
                        if 'new_percentage' in df_gene.columns:
                            stats["mean_new_percentage"] = f"{df_gene['new_percentage'].mean():.1f}"
                        
                        if 'known_annotations' in df_gene.columns and 'new_annotations' in df_gene.columns:
                            total_annotations = df_gene['known_annotations'] + df_gene['new_annotations']
                            participation_percentage = (df_gene['known_annotations'] / total_annotations) * 100
                            stats["mean_participation_percentage"] = f"{participation_percentage.mean():.1f}"
                        
                        known_ic_values = []
                        new_ic_values = []
                        for _, row in df_gene.iterrows():
                            known_ic = row.get('IC known(CI95%)', '')
                            try:
                                if isinstance(known_ic, str) and '(' in known_ic:
                                    known_ic = known_ic.split('(')[0].strip()
                                known_ic = float(known_ic)
                                if not pd.isna(known_ic):
                                    known_ic_values.append(known_ic)
                            except (ValueError, TypeError):
                                pass
                            
                            new_ic = row.get('IC new(CI95%)', '')
                            try:
                                if isinstance(new_ic, str) and '(' in new_ic:
                                    new_ic = new_ic.split('(')[0].strip()
                                new_ic = float(new_ic)
                                if not pd.isna(new_ic):
                                    new_ic_values.append(new_ic)
                            except (ValueError, TypeError):
                                pass

                        if known_ic_values:
                            stats["mean_known_ic"] = f"{sum(known_ic_values) / len(known_ic_values):.2f}"
                        if new_ic_values:
                            stats["mean_new_ic"] = f"{sum(new_ic_values) / len(new_ic_values):.2f}"
                else:
                    if len(results) > 0:
                        first_row = results[0]
                        
                        stats["total_annotations"] = first_row.get('total_annotations', 'N/A')
                        stats["known_functions"] = first_row.get('known_functions', 'N/A')
                        stats["known_percentage"] = first_row.get('known_percentage', 'N/A')
                        stats["new_functions"] = first_row.get('new_functions', 'N/A')
                        stats["new_percentage"] = first_row.get('new_percentage', 'N/A')
                        
                        known_ic = first_row.get('IC known', '')
                        new_ic = first_row.get('IC new', '')
                        
                        try:
                            if isinstance(known_ic, str) and '(' in known_ic:
                                known_ic = known_ic.split('(')[0].strip()
                            stats["mean_known_ic_bulk"] = f"{float(known_ic):.2f}" if known_ic else "N/A"
                        except (ValueError, TypeError):
                            stats["mean_known_ic_bulk"] = "N/A"
                            
                        try:
                            if isinstance(new_ic, str) and '(' in new_ic:
                                new_ic = new_ic.split('(')[0].strip()
                            stats["mean_new_ic_bulk"] = f"{float(new_ic):.2f}" if new_ic else "N/A"
                        except (ValueError, TypeError):
                            stats["mean_new_ic_bulk"] = "N/A"

        if data_source == 'scRNA':
            annotations_data = load_csv(annotations_file)
            modules_data = []
            
            # Buscar el gen en todos los módulos
            for file_name in os.listdir(modules_dir):
                if file_name.endswith(".csv"):
                    df_module = pd.read_csv(os.path.join(modules_dir, file_name))
                    gene_data = df_module[df_module['gene'] == search_term]
                    if not gene_data.empty:
                        cell_type = extract_cell_type(file_name)
                        iteration = extract_iteration(file_name)
                        for _, row in gene_data.iterrows():
                            modules_data.append({
                                'cell_type': cell_type,
                                'iteration': iteration,
                                'cluster': row['subcluster'],
                                'module': row['module']
                            })
            
            # Buscar anotaciones correspondientes
            if annotations_data is not None and modules_data:
                for module_info in modules_data:
                    mask = (
                        (annotations_data['cell_type'] == module_info['cell_type'].replace(' ', '_')) &
                        (annotations_data['iteration'] == module_info['iteration']) &
                        (annotations_data['cluster'] == module_info['cluster']) &
                        (annotations_data['module'] == module_info['module']) &
                        (~annotations_data['intersection'].str.contains(search_term, na=False)))
                    
                    filtered = annotations_data[mask].copy()
                    if not filtered.empty:
                        # Formatear valores
                        filtered['p_value'] = filtered['p_value'].apply(format_p_value)
                        filtered['intersection'] = filtered['intersection'].apply(format_intersection)
                        filtered['IC'] = filtered['IC'].apply(format_ic)
                        filtered['cell_type'] = filtered['cell_type'].str.replace('_', ' ')
                        
                        # Aplicar filtro de cell type
                        if cell_type_filter:
                            filtered = filtered[filtered['cell_type'] == cell_type_filter]
                        
                        new_annotations.extend(filtered.to_dict('records'))

        elif data_source == 'bulk':
            annotations_df = pd.read_csv(bulk_annotations_file)
            modules_df = pd.read_csv(bulk_modules_file)

            annotations_df = annotations_df[annotations_df['cutoff'] == 10]
            modules_df = modules_df[modules_df['cutoff'] == 10]

            # Buscar módulos donde aparece el gen
            gene_modules = modules_df[modules_df['gene'] == search_term]
            new_annotations = []

            if not gene_modules.empty:
                for _, module_row in gene_modules.iterrows():
                    module = module_row['module']

                    # Filtro como en predict_got1_functions()
                    filtered_annot = annotations_df[
                        (annotations_df['target'] == 'APP') &
                        (annotations_df['tissue'] == 'DLPFC') &
                        (annotations_df['phenotype'] == 'AD') &
                        (annotations_df['module'] == module)
                    ]

                    for _, annot_row in filtered_annot.iterrows():
                        # Extraer genes de la intersección
                        genes = str(annot_row['intersection']).replace('[','').replace(']','').replace("'", '').split(',')
                        genes = [g.strip() for g in genes]
                        
                        # Si el gen buscado NO está en la intersección, es "new"
                        if search_term not in genes:
                            new_annotations.append({
                                "Cutoff": annot_row.get("cutoff"),
                                "Target": annot_row.get("target", "APP"),
                                "Tissue": annot_row.get("tissue", "DLPFC"),
                                "Phenotype": annot_row.get("phenotype", "AD"),
                                "Module": annot_row.get("module"),
                                "Term id": annot_row.get("term_id"),
                                "Term name": annot_row.get("term_name"),
                                "P-value": format_p_value(annot_row.get("p_value")),
                                "Intersection": format_intersection(annot_row.get("intersection")),
                                "Length of Intersection": annot_row.get("length_intersection", ""),
                                "Source": annot_row.get("source", ""),
                                "IC": format_ic(annot_row.get("IC", ""))
                            })


    
    # Determinar si mostrar anotaciones
    show_annotations = request.form.get('show_annotations', 'false') == 'true'

    # Configurar headers según data_source
    if data_source == 'scRNA':
        headers = [
            "Cell type", "Module size", 
            "Module membership",
            "Module membership (percentile %)", "All module terms",
            f"Term that include {search_term} ", f"Term that include {search_term} (%)", f"IC of terms that include {search_term} (95% CI)", 
            f"New terms associated with {search_term} ", f"New terms associated with {search_term} (%)", f"IC of new terms associated with {search_term} (95% CI)"
        ]
        annotation_headers = [
            "Iteration", "Cell type", "Cluster", "Module", "Term id", "Term name", 
            "P-value", "Intersection", "Length of Intersection", "Source", 
            "Subgraph ID", "Subgraph size", "IC"
        ]
        annotation_column_descriptions = {
            "Iteration": "Number of times the pseudo-cell algorithm has been executed (e.g., T0, T5, etc.)",
            "Cell type": "Cell type where the gene is expressed",
            "Cluster": "Subgroup within the cell type",
            "Module": "Name of the module",
            "Term id": "Identifier of an annotation in the GO database",
            "Term name": "Full name of a GO annotation",
            "P-value": "Significance of an annotation given a set of genes in a co-expression module (lower = more significant)",
            "Intersection": "Name of the genes in a module that intersect with a GO database annotation",
            "Length of Intersection": "Number of genes in a module that intersect with a GO database annotation",
            "Source": "Database used to obtain the annotations (biological_process, cellular_component and molecular_function)",
            "Subgraph ID": "GO annotation subgraph identifier",
            "Subgraph size": "Number of GO annotations that make up a subgraph",
            "IC": "Index that tells us how informative an annotation is. The higher the CI, the more specific and informative the annotation is."
        }

    else:
        headers = [
            "Target", "Total annotations", 
            "Known functions", "Known functions (%)", "IC known", 
            "New functions", "New functions (%)", "IC new"
        ]
        annotation_headers = [
            "Cutoff", "Target", "Tissue", "Phenotype", "Module", "Term id", "Term name", 
            "P-value", "Intersection", "Length of Intersection", "Source", "IC"
        ]
        annotation_column_descriptions = {
            "Cutoff": "Module inclusion cutoff value used for TGCN.",
            "Target": "Main variable of interest (e.g., gene or covariate).",
            "Tissue": "Brain tissue from which the data was derived.",
            "Phenotype": "Condition or group studied (e.g., AD, control).",
            "Module": "Gene module associated with the annotation.",
            "Term id": "GO term ID associated with the module.",
            "Term name": "GO term name associated with the module.",
            "P-value": "Statistical significance of the enrichment.",
            "Intersection": "Genes shared between the module and the GO term.",
            "Length of Intersection": "Number of intersecting genes.",
            "Source": "Database source of the GO term.",
            "IC": "Information Content of the annotation (higher means more specific)."
        }
    
    table_rows = ""
    for result in results:
        if data_source == 'scRNA':
            mm_col = f"{search_term}_MM" if search_term else 'Gene_MM'
            percentile_col = f"{search_term}_MM_percentile" if search_term else 'Gene_MM_percentile'
            
            row = "".join([
                f"<td class='text-center'>{result.get(col, '')}</td>"
                for col in [
                    'tipo_celular', 'module_size', mm_col, percentile_col,
                    'total_annotations', 'known_annotations', 'known_percentage',
                    'IC known(CI95%)', 'new_annotations', 'new_percentage', 'IC new(CI95%)'
                ]
            ])
        else:
            row = "".join([
                f"<td class='text-center'>{result.get(col, '')}</td>"
                for col in [
                    'target', 'total_annotations',
                    'known_functions', 'known_percentage', 'IC known',
                    'new_functions', 'new_percentage', 'IC new'
                ]
            ])
        table_rows += f"<tr>{row}</tr>"

    cell_types = set()
    if data_source == 'scRNA':
        for file_name in os.listdir(modules_dir):
            if file_name.endswith(".csv"):
                cell_type = extract_cell_type(file_name)
                if cell_type:
                    cell_types.add(cell_type)

    if data_source == 'scRNA':
        column_descriptions = {
            "Cell type": "Cell type where the gene is expressed", 
            "Module size": "Number of genes that make up the module where the gene symbol is located",
            "Module membership": "Connectivity of a gene in a module and its relevance in that module",
            "Module membership (percentile %)": "Percentile of the module membership",
            "All module terms": "Total number of annotations of the module",
            f"Term that include {search_term} ": "Number of known annotations",
            f"Term that include {search_term} (%)": "Percentage of known annotations",
            f"IC of terms that include {search_term} (95% CI)": "Index that tells us how informative the known annotations are (confidence interval 95%)", 
            f"New terms associated with {search_term} ": "Number of new annotations",
            f"New terms associated with {search_term} (%)": "Percentage of known annotations", 
            f"IC of new terms associated with {search_term} (95% CI)": "Index that tells us how informative the new annotations are (confidence interval 95%)"
        }
    else:
        column_descriptions = {
            "Target": "Main variable of interest (gene or covariate)",
            "Total annotations": "Total number of annotations for this module",
            "Known functions": "Number of known functional annotations",
            "Known functions (%)": "Percentage of known annotations from total",
            "IC known": "Information content of known annotations",
            "New functions": "Number of predicted new functional annotations",
            "New functions (%)": "Percentage of new annotations from total",
            "IC new": "Information content of new annotations"
        }

    annotation_rows = ""
    if new_annotations:
        for ann in new_annotations:
            if data_source == 'scRNA':
                # Columnas para scRNA-seq
                annotation_rows += "<tr>" + "".join([
                    f"<td class='text-center'>{ann.get(key, '')}</td>"
                    for key in [
                        "iteration", "cell_type", "cluster", "module", "term_id",
                        "term_name", "p_value", "intersection", "length_intersection",
                        "source", "subgraph_id", "subgraph_size", "IC"
                    ]
                ]) + "</tr>"
            else:
                # Columnas para bulk RNA-seq
                annotation_rows += "<tr>" + "".join([
                    f"<td class='text-center'>{ann.get(key, '')}</td>"
                    for key in [
                        "Cutoff", "Target", "Tissue", "Phenotype", "Module",
                        "Term id", "Term name", "P-value", "Intersection",
                        "Length of Intersection", "Source", "IC"
                    ]
                ]) + "</tr>"
    else:
        colspan = len(annotation_headers)
        annotation_rows = f"<tr><td colspan='{colspan}' class='text-center'>No new annotations found</td></tr>"

    return render_template_string(f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Predict Gene Functions</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: grid;
            grid-template-columns: 425px 1fr 475px;
            grid-template-rows: auto auto 1fr;
            height: 100vh;
        }}
        .main-header {{
            grid-column: 1 / -1;
            background-color: white;
            color: white;
            padding: 0;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            text-align: center;
            margin-bottom: 25px;
        }}
        .main-header a {{
            background-color: #1a2a3f;
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 20px;
            display: block;
            font-size: 1.4rem;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.8);
        }}
        .main-header a:hover {{
            background-color: #6c84b4;
        }}
        .main-header a.active {{
            background-color: #6c84b4;
        }}
        .main-header a:nth-child(1) {{
            margin-right: 30px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(2) {{
            margin: 0 30px;
            border-bottom-left-radius: 8px;
            border-bottom-right-radius: 8px;
        }}
        .main-header a:nth-child(3) {{
            margin-left: 30px;
            border-bottom-left-radius: 8px;
        }}
        .sub-header {{
            grid-column: 1 / -1;
            background-color: white;
            display: grid;
            grid-template-columns: 1fr minmax(150px, 500px) 100px minmax(150px, 500px) 100px minmax(150px, 500px) 1fr;
            text-align: center;
            padding: 0;
            align-items: stretch;
            margin-bottom: 10px;
        }}
        .sub-header a {{
            color: white;
            text-decoration: none;
            font-weight: bold;
            padding: 15px 5px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 1.25rem;
            background-color: #354e7c;
            min-width: 0;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 1);
        }}
        .sub-header a:nth-child(1) {{ grid-column: 2; }}
        .sub-header a:nth-child(2) {{ grid-column: 4; }}
        .sub-header a:nth-child(3) {{ grid-column: 6; }}
        .sub-header a:hover {{
            background-color: #6c84b4;
        }}
        .sub-header a.active {{
            background-color: #6c84b4;
        }}
        .left-panel {{
            grid-column: 1;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .center-panel {{
            grid-column: 2;
            padding: 0px;
            overflow-y: auto;
        }}
        .right-panel {{
            grid-column: 3;
            background-color: #f8f9fa;
            padding: 20px;
            overflow-y: auto;
        }}
        .table {{
            width: 100%;
            border-collapse: collapse;
        }}
        .table th, .table td {{
            padding: 10px;
            border: 1px solid #ddd;
            text-align: center;
        }}
        .table th {{
            background-color: #a0c4ff;
            color: black;
            font-weight: bold;
        }}
        .table td {{
            white-space: normal !important;
            word-wrap: break-word;
        }}
        .results-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0;
            padding: 10px;
            background-color: white;
            position: sticky;
            top: 0;
            z-index: 1000;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        .table-container {{
            margin-top: 0;
        }}
        .metric-container {{
            margin-top: 20px;
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric-container h3 {{
            color: #1a2a4f;
            border-bottom: 2px solid #a0c4ff;
            padding-bottom: 10px;
        }}
        .metric-container p {{
            margin-bottom: 15px;
            font-size: 1.1rem;
            line-height: 1.6;
        }}
        .filter-dropdown {{
            margin-top: 10px;
        }}
        .tooltip-icon {{
            display: inline-block;
            width: 18px;
            height: 18px;
            background-color: #1a2a4f;
            color: white;
            border-radius: 50%;
            text-align: center;
            font-size: 12px;
            line-height: 18px;
            cursor: help;
            margin-left: 5px;
        }}
        .tooltip-text {{
            visibility: hidden;
            width: 220px;
            background-color: white;
            color: #333;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 12px;
            position: absolute;
            z-index: 1000;
            box-shadow: 0 3px 12px rgba(0, 0, 0, 0.15);
            opacity: 0;
            transition: opacity 0.2s;
            font-size: 14px;
            text-align: left;
        }}
        th:hover .tooltip-text {{
            visibility: visible;
            opacity: 1;
        }}
        .ui-autocomplete {{
            max-height: 300px;
            overflow-y: auto;
            overflow-x: hidden;
            background-color: white;
            border: 1px solid #ddd;
            padding: 5px;
        }}
        .ui-menu-item {{
            padding: 5px;
        }}
        .ui-menu-item:hover {{
            background-color: #f0f0f0;
            cursor: pointer;
        }}
        .run-example-btn {{
            background-color: #e3a42b;  
            color: #ffffff;
            margin-left: auto;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .run-example-btn:hover {{
            background-color: #f9aa33;
            color: white;
        }}
        .annotation-toggle {{
            margin: 10px 0;
            background-color: #4CAF50;
            color: white;
            padding: 10px 20px;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            width: 230px;
        }}
        .annotation-toggle:hover {{
            background-color: #45a049;
            color: white;
        }}
        .hidden-table {{
            display: none;
        }}
    </style>
</head>
<body>
    <div class="main-header">
        <a href="/gene_symbol" class="active">Gene symbol</a>
        <a href="/gene_ontology_terms">Gene ontology terms</a>
        <a href="/cell_type">Cell type</a>
    </div>
    <div class="sub-header">
        <a href="/gene_relevance">Gene relevance</a>
        <a href="/gene_functions">Gene functions</a>
        <a href="/new_gene_functions" class="active">New gene functions</a>
    </div>
    <div class="left-panel">
        <form method="POST">
            <div class="mb-3">
                <label for="gene_name" class="form-label">Search for gene symbol:</label>
                <input type="text" class="form-control" id="gene_name" name="gene_name" 
                    value="{search_term}" placeholder="e.g., SNCA, AATF" 
                    autocomplete="off" required>
            </div>
            <div class="mb-3">
                <label for="data_source" class="form-label">Database:</label>
                <select class="form-select" id="data_source" name="data_source">
                    <option value="scRNA" {"selected" if data_source == 'scRNA' else ""}>scRNA-seq PD scCoExpNets</option>
                    <option value="bulk" {"selected" if data_source == 'bulk' else ""}>bulk RNA-seq AD TGCNs</option>
                </select>
            </div>
            <div class="d-flex justify-content-between align-items-center">
                <button type="submit" class="btn btn-primary">Search</button>
                <button type="button" class="btn run-example-btn" onclick="runExample()">Run example</button>
            </div>
            {"<input type='hidden' name='show_annotations' value='true'/>" if show_annotations else ""}
            
        </form>
    </div>
    <div class="center-panel">
        <div class="results-header">
            <h2>
                {"Predicted functions for" if data_source == "scRNA" else "Functional predictions for"} {search_term} {"based on co-expression" if data_source == "scRNA" else ""}
            </h2>
            <button type="button" class="btn annotation-toggle" onclick="toggleTables()">
                {"Show Main Table" if show_annotations else "Show Annotations Details"}
            </button>
            <form method="POST" action="/download_predict">
                <input type="hidden" name="gene_name" value="{search_term}">
                <input type="hidden" name="cell_type_filter" value="{cell_type_filter}">
                <input type="hidden" name="data_source" value="{data_source}">
                <div class="mb-3">
                    <label for="download_format" class="form-label">Choose format:</label>
                    <select class="form-select d-inline-block w-auto" id="download_format" name="download_format">
                        <option value="csv">CSV</option>
                        <option value="xlsx">XLSX</option>
                        <option value="html">HTML</option>
                    </select>
                    <button type="submit" class="btn btn-success">Download</button>
                </div>
            </form>
        </div>

        <div id="main-table" class="table-container" {'style="display:none;"' if show_annotations else ''}>
            <table class="table table-bordered table-hover">
                <thead class="table-primary">
                    <tr>
                        {''.join([
                            f'<th class="text-center">{header} <div class="tooltip-icon">?</div><div class="tooltip-text">{column_descriptions.get(header, "")}</div></th>'
                            for header in headers
                        ])}
                    </tr>
                </thead>
                <tbody>
                    {table_rows if search_term else f"<tr><td colspan='{len(headers)}' class='text-center'>Enter a gene symbol to search</td></tr>"}
                </tbody>
            </table>
        </div>
        <div id="annotation-table" class="table-container" {'style="display:none;"' if not show_annotations else ''}>
            <table class="table table-bordered table-hover">
                <thead class="table-primary">
                    <tr>
                        {''.join([
                            f'<th class="text-center">{header} <div class="tooltip-icon">?</div><div class="tooltip-text">{annotation_column_descriptions.get(header, "")}</div></th>'
                            for header in annotation_headers
                        ])}
                    </tr>
                </thead>
                <tbody>
                    {annotation_rows if new_annotations else f'''
                    <tr>
                        <td colspan="{len(annotation_headers)}" class="text-center">
                            No new annotations found
                        </td>
                    </tr>'''}
                </tbody>
            </table>
        </div>
    </div>
    <div class="right-panel">
        <div class="metric-container">
            <h3>Functional prediction metrics for {search_term}</h3>
            {f"""
            {"<p>• " + search_term + " is minimally expressed in " + stats["minimal_expression"] + " of 24 cell types (" + stats["minimal_percentage"] + ")</p>" if data_source == "scRNA" and search_term else ""}
            {"<p>• " + search_term + " has a mean module membership of " + stats["mean_mm"] + " (percentile: " + stats["mean_percentile"] + ") per cell type</p>" if data_source == "scRNA" and search_term else ""}
            {"<p>• Participates in " + stats["mean_participation"] + "% of its module's known annotations</p>" if data_source == "scRNA" and search_term else ""}
            {"<p>• " + search_term + " is involved in " + stats["mean_participation"] + " annotations on average (" + stats["mean_participation_percentage"] + "%) per cell type, with a mean IC of " + stats["mean_known_ic"] + "</p>" if data_source == "scRNA" and search_term else ""}
            {"<p>• Associated with " + stats["mean_new_annotations"] + " new annotations per cell type (" + stats["mean_new_percentage"] + "%)</p>" if data_source == "scRNA" and search_term else ""}
            {"<p>• " + search_term + " has been associated with " + stats["mean_new_annotations"] + " new annotations on average per cell type (" + stats["mean_new_percentage"] + "%), with a mean IC of " + stats["mean_new_ic"] + "</p>" if data_source == "scRNA" and search_term else ""}
            """ if search_term else ""}
            {f"""
            {"<p>• " + search_term + " is involved in " + str(stats.get("known_functions", "N/A")) + " annotations (" + str(stats.get("known_percentage", "N/A")) + "% of the total annotations of its module), which showed a mean IC of " + str(stats.get("mean_known_ic_bulk", "N/A")) + "</p>" if data_source == "bulk" and search_term else ""}
            {"<p>• " + search_term + " is associated with " + str(stats.get("new_functions", "N/A")) + " new annotations (" + str(stats.get("new_percentage", "N/A")) + "% of the total annotations of this module), which showed a mean IC of " + str(stats.get("mean_new_ic_bulk", "N/A")) + "</p>" if data_source == "bulk" and search_term else ""}
            """ if search_term else ""}
        </div>
    </div>
    <script>
        function toggleFilters() {{
            const filterDropdown = document.getElementById('filter-dropdown');
            filterDropdown.style.display = filterDropdown.style.display === 'none' ? 'block' : 'none';
        }}
    </script>
    <script src="/static/sort_table.js"></script>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://code.jquery.com/ui/1.13.1/jquery-ui.min.js"></script>
    <script>
    $(function() {{
        $("#gene_name").autocomplete({{
            source: function(request, response) {{
                $.getJSON("/api/genes", {{
                    term: request.term
                }}, response);
            }},
            minLength: 2,
            select: function(event, ui) {{
                $("#gene_name").val(ui.item.value);
                $(this).closest("form").submit();
            }},
            focus: function(event, ui) {{
                $("#gene_name").val(ui.item.value);
                return false;
            }}
        }});
    }});
    </script>
    <script>
    function runExample() {{
        const form = document.querySelector('.left-panel form');
        form.gene_name.value = "SNCA";
        form.data_source.value = "scRNA";
        form.submit();
    }}
    </script>
    <script>
    function toggleTables() {{
        const mainTable = document.getElementById('main-table');
        const annotationTable = document.getElementById('annotation-table');
        const toggleButton = document.querySelector('.annotation-toggle');

        // Alternar visibilidad
        const showAnnotations = mainTable.style.display !== 'none';

        if (showAnnotations) {{
            mainTable.style.display = 'none';
            annotationTable.style.display = 'block';
            toggleButton.innerText = 'Show Main Table';
        }} else {{
            mainTable.style.display = 'block';
            annotationTable.style.display = 'none';
            toggleButton.innerText = 'Show Annotations Details';
        }}
    }}
    </script>

</body>
</html>
""", 
    search_term=search_term,
    data_source=data_source,
    results=results,
    new_annotations=new_annotations,
    stats=stats,
    headers=headers,
    annotation_headers=annotation_headers,
    annotation_rows=annotation_rows,
    show_annotations=show_annotations)


# Distintas funciones de download para los diferentes apartados
@app.route('/download_predict', methods=['POST'])
def download_predict():
    search_term = request.form.get('gene_name', '').strip().upper()
    cell_type_filter = request.form.get('cell_type_filter', '').strip()
    download_format = request.form.get('download_format', 'csv')
    
    predicts_dir = "./DistintosPredicts/Predicts/"
    predict_file = f"predict_{search_term}.csv"
    predict_path = os.path.join(predicts_dir, predict_file)
    
    if not os.path.exists(predict_path):
        return "File not found", 404
    
    df = pd.read_csv(predict_path)
    
    if cell_type_filter:
        df = df[df['tipo_celular'].str.contains(cell_type_filter, case=False, na=False)]
    
    if download_format == 'csv':
        response = make_response(df.to_csv(index=False))
        response.headers["Content-Disposition"] = f"attachment; filename=predict_{search_term}.csv"
        response.headers["Content-Type"] = "text/csv"
        return response
    
    elif download_format == 'xlsx':
        output = BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            df.to_excel(writer, index=False, sheet_name='Results')
        response = make_response(output.getvalue())
        response.headers["Content-Disposition"] = f"attachment; filename=predict_{search_term}.xlsx"
        response.headers["Content-Type"] = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        return response
    
    elif download_format == 'html':
        response = make_response(df.to_html(index=False))
        response.headers["Content-Disposition"] = f"attachment; filename=predict_{search_term}.html"
        response.headers["Content-Type"] = "text/html"
        return response
    
@app.route('/download_gene_functions', methods=['POST'])
def download_gene_functions():
    search_term = request.form.get('gene_name', '').strip().upper()
    cell_type_filter = request.form.get('cell_type_filter', '').strip()
    iteration_filter = request.form.get('iteration_filter', '').strip()
    cluster_filter = request.form.get('cluster_filter', '').strip()
    module_filter = request.form.get('module_filter', '').strip()
    download_format = request.form.get('download_format', 'csv')
    
    annotations_data = load_csv(annotations_file)
    if annotations_data is None:
        return "Data not available", 404
    
    mask = annotations_data['intersection'].str.contains(search_term, case=False, na=False)
    filtered_data = annotations_data[mask].copy()
    
    if cell_type_filter:
        cell_filters = [f.strip().lower().replace("_", " ") for f in cell_type_filter.split(',')]
        filtered_data = filtered_data[
            filtered_data['cell_type'].str.replace("_", " ").str.lower().isin(cell_filters)
        ]
    
    if iteration_filter:
        iter_filters = [f.strip() for f in iteration_filter.split(',')]
        filtered_data = filtered_data[filtered_data['iteration'].isin(iter_filters)]
    
    if cluster_filter:
        cluster_filters = [f.strip() for f in cluster_filter.split(',')]
        filtered_data = filtered_data[filtered_data['cluster'].astype(str).isin(cluster_filters)]
    
    if module_filter:
        module_filters = [f.strip() for f in module_filter.split(',')]
        filtered_data = filtered_data[filtered_data['module'].astype(str).isin(module_filters)]
    
    if download_format == 'csv':
        response = make_response(filtered_data.to_csv(index=False))
        response.headers["Content-Disposition"] = f"attachment; filename={search_term}_functions.csv"
        response.headers["Content-Type"] = "text/csv"
        return response
    
    elif download_format == 'xlsx':
        output = BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            filtered_data.to_excel(writer, index=False, sheet_name='Results')
        response = make_response(output.getvalue())
        response.headers["Content-Disposition"] = f"attachment; filename={search_term}_functions.xlsx"
        response.headers["Content-Type"] = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        return response
    
    elif download_format == 'html':
        response = make_response(filtered_data.to_html(index=False))
        response.headers["Content-Disposition"] = f"attachment; filename={search_term}_functions.html"
        response.headers["Content-Type"] = "text/html"
        return response

@app.route('/download_go_terms', methods=['POST'])
def download_go_terms():
    data_source = request.form.get('data_source', 'scRNA').strip()
    search_term = request.form.get('search_term', '').strip()
    download_format = request.form.get('download_format', 'csv').strip()

    if data_source == 'scRNA':
        filters = {
            'cell_type_filter': request.form.get('cell_type_filter', '').strip(),
            'iteration_filter': request.form.get('iteration_filter', '').strip(),
            'cluster_filter': request.form.get('cluster_filter', '').strip(),
            'module_filter': request.form.get('module_filter', '').strip()
        }
        file_path = annotations_file 
    else:
        filters = {
            'target_filter': request.form.get('target_filter', '').strip(),
            'tissue_filter': request.form.get('tissue_filter', '').strip(),
            'cutoff_filter': request.form.get('cutoff_filter', '').strip(),
            'module_filter': request.form.get('module_filter', '').strip()
        }
        file_path = bulk_annotations_file
    
    try:
        df = pd.read_csv(file_path)
        
        if search_term:
            df = df[(df['term_id'].str.lower() == search_term.lower()) | 
                   (df['term_name'].str.lower() == search_term.lower())]
        
        if data_source == 'scRNA':
            if filters['cell_type_filter']:
                cell_filters = [f.strip().lower().replace("_", " ") for f in filters['cell_type_filter'].split(',')]
                df = df[df['cell_type'].str.replace("_", " ").str.lower().isin(cell_filters)]
            if filters['iteration_filter']:
                df = df[df['iteration'].isin(filters['iteration_filter'].split(','))]
            if filters['cluster_filter']:
                df = df[df['cluster'].astype(str).isin(filters['cluster_filter'].split(','))]
        else:
            if filters['target_filter']:
                df = df[df['target'] == filters['target_filter']]
            if filters['tissue_filter']:
                df = df[df['tissue'] == filters['tissue_filter']]
            if filters['cutoff_filter']:
                df = df[df['cutoff'].astype(str) == filters['cutoff_filter']]
        
        if filters['module_filter']:
            df = df[df['module'].astype(str).isin(filters['module_filter'].split(','))]
        
        df = df.sort_values(by='p_value')
        
        if download_format == 'csv':
            output = df.to_csv(index=False)
            response = make_response(output)
            response.headers["Content-Disposition"] = f"attachment; filename=go_term_{search_term}.csv"
            response.headers["Content-Type"] = "text/csv"
            return response
            
        elif download_format == 'xlsx':
            output = BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df.to_excel(writer, index=False, sheet_name='GO_Term_Results')
            response = make_response(output.getvalue())
            response.headers["Content-Disposition"] = f"attachment; filename=go_term_{search_term}.xlsx"
            response.headers["Content-Type"] = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            return response
            
        elif download_format == 'html':
            output = df.to_html(index=False)
            response = make_response(output)
            response.headers["Content-Disposition"] = f"attachment; filename=go_term_{search_term}.html"
            response.headers["Content-Type"] = "text/html"
            return response
            
    except Exception as e:
        return f"Error generating download: {str(e)}", 500

@app.route('/download_exclusive_go_terms', methods=['POST'])
def download_exclusive_go_terms():
    cell_type_filter = request.form.get('cell_type_filter', '').strip()
    download_format = request.form.get('download_format', 'csv').strip()
    
    try:
        annotations_data = pd.read_csv(annotations_file)
        
        cell_type_mask = (annotations_data["cell_type"].str.replace("_", " ") + " " + 
                         annotations_data["cluster"].astype(str)) == cell_type_filter
        cell_type_data = annotations_data[cell_type_mask]
        
        unique_terms = cell_type_data[["term_id", "term_name"]].drop_duplicates()
        exclusive_terms = []
        
        for _, term_row in unique_terms.iterrows():
            term_id = term_row["term_id"]
            other_occurrences = annotations_data[
                (annotations_data["term_id"] == term_id) & ~cell_type_mask
            ]
            if other_occurrences.empty:
                exclusive_terms.append(term_id)
        
        df = cell_type_data[cell_type_data["term_id"].isin(exclusive_terms)]
        df = df.sort_values(by='p_value')
        
        df = df[["iteration", "cell_type", "cluster", "module", "term_id", "term_name",
                "p_value", "intersection", "length_intersection", "source",
                "subgraph_id", "subgraph_size", "IC"]]
        
        if download_format == 'csv':
            output = df.to_csv(index=False)
            response = make_response(output)
            response.headers["Content-Disposition"] = f"attachment; filename=exclusive_go_terms_{cell_type_filter}.csv"
            response.headers["Content-Type"] = "text/csv"
            return response
            
        elif download_format == 'xlsx':
            output = BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df.to_excel(writer, index=False, sheet_name='Exclusive_GO_Terms')
            response = make_response(output.getvalue())
            response.headers["Content-Disposition"] = f"attachment; filename=exclusive_go_terms_{cell_type_filter}.xlsx"
            response.headers["Content-Type"] = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            return response
            
        elif download_format == 'html':
            output = df.to_html(index=False)
            response = make_response(output)
            response.headers["Content-Disposition"] = f"attachment; filename=exclusive_go_terms_{cell_type_filter}.html"
            response.headers["Content-Type"] = "text/html"
            return response
            
    except Exception as e:
        return f"Error generating download: {str(e)}", 500

@app.route('/download_exclusive_genes', methods=['POST'])
def download_exclusive_genes():
    cell_type_filter = request.form.get('cell_type_filter', '').strip()
    download_format = request.form.get('download_format', 'csv').strip()
    
    if not cell_type_filter:
        return "Please select a cell type first", 400
    
    try:
        df_min = pd.read_csv("./Minimally_Expressed_Statistics.csv")
        df_t0 = pd.read_csv("./Relevant_At_T0_Statistics.csv")
        df_all = pd.read_csv("./Relevant_In_All_Iterations_Statistics.csv")
        
        normalized_filter = cell_type_filter.replace("_", " ")
        
        def get_exclusive_genes(df):
            genes_in_cell = set(
                df[df["Cell Types"]
                .fillna("")
                .str.replace("_", " ")
                .str.contains(normalized_filter)]["Gene"]
            )

            exclusive_genes = []
            for gene in genes_in_cell:
                rows = df[df["Gene"] == gene]
                cell_types_series = rows["Cell Types"].dropna().astype(str)
                all_cell_types = (
                    cell_types_series
                    .str.replace("_", " ")
                    .str.split("; ")
                    .explode()
                    .dropna()
                    .str.strip()
                )

                if all(ct == normalized_filter for ct in all_cell_types):
                    exclusive_genes.append(gene)

            return exclusive_genes
        
        minimal_genes = get_exclusive_genes(df_min)
        t0_genes = get_exclusive_genes(df_t0)
        all_iter_genes = get_exclusive_genes(df_all)
        
        data = []
        for gene in minimal_genes:
            data.append({"Gene": gene, "Criteria": "Minimally expressed", "Cell Type": cell_type_filter})
        for gene in t0_genes:
            data.append({"Gene": gene, "Criteria": "Relevant at T0", "Cell Type": cell_type_filter})
        for gene in all_iter_genes:
            data.append({"Gene": gene, "Criteria": "Relevant in all iterations", "Cell Type": cell_type_filter})
        
        df = pd.DataFrame(data)
        
        if download_format == 'csv':
            output = df.to_csv(index=False)
            response = make_response(output)
            response.headers["Content-Disposition"] = f"attachment; filename=exclusive_genes_{cell_type_filter}.csv"
            response.headers["Content-Type"] = "text/csv"
            return response
            
        elif download_format == 'xlsx':
            output = BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df.to_excel(writer, index=False, sheet_name='Exclusive_Genes')
            response = make_response(output.getvalue())
            response.headers["Content-Disposition"] = f"attachment; filename=exclusive_genes_{cell_type_filter}.xlsx"
            response.headers["Content-Type"] = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            return response
            
        elif download_format == 'html':
            output = df.to_html(index=False)
            response = make_response(output)
            response.headers["Content-Disposition"] = f"attachment; filename=exclusive_genes_{cell_type_filter}.html"
            response.headers["Content-Type"] = "text/html"
            return response
            
    except Exception as e:
        return f"Error generating download: {str(e)}", 500

@app.route('/download', methods=['POST'])
def download():
    search_term = request.form.get('gene_name', '').strip()
    cell_type_filter = request.form.get('cell_type_filter', '').strip()
    iteration_filter = request.form.get('iteration_filter', '').strip()
    cluster_filter = request.form.get('cluster_filter', '').strip()
    module_filter = request.form.get('module_filter', '').strip()
    percentile_filter = request.form.get('percentile_filter', '').strip()
    download_format = request.form.get('download_format', 'csv')

    percentile_filter = float(percentile_filter) if percentile_filter else None
    results = query_dataset("modules", search_term, cell_type_filter, iteration_filter, 
                          cluster_filter, module_filter, percentile_filter)

    df = pd.DataFrame(results)

    if download_format == 'csv':
        response = make_response(df.to_csv(index=False))
        response.headers["Content-Disposition"] = "attachment; filename=results.csv"
        response.headers["Content-Type"] = "text/csv"
        return response

    elif download_format == 'xlsx':
        output = BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            df.to_excel(writer, index=False, sheet_name='Results')
        response = make_response(output.getvalue())
        response.headers["Content-Disposition"] = "attachment; filename=results.xlsx"
        response.headers["Content-Type"] = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        return response

    elif download_format == 'html':
        response = make_response(df.to_html(index=False))
        response.headers["Content-Disposition"] = "attachment; filename=results.html"
        response.headers["Content-Type"] = "text/html"
        return response

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')