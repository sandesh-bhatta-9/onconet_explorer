import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from collections import Counter
from bioservices import KEGG
import concurrent.futures
import base64
from pathlib import Path

# --- 1. Page Configuration ---
st.set_page_config(page_title="Cancer Pathway Analyzer", layout="wide")
st.title("🔬 Cancer-Specific Pathway Analyzer")
st.write("Compare cancer pathways against chemotherapy and natural product pathways.")

# --- 2. Icon & KEGG Utilities ---

ICON_NAMES = ["Gene Cards", "NCBI", "ENSEMBL", "KEGG", "GEO"]
ICONS_DIR = Path(__file__).parent / "Icons"

@st.cache_data
def load_icons_b64():
    """Loads icons from the 'Icons' folder and encodes them."""
    icon_map = {}
    for name in ICON_NAMES:
        fp = ICONS_DIR / f"{name}.png"
        if fp.is_file():
            with open(fp, "rb") as f:
                b64 = base64.b64encode(f.read()).decode("utf-8")
                icon_map[name] = f"data:image/png;base64,{b64}"
        else:
            icon_map[name] = None
    return icon_map

ICON_B64 = load_icons_b64()

@st.cache_data
def get_all_kegg_pathways():
    """Fetches all human (hsa) pathways from KEGG."""
    k = KEGG()
    pathways_raw = k.list("pathway/hsa")
    d = {}
    for line in pathways_raw.strip().split('\n'):
        pid, name_desc = line.split('\t')
        name = name_desc.split(' - ')[0]
        d[name] = pid.replace('path:', '')
    return d

@st.cache_data
def get_genes_from_pathway(pathway_id: str):
    """
    Fetches and parses genes for a single KEGG pathway using the robust bioservices parser.
    """
    k = KEGG()
    k.organism = "hsa"
    genes = set()
    try:
        data = k.get(pathway_id)
        if not data:
            return set()
        
        parsed = k.parse(data) # Use the built-in parser
        if 'GENE' not in parsed:
            return set()

        gene_data = parsed['GENE']
        
        if isinstance(gene_data, dict):
            for description in gene_data.values():
                symbol = description.split(';')[0]
                genes.add(symbol)
        elif isinstance(gene_data, list):
            for entry in gene_data:
                parts = entry.split()
                if len(parts) > 1:
                    symbol = parts[1].replace(';', '')
                    genes.add(symbol)
        return genes
    except Exception as e:
        st.error(f"Error parsing pathway {pathway_id}: {e}")
        return set()

def generate_url_links(gene_name: str) -> dict:
    """Generates plain URLs for the CSV download."""
    return {
        "GeneCards_URL": f"https://www.genecards.org/Search/Keyword?queryString={gene_name}",
        "NCBI_URL": f"https://www.ncbi.nlm.nih.gov/gene/?term={gene_name}",
        "ENSEMBL_URL": f"https://useast.ensembl.org/Search/Results?q={gene_name}",
        "GEO_URL": f"https://www.ncbi.nlm.nih.gov/gds/?term={gene_name}",
    }

def generate_icon_links(gene_name: str) -> dict:
    """
    Generates HTML links using the pre-loaded B64 icons.
    """
    urls = {
        "Gene Cards": f"https://www.genecards.org/Search/Keyword?queryString={gene_name}",
        "NCBI": f"https://www.ncbi.nlm.nih.gov/gene/?term={gene_name}",
        "ENSEMBL": f"https://useast.ensembl.org/Search/Results?q={gene_name}",
        "GEO": f"https://www.ncbi.nlm.nih.gov/gds/?term={gene_name}",
    }
    html = {}
    for db, url in urls.items():
        b64src = ICON_B64.get(db)
        if b64src:
            html[db] = (f'<a href="{url}" target="_blank"><img src="{b64src}" width="24" height="24" style="margin:2px; border-radius:4px; border:1px solid #ccc;"></a>')
        else:
            html[db] = f'<a href="{url}" target="_blank">{db}</a>'
    return html

# --- 3. Sidebar ---
all_paths = get_all_kegg_pathways()
st.sidebar.header("Select Pathways for Comparison")

CANCER_KEYWORDS = [
    'cancer', 'glioma', 'leukemia', 'lymphoma', 'melanoma', 
    'gastric', 'colorectal', 'prostate', 'breast'
]
CHEMO_KEYWORDS = ['fluoropyrimidine', 'folate', 'platinum']

NATURAL_PRODUCT_KEYWORDS = [
    'metabolism of xenobiotics by cytochrome p450',
    'drug metabolism - cytochrome p450',
    'drug metabolism - other enzymes',
    'steroid hormone biosynthesis',
    'retinol metabolism',
    'metabolism',
]

cancer_options = {n:pid for n,pid in all_paths.items() if any(k in n.lower() for k in CANCER_KEYWORDS)}
chemo_options = {n:pid for n,pid in all_paths.items() if any(k in n.lower() for k in CHEMO_KEYWORDS)}
natural_options = {n:pid for n,pid in all_paths.items() if any(k in n.lower() for k in NATURAL_PRODUCT_KEYWORDS)}


# --- *** NEW: EXAMPLE BUTTONS & SESSION STATE *** ---
st.sidebar.markdown("---")
st.sidebar.subheader("Load Examples")

def set_example_state(simple=True):
    """Callback function to set the session state for example buttons."""
    if simple:
        # Simple Example: 1 Cancer vs 1 Chemo
        st.session_state.cancer_key = [k for k in cancer_options.keys() if "Melanoma" in k]
        st.session_state.chemo_key = [k for k in chemo_options.keys() if "Platinum" in k]
        st.session_state.natural_key = []
    else:
        # Complex Example: 3 Cancers vs 1 Chemo vs 1 Natural
        st.session_state.cancer_key = [k for k in cancer_options.keys() if k in ["Melanoma", "Renal cell carcinoma", "Gastric cancer"]]
        st.session_state.chemo_key = [k for k in chemo_options.keys() if "Platinum" in k]
        st.session_state.natural_key = [k for k in natural_options.keys() if "cytochrome P450" in k]

# Initialize session state keys for multiselect widgets
# This sets the *initial* default when the app first loads
if "cancer_key" not in st.session_state:
    st.session_state.cancer_key = [k for k in cancer_options.keys() if k in ["Melanoma", "Renal cell carcinoma"]]
if "chemo_key" not in st.session_state:
    st.session_state.chemo_key = []
if "natural_key" not in st.session_state:
    st.session_state.natural_key = []

st.sidebar.button("Load Simple Example (1 vs 1)", on_click=set_example_state, args=(True,))
st.sidebar.button("Load Complex Example (3 vs 2)", on_click=set_example_state, args=(False,))
st.sidebar.markdown("---")
# --- *** END OF NEW SECTION *** ---


st.sidebar.subheader("Cancer Pathways")
sel_cancer = st.sidebar.multiselect(
    "Select cancer types:", 
    list(cancer_options.keys()), 
    key="cancer_key" # <-- MODIFIED: Use session state key
)

st.sidebar.subheader("Therapeutic Pathways")
sel_chemo = st.sidebar.multiselect(
    "Select chemotherapy pathways:", 
    list(chemo_options.keys()), 
    key="chemo_key" # <-- MODIFIED: Use session state key
)
sel_natural = st.sidebar.multiselect(
    "Select natural product pathways:", 
    list(natural_options.keys()), 
    key="natural_key" # <-- MODIFIED: Use session state key
)

combined_options = {**cancer_options, **chemo_options, **natural_options}
selected_pathways = sel_cancer + sel_chemo + sel_natural

# --- 4. Main App Logic ---
if selected_pathways:
    
    # --- 4a. Parallel Data Fetching ---
    p2g = {} # Pathway-to-Genes map
    with st.spinner(f"Fetching {len(selected_pathways)} pathways in parallel..."):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_to_pathway = {
                executor.submit(get_genes_from_pathway, combined_options[p_name]): p_name 
                for p_name in selected_pathways
            }
            for future in concurrent.futures.as_completed(future_to_pathway):
                p_name = future_to_pathway[future]
                p2g[p_name] = future.result()

    # --- 4b. Gene Mapping for New Coloring Logic ---
    flat_genes = [g for genes in p2g.values() for g in genes]
    gene_counts = Counter(flat_genes)
    
    g2c = {} # Gene-to-Cancers map
    for gene in gene_counts:
        g2c[gene] = set()
        for cancer_name in sel_cancer:
            if gene in p2g[cancer_name]:
                g2c[gene].add(cancer_name)

    # --- 4c. Network Visualization (MOVED FIRST) ---
    st.header("🕸️ Pathway Network")
    st.info(
        "**Note on Layout:** The distance between nodes (e.g., 'Melanoma' and 'Renal cell carcinoma') "
        "is determined by the layout algorithm based on shared connections (genes). "
        "It is a visualization artifact, not a direct measure of biological similarity."
    )
    
    G = nx.Graph()
    for pathway in p2g: 
        G.add_node(pathway, type="pathway")
    
    for gene, cnt in gene_counts.items(): 
        G.add_node(gene, type="gene", count=cnt)
    
    for pathway, genes in p2g.items():
        for gene in genes: 
            G.add_edge(pathway, gene)
    
    pos = nx.spring_layout(G, k=0.5, seed=42)
    
    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]; x1, y1 = pos[v]
        edge_x += [x0, x1, None]; edge_y += [y0, y1, None]
    edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(color='#888', width=0.5), hoverinfo='none')

    node_x, node_y, text, color = [], [], [], []
    total_selected_cancers = len(sel_cancer)

    for node, attr in G.nodes(data=True):
        x, y = pos[node]
        node_x.append(x); node_y.append(y)
        
        if attr["type"] == "pathway":
            if node in sel_cancer: color.append("crimson")
            elif node in sel_chemo: color.append("mediumblue")
            elif node in sel_natural: color.append("forestgreen")
            else: color.append("grey")
            text.append(node)
        
        else:
            cancers_for_gene = g2c.get(node, set())
            num_cancers = len(cancers_for_gene)
            
            if num_cancers == total_selected_cancers and total_selected_cancers > 1:
                color.append("hotpink")
                text.append(f"{node} (Shared by ALL {num_cancers} cancers)")
            elif 1 < num_cancers < total_selected_cancers:
                color.append("cyan")
                text.append(f"{node} (Shared by {num_cancers} cancers)")
            elif num_cancers == 1:
                color.append("orange")
                text.append(f"{node} (Unique to {list(cancers_for_gene)[0]})")
            else:
                color.append("lightgray")
                text.append(f"{node} (Not in selected cancers)")
            
    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers', marker=dict(size=10, color=color, line_width=2), text=text, hoverinfo='text')
    fig_title = (
        "Cancer, Chemotherapy, and Natural Product Pathway Network<br>"
        "<sup>Colors: crimson (Cancer), mediumblue (Chemo), forestgreen (Natural) | "
        "hotpink (Gene: All Cancers), cyan (Gene: Some Cancers), orange (Gene: One Cancer)</sup>"
    )
    fig = go.Figure(
        data=[edge_trace, node_trace], 
        layout=go.Layout(
            title=fig_title, 
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False), 
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False), 
            margin=dict(b=20, l=5, r=5, t=60), 
            hovermode='closest'
        )
    )
    st.plotly_chart(fig, use_container_width=True)
    st.markdown("---")


    # --- 4d. Shared Gene Analysis (MOVED SECOND) ---
    st.header("📊 Shared Gene Analysis")
    
    # --- *** NEW: GENE SEARCH BOX *** ---
    gene_filter = st.text_input("Search for a specific gene in the tables below (e.g., 'TP53')").strip().upper()
    # --- *** END OF NEW SECTION *** ---

    genes_all_cancers = sorted([g for g, c in g2c.items() if len(c) == total_selected_cancers])
    genes_some_cancers = sorted([g for g, c in g2c.items() if 1 < len(c) < total_selected_cancers])
    genes_one_cancer = sorted([g for g, c in g2c.items() if len(c) == 1])

    analysis_tabs = st.tabs(["🧬 Shared by ALL Cancers", "🧬 Shared by SOME Cancers", "🧬 Unique to ONE Cancer"])

    with analysis_tabs[0]:
        st.subheader(f"Genes in {total_selected_cancers} of {total_selected_cancers} selected cancers")
        if not genes_all_cancers and total_selected_cancers > 1:
            st.warning("No genes were found to be common to ALL selected cancer pathways.")
        elif total_selected_cancers <= 1:
            st.info("Select at least two cancer pathways to see this comparison.")
        else:
            st.write(f"Found {len(genes_all_cancers)} genes common to all selected cancers.")
            display_data = [{'Gene': g, 'Pathways': ', '.join([n for n, gs in p2g.items() if g in gs]), **generate_icon_links(g)} for g in genes_all_cancers]
            df_display = pd.DataFrame(display_data).sort_values(by="Gene")
            
            # --- MODIFIED: Apply Filter ---
            if gene_filter:
                df_display = df_display[df_display['Gene'].str.upper().str.contains(gene_filter)]
            
            cols_display = ['Gene','Pathways','Gene Cards','NCBI','ENSEMBL','GEO']
            st.write(df_display[cols_display].to_html(escape=False,index=False), unsafe_allow_html=True)
            
    with analysis_tabs[1]:
        st.subheader(f"Genes in 2 to {total_selected_cancers - 1} selected cancers")
        if not genes_some_cancers:
            st.info("No genes were found to be shared by *some* (but not all) selected cancers.")
        else:
            st.write(f"Found {len(genes_some_cancers)} genes shared by some selected cancers.")
            display_data = [{'Gene': g, 'Cancers': ', '.join(g2c[g]), **generate_icon_links(g)} for g in genes_some_cancers]
            df_display = pd.DataFrame(display_data).sort_values(by="Gene")

            # --- MODIFIED: Apply Filter ---
            if gene_filter:
                df_display = df_display[df_display['Gene'].str.upper().str.contains(gene_filter)]

            cols_display = ['Gene','Cancers','Gene Cards','NCBI','ENSEMBL','GEO']
            st.write(df_display[cols_display].to_html(escape=False,index=False), unsafe_allow_html=True)

    with analysis_tabs[2]:
        st.subheader("Genes unique to a single selected cancer pathway")
        if not genes_one_cancer:
            st.info("No genes were found to be unique to a single selected cancer pathway.")
        else:
            st.write(f"Found {len(genes_one_cancer)} genes unique to one selected cancer.")
            display_data = [{'Gene': g, 'Cancer Pathway': list(g2c[g])[0], **generate_icon_links(g)} for g in genes_one_cancer]
            df_display = pd.DataFrame(display_data).sort_values(by="Cancer Pathway")
            
            # --- MODIFIED: Apply Filter ---
            if gene_filter:
                df_display = df_display[df_display['Gene'].str.upper().str.contains(gene_filter)]
            
            cols_display = ['Gene','Cancer Pathway','Gene Cards','NCBI','ENSEMBL','GEO']
            st.write(df_display[cols_display].to_html(escape=False,index=False), unsafe_allow_html=True)


    # --- 4e. Download Section ---
    st.markdown("---")
    st.subheader("📥 Download Data")
    if genes_all_cancers:
        download_data = []
        for g in genes_all_cancers:
            row = {
                'Gene Name': g,
                'Pathways': ', '.join([n for n, gs in p2g.items() if g in gs]),
                'Cancers': ', '.join(g2c[g]),
                'Frequency (All Pathways)': gene_counts[g],
                **generate_url_links(g) # Download CSV uses plain text URLs
            }
            download_data.append(row)
        df_download = pd.DataFrame(download_data)
        
        column_order = ['Gene Name', 'Pathways', 'Cancers', 'Frequency (All Pathways)', 'GeneCards_URL', 'NCBI_URL', 'ENSEMBL_URL', 'GEO_URL']
        df_download = df_download[column_order].sort_values(by="Gene Name")

        csv_string = df_download.to_csv(index=False).encode('utf-8')
        st.download_button(
            label=f"Download 'Shared by ALL' List ({len(df_download)} genes)", 
            data=csv_string, 
            file_name="shared_genes_ALL_cancers.csv", 
            mime="text/csv",
            key="dl_all_cancers"
        )
    else:
        st.write("No 'Shared by ALL' data to download.")

else:
    st.info("☝️ Please select at least one pathway from the sidebar to begin.")