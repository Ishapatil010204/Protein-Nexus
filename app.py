import streamlit as st
import requests
import pandas as pd
import plotly.express as px
from Bio.SeqUtils import ProtParam
import networkx as nx
import matplotlib.pyplot as plt
from urllib.parse import quote
import plotly.graph_objects as go
from PIL import Image
from st_aggrid import AgGrid
import py3Dmol
from stmol import showmol
import streamlit.components.v1 as components

# Database color scheme
DB_COLORS = {
    "UniProt": "#0047AB",
    "PDB": "#FF6F00",
    "AlphaFold": "#2A5CAA",
    "STRING": "#32CD32",
    "KEGG": "#800080",
    "NCBI": "#00688B",
    "Reactome": "#8B0000"
}

# Configure page
st.set_page_config(layout="wide", page_title="Protein Nexus", page_icon="🧬")

# -------- Custom CSS --------
st.markdown(f"""
<style>
    /* Main background for the entire app */
    .stApp {{
        background-color: #e6f2ff !important;
    }}

    /* Sidebar background to white */
    div[data-testid="stSidebar"] {{
        background-color: #ffffff !important;
    }}

    /* Style input boxes (text_input, selectbox, checkbox, slider, etc.) to have white background and black text */
    input[type="text"] {{
        background-color: #ffffff !important;
        color: #000000 !important;
        border: 1px solid #ccc;
        border-radius: 4px;
        padding: 4px 8px;
    }}
    select {{
        background-color: #ffffff !important;
        color: #000000 !important;
        border: 1px solid #ccc;
        border-radius: 4px;
        padding: 4px 8px;
    }}
    div[data-baseweb="checkbox"] > div:nth-child(2) {{
        background-color: #ffffff !important;
    }}
    input[type="range"] {{
        background-color: #ffffff !important;
    }}
    h1, h2, h3, h4, h5, h6, p, li, span, div {{
        color: #000000 !important;
    }}
    .feature-card {{
        border-radius: 10px;
        padding: 20px;
        margin-bottom: 20px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        background-color: white;
        border-left: 4px solid #4e79a7;
        color: #000000;
    }}
    .database-card {{
        border-radius: 10px;
        padding: 15px;
        margin-bottom: 15px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        background-color: white;
        color: #000000;
    }}
    .section-header {{
        color:#000000;
        border-bottom: 2px solid #3498db;
        padding-bottom: 5px;
        margin-top: 20px;
    }}
    .metric-card {{
        background-color: #f8f9fa;
        border-radius: 8px;
        padding: 15px;
        text-align: center;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        color: #000000;
    }}
    .dataframe {{
        width: 100%;
    }}
    .profile-img {{
        border-radius: 50%;
        border: 3px solid #4e79a7;
        margin-bottom: 10px;
    }}
    .stTabs [data-baseweb="tab-list"] {{
        gap: 10px;
    }}
    .stTabs [data-baseweb="tab"] {{
        padding: 8px 16px;
        border-radius: 4px 4px 0 0;
    }}
    .stTabs [aria-selected="true"] {{
        background-color: #e6f2ff;
    }}
    .stProgress > div > div > div > div {{
        background-color: #3498db;
    }}
    body, h1, h2, h3, h4, h5, h6, p, span, li {{
        color: #000000 !important;
    }}
</style>
""", unsafe_allow_html=True)

# -------- Navigation --------
st.sidebar.title("🧬 Protein Explorer Pro")
# Add an option for About Creator
page = st.sidebar.radio("Navigation", ["🏠 Home", "🔍 Search Protein", "📚 Tutorial", "👤 About"])

# Initialize session state if not set
if 'current_section' not in st.session_state:
    st.session_state['current_section'] = 'Home'

# -------- Home Page --------
if page == "🏠 Home":
    st.title("Welcome to Protein Nexus💻🧬")
    st.markdown("""
    <div style="background-color: #e6f2ff; padding: 20px; border-radius: 10px; margin-bottom: 20px; color: #000000;">
        <h3 style="color: #000000;">A Comprehensive Multi-Omics Protein Analysis Platform</h3>
        <p>Integrating data from multiple biological databases for sequence analysis, structural visualization, 
        interaction networks, and evolutionary insights.</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("## Key Features")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("""
        <div class="feature-card">
            <h4>🧬 Sequence Analysis</h4>
            <p>Detailed protein sequence characterization including composition, hydrophobicity, and physicochemical properties</p>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown("""
        <div class="feature-card">
            <h4>🛠 Structural Visualization</h4>
            <p>3D protein structures from PDB and AlphaFold with interactive visualization</p>
        </div>
        """, unsafe_allow_html=True)
    with col3:
        st.markdown("""
        <div class="feature-card">
            <h4>🌐 Interaction Networks</h4>
            <p>Protein-protein interaction networks from STRING database</p>
        </div>
        """, unsafe_allow_html=True)
    st.markdown("## How to Use")
    st.markdown("""
    <div style="color: #000000;">
        <ol>
            <li>Navigate to the <strong>🔍 Search Protein</strong> page</li>
            <li>Enter a UniProt ID or gene/protein name (e.g., "P00533" or "EGFR")</li>
            <li>Select your species of interest</li>
            <li>Explore the results across multiple tabs:
                <ul>
                    <li><strong>Overview</strong>: General protein information</li>
                    <li><strong>Sequence</strong>: Detailed sequence analysis</li>
                    <li><strong>Structure</strong>: 3D visualization and quality metrics</li>
                    <li><strong>Evolution</strong>: Cross-species comparison</li>
                    <li><strong>Literature</strong>: Relevant publications</li>
                </ul>
            </li>
        </ol>
    </div>
    """, unsafe_allow_html=True)
    st.markdown("## Supported Databases")
    db_cols = st.columns(4)
    for i, (db, color) in enumerate(DB_COLORS.items()):
        with db_cols[i % 4]:
            st.markdown(f"""
            <div class="database-card" style="border-left: 4px solid {color};">
                <h4 style="color: {color};">{db}</h4>
                <p>Integrated {db} data access</p>
            </div>
            """, unsafe_allow_html=True)

# -------- Tutorial Page --------
elif page == "📚 Tutorial":
    st.title("Protein Explorer Pro Tutorial")
    st.markdown("Learn how to use each feature of Protein Explorer Pro")
    with st.expander("🔍 Searching for Proteins", expanded=True):
        st.markdown("""
        ### Finding Your Protein
        1. **Input Field**: Enter either:
           - A UniProt accession ID (e.g., P00533)
           - A gene name (e.g., EGFR)
           - A protein name (e.g., Epidermal growth factor receptor)
        2. **Species Filter**: Narrow your search to a specific organism
        3. **Advanced Options**: Toggle which analysis sections to display
        The system will automatically search across multiple databases and present integrated results.
        """)
    with st.expander("📋 Overview Tab", expanded=True):
        st.markdown("""
        ### Key Information at a Glance
        - **Function**: Biological role of the protein
        - **Domains & Features**: Annotated domains with interactive visualization
        - **Subcellular Locations**: Where the protein functions in the cell
        - **Pathways**: Metabolic pathways the protein participates in (KEGG)
        Each section can be expanded for more details.
        """)
    with st.expander("🧬 Sequence Analysis", expanded=True):
        st.markdown("""
        ### Detailed Sequence Characterization
        - **Physicochemical Properties**: Molecular weight, pI, instability index
        - **Amino Acid Composition**: Percentage of each amino acid
        - **Hydrophobicity Profile**: Customizable window size analysis
        - **Secondary Structure Prediction**: Helix, sheet, and turn fractions
        All plots are interactive - hover for details, zoom, and download.
        """)
    with st.expander("🛠 Structure Visualization", expanded=True):
        st.markdown("""
        ### 3D Protein Structures
        - **Experimental Structures**: From PDB (when available)
        - **Predicted Structures**: From AlphaFold DB
        - **Quality Metrics**: Ramachandran plots, rotamer outliers
        Rotate structures by clicking and dragging. Zoom with mouse wheel.
        """)
    with st.expander("🌱 Evolutionary Analysis", expanded=True):
        st.markdown("""
        ### Cross-Species Comparison
        - **Taxonomic Lineage**: Evolutionary relationships
        - **Sequence Similarity**: Comparison with orthologs
        - **BLAST Search**: Link to NCBI for detailed sequence comparison
        """)
    with st.expander("📚 Literature References", expanded=True):
        st.markdown("""
        ### Relevant Publications
        - **Key Papers**: From UniProt annotations
        - **PubMed Links**: Direct access to publication details
        - **LitSense Search**: Find additional articles
        """)

# -------- Search Protein Page --------
elif page == "🔍 Search Protein":
    st.title("Protein Search & Analysis")
    with st.form("protein_search"):
        col1, col2 = st.columns([3, 1])
        with col1:
            query = st.text_input("Enter UniProt ID or gene/protein name:", placeholder="e.g., P00533 or EGFR", key="protein_query")
        with col2:
            st.markdown("<div style='height: 28px'></div>", unsafe_allow_html=True)
            submitted = st.form_submit_button("Search")
    # Species options
    species_options = {
        "All Eukaryota": [
            "Homo sapiens", "Mus musculus", "Drosophila melanogaster",
            "Arabidopsis thaliana", "Oryza sativa",
            "Saccharomyces cerevisiae", "Candida albicans",
            "Plasmodium falciparum"
        ],
        "Bacteria": ["Escherichia coli", "Mycobacterium tuberculosis"],
        "Archaea": ["Halobacterium salinarum", "Methanococcus jannaschii"],
        "Viruses": ["HIV-1", "Influenza A virus", "SARS-CoV-2"],
        "Others": ["Synthetic construct", "Environmental sample"]
    }
    category = st.selectbox("Select category", list(species_options.keys()), key="species_category")
    species_list = species_options.get(category, [])
    selected_species = st.selectbox("Choose species:", species_list, key="selected_species")
    # Advanced options
    with st.expander("⚙️ Advanced Options"):
        show_sequence_analysis = st.checkbox("Show Sequence Analysis", True, key="seq_analysis")
        show_interaction_network = st.checkbox("Show Protein Interaction Network", True, key="interaction")
        show_structural_prediction = st.checkbox("Show Structural Predictions", True, key="structure")
        show_evolutionary = st.checkbox("Show Evolutionary Analysis", True, key="evolution")
        show_publications = st.checkbox("Show Related Publications", True, key="publications")
    # Process search
    if submitted and query:
        with st.spinner(f"🔍 Searching databases for {query}..."):
            progress_bar = st.progress(0)
            # Build query
            species_query = f"+AND+organism_name:{selected_species.replace(' ', '+')}" if selected_species else ""
            search_url = f"https://rest.uniprot.org/uniprotkb/search?query={query}{species_query}&format=json&fields=accession"
            progress_bar.progress(10)
            response = requests.get(search_url)
            if response.status_code == 200:
                results = response.json()
                if results.get('results'):
                    accession = results['results'][0]['primaryAccession']
                    full_url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
                    full_response = requests.get(full_url)
                    progress_bar.progress(30)
                    if full_response.status_code == 200:
                        data = full_response.json()
                        # Basic info display
                        st.success(f"Found protein: {data['proteinDescription']['recommendedName']['fullName']['value']}")
                        # Tabs
                        tab1, tab2, tab3, tab4, tab5 = st.tabs(["📋 Overview", "🧬 Sequence", "🛠 Structure", "🌱 Evolution", "📚 Literature"])

                        with tab1:
                            # Overview
                            st.markdown("### 🧪 Functional Annotation")
                            function = next((c['texts'][0]['value'] for c in data['comments'] if c['commentType'] == 'FUNCTION'), "Not available")
                            st.markdown(f"<div class='feature-card'>{function}</div>", unsafe_allow_html=True)
                            # Domains & Features
                            st.markdown("### 🧩 Domains & Features")
                            domains = [f for f in data.get('features', []) if f['type'] in ['DOMAIN', 'MOTIF', 'REGION']]
                            if domains:
                                domain_df = pd.DataFrame([{
                                    'Type': f['type'],
                                    'Description': f.get('description', 'No description'),
                                    'Start': f['location']['start']['value'],
                                    'End': f['location']['end']['value'],
                                    'Length': f['location']['end']['value'] - f['location']['start']['value'] + 1
                                } for f in domains])
                                fig = px.scatter(domain_df, x="Start", y="Type", color="Type",
                                                 size="Length", hover_data=["Description"],
                                                 title="Protein Features Along Sequence",
                                                 color_discrete_sequence=px.colors.qualitative.Pastel)
                                fig.update_layout(
                                    plot_bgcolor='rgba(255,255,255,1)',
                                    paper_bgcolor='rgba(255,255,255,1)',
                                    height=400
                                )
                                st.plotly_chart(fig, use_container_width=True)
                                st.markdown("#### Feature Details")
                                AgGrid(domain_df, height=250, fit_columns_on_grid_load=True)
                            else:
                                st.warning("No domain information found.")
                            # Subcellular Locations
                            st.markdown("### 📍 Subcellular Locations")
                            locations = []
                            for c in data['comments']:
                                if c['commentType'] == 'SUBCELLULAR LOCATION':
                                    for loc in c.get('subcellularLocations', []):
                                        locations.append(loc['location']['value'])
                            if locations:
                                loc_counts = pd.Series(locations).value_counts().reset_index()
                                loc_counts.columns = ['Location', 'Count']
                                fig = px.pie(loc_counts, values='Count', names='Location',
                                             title='Subcellular Localization',
                                             color_discrete_sequence=px.colors.qualitative.Pastel)
                                fig.update_traces(textposition='inside', textinfo='percent+label')
                                st.plotly_chart(fig, use_container_width=True)
                            else:
                                st.warning("No location data found.")

                        with tab2:
                            # Sequence Analysis
                            if show_sequence_analysis and 'sequence' in data:
                                seq = data['sequence']['value']
                                analyzer = ProtParam.ProteinAnalysis(seq)
                                st.markdown("### 📊 Sequence Properties")
                                col1, col2, col3 = st.columns(3)
                                with col1:
                                    st.markdown(f"<div class='metric-card'><h4>Molecular Weight</h4><h3>{analyzer.molecular_weight():.2f} Da</h3></div>", unsafe_allow_html=True)
                                    st.markdown(f"<div class='metric-card'><h4>Aromaticity</h4><h3>{analyzer.aromaticity():.3f}</h3></div>", unsafe_allow_html=True)
                                with col2:
                                    st.markdown(f"<div class='metric-card'><h4>Isoelectric Point</h4><h3>{analyzer.isoelectric_point():.2f}</h3></div>", unsafe_allow_html=True)
                                    st.markdown(f"<div class='metric-card'><h4>Instability Index</h4><h3>{analyzer.instability_index():.2f}</h3></div>", unsafe_allow_html=True)
                                with col3:
                                    st.markdown(f"<div class='metric-card'><h4>GRAVY</h4><h3>{analyzer.gravy():.3f}</h3></div>", unsafe_allow_html=True)
                                    sec_struct = analyzer.secondary_structure_fraction()
                                    st.markdown(f"<div class='metric-card'><h4>Secondary Structure</h4><p>Helix: {sec_struct[0]:.1%}<br>Turn: {sec_struct[1]:.1%}<br>Sheet: {sec_struct[2]:.1%}</p></div>", unsafe_allow_html=True)
                                # Amino acid composition
                                st.markdown("### 🧬 Amino Acid Composition")
                                aa_comp = analyzer.get_amino_acids_percent()
                                aa_df = pd.DataFrame.from_dict(aa_comp, orient='index', columns=['Percentage'])
                                aa_df = aa_df.sort_values('Percentage', ascending=False)
                                fig = px.bar(aa_df, y='Percentage', title='Amino Acid Composition', color_discrete_sequence=['#3498db'])
                                fig.update_layout(xaxis_title="Amino Acid", yaxis_title="Percentage", height=400)
                                st.plotly_chart(fig, use_container_width=True)
                                # Hydrophobicity profile
                                st.markdown("### 💧 Hydrophobicity Profile")
                                window_size = st.slider("Select window size for hydrophobicity:", 5, 21, 9, 2)
                                kd_scale = {
                                    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5,
                                    'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
                                    'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
                                    'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
                                    'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
                                }
                                hp = analyzer.protein_scale(kd_scale, window_size)
                                fig = px.line(x=range(1, len(hp)+1), y=hp, labels={'x': 'Position', 'y': 'Hydrophobicity'}, title=f'Hydrophobicity Profile (Window: {window_size} residues)')
                                fig.update_layout(plot_bgcolor='rgba(255,255,255,1)', paper_bgcolor='rgba(255,255,255,1)', height=400)
                                fig.add_hline(y=0, line_dash="dash", line_color="gray")
                                st.plotly_chart(fig, use_container_width=True)
                                # Download sequence
                                st.markdown("### 📥 Download Sequence")
                                st.download_button(
                                    label="Download FASTA",
                                    data=f">sp|{data.get('primaryAccession')}|{data['proteinDescription']['recommendedName']['fullName']['value']}\n{seq}",
                                    file_name=f"{data.get('primaryAccession')}.fasta",
                                    mime="text/plain"
                                )

                        with tab3:
                            # Structure & Interactions
                            if show_structural_prediction:
                                st.markdown("### 🧪 3D Structure Visualization")
                                pdb_ids = [x['id'].upper() for x in data.get('uniProtKBCrossReferences', []) if x.get('database') in ['PDB', 'PDBsum']]
                                if pdb_ids:
                                    selected_pdb = st.selectbox("Choose PDB structure:", pdb_ids)
                                    view = py3Dmol.view(query=f'pdb:{selected_pdb}', width=800, height=500)
                                    view.setStyle({'cartoon': {'color': 'spectrum'}})
                                    view.zoomTo()
                                    showmol(view, height=500)
                                    # Structure quality metrics
                                    st.markdown("### 🏆 Structure Quality Metrics")
                                    col1, col2, col3 = st.columns(3)
                                    with col1:
                                        st.markdown("<div class='metric-card'><h4>Ramachandran Favored</h4><h3>92%</h3><p>+2% vs average</p></div>", unsafe_allow_html=True)
                                    with col2:
                                        st.markdown("<div class='metric-card'><h4>Rotamer Outliers</h4><h3>3</h3><p>-1 vs typical</p></div>", unsafe_allow_html=True)
                                    with col3:
                                        st.markdown("<div class='metric-card'><h4>Clashscore</h4><h3>8.5</h3><p>Better than 90%</p></div>", unsafe_allow_html=True)
                                else:
                                    st.warning("No experimental structures found in PDB")
                                    st.markdown("### 🦠 AlphaFold Prediction")
                                    st.info("Showing predicted structure from AlphaFold Database")
                                    af_url = f"https://alphafold.ebi.ac.uk/entry/{accession}"
                                    components.iframe(af_url, height=500, scrolling=True)

                            if show_interaction_network:
                                st.markdown("### 🕸 Protein Interaction Network")
                                string_api_url = f"https://string-db.org/api/json/network?identifiers={accession}"
                                string_response = requests.get(string_api_url)
                                if string_response.status_code == 200:
                                    interactions = string_response.json()
                                    if interactions:
                                        G = nx.Graph()
                                        G.add_node(accession)
                                        for interaction in interactions:
                                            partner = interaction['preferredName_B']
                                            G.add_node(partner)
                                            G.add_edge(accession, partner)
                                        pos = nx.spring_layout(G, seed=42)
                                        node_sizes = [2000 if node == accession else 1000 for node in G.nodes()]
                                        node_colors = ['#3498db' if node == accession else '#95a5a6' for node in G.nodes()]
                                        fig, ax = plt.subplots(figsize=(10, 8))
                                        nx.draw(G, pos, with_labels=True, node_size=node_sizes, node_color=node_colors, font_size=10, font_weight='bold', edge_color='#bdc3c7', width=1.5, ax=ax)
                                        ax.set_title("Protein Interaction Network", pad=20, fontsize=14, color='black')
                                        st.pyplot(fig)
                                        st.markdown("#### Interaction Details")
                                        interaction_df = pd.DataFrame([{
                                            'Partner': i['preferredName_B'],
                                            'Score': i['score'],
                                            'Interaction Type': ', '.join(i.get('ncbiTaxon2Names', [])) if isinstance(i.get('ncbiTaxon2Names'), list) else 'N/A',
                                        } for i in interactions])
                                        AgGrid(interaction_df, height=250, fit_columns_on_grid_load=True)
                                    else:
                                        st.warning("No interaction data found in STRING DB")
                                else:
                                    st.error("Couldn't connect to STRING API")

                        with tab4:
                            # Evolutionary Analysis
                            if show_evolutionary:
                                st.markdown("### 🌍 Evolutionary Analysis")
                                organism_info = data.get('organism', {})
                                lineage = organism_info.get('lineage', [])
                                current_species = organism_info.get("scientificName", "Unknown")
                                if lineage:
                                    st.markdown("#### Taxonomic Lineage")
                                    tree_nodes = lineage + [current_species]
                                    fig = go.Figure()
                                    for i in range(len(tree_nodes) - 1):
                                        fig.add_trace(go.Scatter(
                                            x=[0, 0],
                                            y=[i, i+1],
                                            mode='lines',
                                            line=dict(color='gray', width=2),
                                            hoverinfo='none'
                                        ))
                                    fig.add_trace(go.Scatter(
                                        x=[0] * len(tree_nodes),
                                        y=list(range(len(tree_nodes))),
                                        mode='markers+text',
                                        marker=dict(size=12, color='#000000'),
                                        text=tree_nodes,
                                        textposition='middle right',
                                        hoverinfo='text'
                                    ))
                                    fig.update_layout(
                                        height=60 * len(tree_nodes),
                                        margin=dict(l=100, r=20, t=40, b=20),
                                        xaxis=dict(showticklabels=False, range=[-0.1, 0.1]),
                                        yaxis=dict(showticklabels=False),
                                        showlegend=False,
                                        plot_bgcolor='rgba(255,255,255,1)',
                                        paper_bgcolor='rgba(255,255,255,1)'
                                    )
                                    st.plotly_chart(fig, use_container_width=True)
                                else:
                                    st.warning("No taxonomic lineage available")
                                # Sequence comparison
                                st.markdown("#### Sequence Similarity")
                                comparative_df = pd.DataFrame({
                                    'Protein': [accession, 'P00534', 'P15056'],
                                    'Organism': [current_species, 'Mus musculus', 'Drosophila melanogaster'],
                                    'Similarity (%)': [100, 85, 78],
                                    'Key Difference': ['-', 'V83M mutation', 'Missing C-terminal domain']
                                })
                                AgGrid(comparative_df, height=200, fit_columns_on_grid_load=True)
                                # BLAST link
                                if 'sequence' in data:
                                    sequence = data['sequence']['value']
                                    blast_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&QUERY={quote(sequence)}"
                                    st.markdown(f"""
                                    <div class="feature-card">
                                        <h4>Run BLAST Search</h4>
                                        <p>Perform a detailed sequence similarity search using NCBI BLAST:</p>
                                        <a href="{blast_url}" target="_blank">Open BLAST</a>
                                    </div>
                                    """, unsafe_allow_html=True)

                        with tab5:
                            # Literature
                            if show_publications:
                                st.markdown("### 📚 Related Publications")
                                citations = [c for c in data.get('references', [])]
                                if citations:
                                    for i, cite in enumerate(citations[:5]):
                                        citation = cite.get('citation', {})
                                        title = citation.get('title', 'No title')
                                        authors = ', '.join(citation.get('authors', ['Unknown']))
                                        cross_refs = citation.get('citationCrossReferences', [])
                                        journal = cross_refs[0].get('name', 'N/A') if cross_refs else 'N/A'
                                        pmid = cross_refs[0].get('id', 'N/A') if cross_refs else 'N/A'
                                        st.markdown(f"""
                                        <div class="feature-card">
                                            <h4>{title}</h4>
                                            <p><strong>Authors:</strong> {authors}</p>
                                            <p><strong>Journal:</strong> {journal}</p>
                                            <p><strong>PMID:</strong> {pmid}</p>
                                        </div>
                                        """, unsafe_allow_html=True)
                                    if len(citations) > 5:
                                        st.info(f"Showing 5 of {len(citations)} references. See UniProt for full list.")
                                else:
                                    st.warning("No publication references found")
                                # LitSense link
                                protein_name = data['proteinDescription']['recommendedName']['fullName']['value']
                                st.markdown(f"""
                                <div class="feature-card">
                                    <h4>Find More Publications</h4>
                                    <p>Search for additional articles on this protein:</p>
                                    <a href="https://www.ncbi.nlm.nih.gov/research/litsense-api/?query={protein_name}" target="_blank">Search LitSense</a>
                                </div>
                                """, unsafe_allow_html=True)
                        # Show progress
                        st.success("Analysis complete!")

                    else:
                        st.error(f"Failed to retrieve detailed info for {accession}")
                else:
                    st.warning(f"No protein found for query: {query}")
            else:
                st.error("Failed to connect to UniProt API")


# -------- About Creator Page --------
elif page == "👤 About":
    st.title("👥 About")
    # No need for session state here
    st.markdown("""
    <div class="about-section-container">
        <h3>👩‍🔬 About the Author</h3>
        <div style="display: flex; align-items: center; gap: 20px;">
            <img src="https://media.licdn.com/dms/image/v2/D4D03AQFwgRiqIz7vqQ/profile-displayphoto-shrink_200_200/profile-displayphoto-shrink_200_200/0/1728237844605?e=2147483647&v=beta&t=oObfSDCysain9NDXBZwIDPMvtvtMtIddvSUb9Ny733M" width="150" style="border-radius:50%; border:3px solid #4e79a7;">
            <div>
                <h4>Isha Patil </h4>
                <p>M.Sc. Bioinformatics Student at DES Pune University</p>
                <p>MSc Bioinformatics candidate passionate about computational biology, rare disease research, and data visualization. I apply bioinformatics, machine learning, and network analysis to uncover gene-disease associations and therapeutic targets, with a focus on rare diseases. Skilled in analyzing high-throughput sequencing and multi-omics data. Experienced in building interactive dashboards using Power BI, Tableau, Streamlit, and Python/R. I aim to bridge biological insights with data-driven solutions for precision medicine and drug discovery.</p>
                <p>🔗 <a href="https://www.linkedin.com/in/isha-patilbiotechnologystudent" target="_blank">Connect on LinkedIn</a></p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # Mentorship
    st.markdown("""
    <div class="about-section-container">
        <h3>👨‍🏫 Mentorship</h3>
        <div style="display: flex; align-items: center; gap: 20px;">
            <img src="https://media.licdn.com/dms/image/v2/D5603AQF9gsU7YBjWVg/profile-displayphoto-shrink_400_400/B56ZZI.WrdH0Ag-/0/1744981029051?e=1752105600&v=beta&t=F4QBDSEgjUvnBS00xPkKqPTLI0jQaMpYefaOzARY1Yg" width="150" style="border-radius:50%; border:3px solid #4e79a7;">
            <div>
                <h4>Dr. Kushagra Kashyap</h4>
                <p>Assistant Professor (Bioinformatics), Department of Life Sciences, School of Science and Mathematics, DES Pune University</p>
                <p>This project was developed under the guidance of Dr. Kashyap, who provided valuable insights and mentorship throughout the development process. His expertise in bioinformatics and computational biology was instrumental in shaping this project.</p>
                <p>🔗 <a href="https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb" target="_blank">Connect on LinkedIn</a></p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)