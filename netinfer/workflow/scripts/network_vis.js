// Network visualization and interaction code for NetInfer

class NetworkVisualizer {
    constructor(containerId, data, options) {
        this.container = document.getElementById(containerId);
        this.data = data;
        this.options = Object.assign({
            nodeSize: 'degree',
            edgeWidth: 'methods_count',
            colorBy: 'taxonomy'
        }, options);
        
        this.network = null;
        this.selectedEdge = null;
        
        this.initializeNetwork();
        this.setupEventListeners();
    }
    
    initializeNetwork() {
        // Prepare nodes with size and color mapping
        const nodes = new vis.DataSet(this.data.nodes.map(node => {
            const size = this.calculateNodeSize(node);
            const color = this.getNodeColor(node);
            
            return {
                id: node.id,
                label: node.id,
                size: size,
                color: color,
                title: this.generateNodeTooltip(node)
            };
        }));
        
        // Prepare edges with width mapping
        const edges = new vis.DataSet(this.data.edges.map(edge => {
            const width = this.calculateEdgeWidth(edge);
            
            return {
                from: edge.source,
                to: edge.target,
                width: width,
                title: this.generateEdgeTooltip(edge)
            };
        }));
        
        // Network configuration
        const options = {
            nodes: {
                shape: 'circle',
                font: {
                    size: 12
                },
                borderWidth: 2
            },
            edges: {
                smooth: {
                    type: 'continuous'
                },
                color: {
                    color: '#848484',
                    highlight: '#FF0000'
                }
            },
            physics: {
                solver: 'forceAtlas2Based',
                forceAtlas2Based: {
                    gravitationalConstant: -50,
                    centralGravity: 0.01,
                    springLength: 200,
                    springConstant: 0.08
                },
                stabilization: {
                    iterations: 100
                }
            },
            interaction: {
                hover: true,
                tooltipDelay: 100,
                hideEdgesOnDrag: true,
                multiselect: false
            }
        };
        
        // Create network
        this.network = new vis.Network(
            this.container,
            { nodes: nodes, edges: edges },
            options
        );
    }
    
    setupEventListeners() {
        // Edge selection event
        this.network.on('selectEdge', params => {
            if (params.edges.length > 0) {
                const edgeId = params.edges[0];
                const edge = this.data.edges.find(e => 
                    (e.source === this.network.body.edges[edgeId].fromId &&
                     e.target === this.network.body.edges[edgeId].toId)
                );
                
                if (edge) {
                    this.selectedEdge = edge;
                    this.updateAbundancePlot(edge);
                }
            }
        });
        
        // Add filter controls event listeners
        document.getElementById('filter-methods').addEventListener('change', e => {
            this.filterByMethods(parseInt(e.target.value));
        });
        
        document.getElementById('filter-taxonomy').addEventListener('change', e => {
            this.filterByTaxonomy(e.target.value);
        });
    }
    
    calculateNodeSize(node) {
        const minSize = 10;
        const maxSize = 50;
        
        let value;
        switch (this.options.nodeSize) {
            case 'degree':
                value = node.degree;
                break;
            case 'betweenness':
                value = node.betweenness;
                break;
            case 'closeness':
                value = node.closeness;
                break;
            default:
                value = node.degree;
        }
        
        return minSize + (maxSize - minSize) * 
               (value - Math.min(...this.data.nodes.map(n => n[this.options.nodeSize]))) /
               (Math.max(...this.data.nodes.map(n => n[this.options.nodeSize])) - 
                Math.min(...this.data.nodes.map(n => n[this.options.nodeSize])));
    }
    
    calculateEdgeWidth(edge) {
        const minWidth = 1;
        const maxWidth = 10;
        
        let value;
        switch (this.options.edgeWidth) {
            case 'methods_count':
                value = edge.methods_count;
                break;
            default:
                value = edge.methods_count;
        }
        
        return minWidth + (maxWidth - minWidth) * 
               (value - Math.min(...this.data.edges.map(e => e[this.options.edgeWidth]))) /
               (Math.max(...this.data.edges.map(e => e[this.options.edgeWidth])) - 
                Math.min(...this.data.edges.map(e => e[this.options.edgeWidth])));
    }
    
    getNodeColor(node) {
        if (!node.taxonomy) {
            return '#808080';  // Default gray
        }
        
        // Extract taxonomic level for coloring
        const level = this.options.colorBy;
        const taxonomy = node.taxonomy.split(';');
        const taxonLevel = taxonomy.find(t => t.startsWith(level.charAt(0) + '__'));
        
        if (!taxonLevel) {
            return '#808080';
        }
        
        // Generate consistent color from taxonomy string
        return this.stringToColor(taxonLevel);
    }
    
    stringToColor(str) {
        let hash = 0;
        for (let i = 0; i < str.length; i++) {
            hash = str.charCodeAt(i) + ((hash << 5) - hash);
        }
        
        const c = (hash & 0x00FFFFFF)
            .toString(16)
            .toUpperCase();
        
        return '#' + '00000'.substring(0, 6 - c.length) + c;
    }
    
    generateNodeTooltip(node) {
        return `
            <strong>${node.id}</strong><br>
            Degree: ${node.degree}<br>
            Betweenness: ${node.betweenness.toFixed(4)}<br>
            Closeness: ${node.closeness.toFixed(4)}<br>
            ${node.taxonomy ? `Taxonomy: ${node.taxonomy}` : ''}
        `;
    }
    
    generateEdgeTooltip(edge) {
        return `
            ${edge.source} ↔ ${edge.target}<br>
            Methods: ${edge.methods_count}<br>
            Supporting methods: ${edge.methods}
        `;
    }
    
    updateAbundancePlot(edge) {
        // Update abundance correlation plot using Plotly
        Plotly.newPlot('abundance-plot', [{
            x: abundanceData[edge.source],
            y: abundanceData[edge.target],
            mode: 'markers',
            type: 'scatter',
            name: 'Samples'
        }], {
            title: `Abundance Correlation: ${edge.source} vs ${edge.target}`,
            xaxis: { title: edge.source },
            yaxis: { title: edge.target }
        });
    }
    
    filterByMethods(minMethods) {
        const filteredEdges = this.data.edges.filter(edge => 
            edge.methods_count >= minMethods
        );
        
        this.network.setData({
            nodes: this.data.nodes,
            edges: filteredEdges
        });
    }
    
    filterByTaxonomy(taxon) {
        if (!taxon) {
            this.network.setData({
                nodes: this.data.nodes,
                edges: this.data.edges
            });
            return;
        }
        
        const filteredNodes = this.data.nodes.filter(node =>
            node.taxonomy && node.taxonomy.includes(taxon)
        );
        
        const nodeIds = new Set(filteredNodes.map(node => node.id));
        
        const filteredEdges = this.data.edges.filter(edge =>
            nodeIds.has(edge.source) && nodeIds.has(edge.target)
        );
        
        this.network.setData({
            nodes: filteredNodes,
            edges: filteredEdges
        });
    }
}

// Initialize network visualization
const networkData = JSON.parse(document.getElementById('network-data').textContent);
const abundanceData = JSON.parse(document.getElementById('abundance-data').textContent);

const visualizer = new NetworkVisualizer('network', networkData, {
    nodeSize: 'degree',
    edgeWidth: 'methods_count',
    colorBy: 'phylum'
});