let currentIterations = [];  // Vamos a guardamos aquí todas las iterations disponibles para el cell type

function updateClusters() {
    const cellType = document.getElementById('cell_type_filter').value;
    const clusterSelect = document.getElementById('cluster_filter');
    const iterationSelect = document.getElementById('iteration_filter');

    clusterSelect.innerHTML = '<option value="">All clusters</option>';
    iterationSelect.innerHTML = '<option value="">All iterations</option>';
    currentIterations = [];

    if (!cellType) {
        clusterSelect.disabled = true;
        iterationSelect.disabled = true;
        return;
    }

    // Actualizamos los clusters
    fetch(`/api/clusters?cell_type=${encodeURIComponent(cellType)}`)
        .then(response => {
            if (!response.ok) throw new Error('Network response was not ok');
            return response.json();
        })
        .then(clusters => {
            // Limpiamos y añadimos los clusters únicos
            clusterSelect.innerHTML = '<option value="">All clusters</option>';
            const uniqueClusters = [...new Set(clusters)];
            
            uniqueClusters.forEach(cluster => {
                const option = document.createElement('option');
                option.value = cluster;
                option.textContent = cluster;
                clusterSelect.appendChild(option);
            });
            
            clusterSelect.disabled = false;
            
            // Obtenemos y guardamos todas las iterations para el cell type
            return fetch(`/api/iterations?cell_type=${encodeURIComponent(cellType)}`);
        })
        .then(response => {
            if (!response.ok) throw new Error('Network response was not ok');
            return response.json();
        })
        .then(iterations => {
            currentIterations = iterations;
            updateIterationSelect(); 
        })
        .catch(error => {
            console.error('Error:', error);
            clusterSelect.innerHTML = '<option value="">Error loading clusters</option>';
            iterationSelect.innerHTML = '<option value="">Error loading iterations</option>';
        });
}

function updateIterationSelect() {
    const cluster = document.getElementById('cluster_filter').value;
    const cellType = document.getElementById('cell_type_filter').value;
    const iterationSelect = document.getElementById('iteration_filter');

    iterationSelect.innerHTML = '<option value="">All iterations</option>';
    iterationSelect.disabled = true;

    if (!cellType) return;

    if (cluster) {
        fetch(`/api/iterations?cell_type=${encodeURIComponent(cellType)}&cluster=${encodeURIComponent(cluster)}`)
            .then(response => {
                if (!response.ok) throw new Error('Network response was not ok');
                return response.json();
            })
            .then(iterations => {
                iterationSelect.innerHTML = '<option value="">All iterations</option>';
                
                if (iterations.length === 0) {
                    const option = document.createElement('option');
                    option.value = '';
                    option.textContent = 'No iterations found for this cluster';
                    option.disabled = true;
                    iterationSelect.appendChild(option);
                    
                    currentIterations.forEach(iteration => {
                        const opt = document.createElement('option');
                        opt.value = iteration;
                        opt.textContent = iteration;
                        iterationSelect.appendChild(opt);
                    });
                } else {
                    iterations.forEach(iteration => {
                        const option = document.createElement('option');
                        option.value = iteration;
                        option.textContent = iteration;
                        iterationSelect.appendChild(option);
                    });
                }
                iterationSelect.disabled = false;
            })
            .catch(error => {
                console.error('Error:', error);
                iterationSelect.innerHTML = '<option value="">Error loading iterations</option>';
            });
    } else {
        iterationSelect.innerHTML = '<option value="">All iterations</option>';
        currentIterations.forEach(iteration => {
            const option = document.createElement('option');
            option.value = iteration;
            option.textContent = iteration;
            iterationSelect.appendChild(option);
        });
        iterationSelect.disabled = false;
    }
}

function setupFilterDependencies() {
    const cellTypeSelect = document.getElementById('cell_type_filter');
    const clusterSelect = document.getElementById('cluster_filter');
    
    if (cellTypeSelect && clusterSelect) {
        cellTypeSelect.addEventListener('change', updateClusters);
        clusterSelect.addEventListener('change', updateIterationSelect);
        
        if (cellTypeSelect.value) {
            updateClusters();
            if (clusterSelect.value) {
                updateIterationSelect();
            }
        }
    }
}


document.addEventListener('DOMContentLoaded', function() {
    const cellTypeSelect = document.getElementById('cell_type_filter');
    const clusterSelect = document.getElementById('cluster_filter');

    if (cellTypeSelect) {
        cellTypeSelect.addEventListener('change', updateClusters);
        
        if (cellTypeSelect.value) {
            updateClusters();
        }
    }

    if (clusterSelect) {
        clusterSelect.addEventListener('change', updateIterationSelect);
    }

    setupFilterDependencies();
});