console.log("toggle_subgraph.js cargado correctamente");

// Funciones de debug
if (!Array.isArray(fullData)) {
    console.error("fullData no es un array:", fullData);
    fullData = []; 
}
console.log("fullData:", fullData);

const formatIntersection = (intersection) => {
    if (!intersection) return '';
    const genes = intersection.split(',')
        .map(gene => gene.trim())
        .filter(gene => gene)
        .sort((a, b) => a.localeCompare(b));
    
    return genes.length > 5 
        ? genes.slice(0, 5).join(', ') + ', ...' 
        : genes.join(', ');
};


let activeButton = null;

document.addEventListener('DOMContentLoaded', () => {
    console.log("DOM completamente cargado");

    document.querySelectorAll('.toggle-subgraph').forEach(button => {
        console.log("Botón encontrado:", button);

        button.addEventListener('click', () => {
            console.log("Botón clickeado:", button);

            const subgraphId = button.getAttribute('data-subgraph');
            console.log("Subgraph ID:", subgraphId);

            const currentRow = button.closest('tr'); 
            const isExpanded = button.textContent === '-';

            if (isExpanded) {
                console.log("Ocultando tabla adicional");

                const additionalTable = currentRow.nextElementSibling;
                if (additionalTable && additionalTable.classList.contains('additional-table')) {
                    additionalTable.remove();
                }
                button.textContent = '+';
                button.classList.remove('btn-success');
                button.classList.add('btn-outline-primary');
                activeButton = null;
            } else {
                console.log("Mostrando tabla adicional");

                // Si hay un botón activo lo desactivamos primero
                if (activeButton && activeButton !== button) {
                    activeButton.textContent = '+';
                    activeButton.classList.remove('btn-success');
                    activeButton.classList.add('btn-outline-primary');
                    const previousTable = activeButton.closest('tr').nextElementSibling;
                    if (previousTable && previousTable.classList.contains('additional-table')) {
                        previousTable.remove();
                    }
                }

                const matchingRows = fullData.filter(row => row["subgraph_id"] === subgraphId);
                console.log("Filas coincidentes:", matchingRows);

                if (matchingRows.length > 0) {
                    const tableHTML = `
                        <tr class="additional-table">
                            <td colspan="${currentRow.cells.length}">
                                <h3 class="text-center">Subgraph members</h3>
                                <table class="table table-bordered table-hover">
                                    <thead class="table-success">
                                        <tr>
                                            <th class="text-center">Iteration</th>
                                            <th class="text-center">Cell type</th>
                                            <th class="text-center">Cluster</th>
                                            <th class="text-center">Module</th>
                                            <th class="text-center">Term id</th>
                                            <th class="text-center">Term name</th>
                                            <th class="text-center">P-value</th>
                                            <th class="text-center">Intersection</th>
                                            <th class="text-center">Length of Intersection</th>
                                            <th class="text-center">Source</th>
                                            <th class="text-center">Subgraph ID</th>
                                            <th class="text-center">Subgraph size</th>
                                            <th class="text-center">IC</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        ${matchingRows.map(row => `
                                            <tr>
                                                <td class="text-center">${row.iteration}</td>
                                                <td class="text-center">${row.cell_type.replace(/_/g, ' ')}</td>
                                                <td class="text-center">${row.cluster}</td>
                                                <td class="text-center">${row.module}</td>
                                                <td class="text-center">${row.term_id}</td>
                                                <td class="text-center">${row.term_name}</td>
                                                <td class="text-center">${row.p_value}</td>
                                                <td class="text-center">${formatIntersection(row.intersection)}</td>
                                                <td class="text-center">${row.length_intersection}</td>
                                                <td class="text-center">${row.source}</td>
                                                <td class="text-center">${row.subgraph_id}</td>
                                                <td class="text-center">${row.subgraph_size}</td>
                                                <td class="text-center">${row.IC}</td>
                                            </tr>
                                        `).join('')}
                                    </tbody>
                                </table>
                            </td>
                        </tr>
                    `;
                    currentRow.insertAdjacentHTML('afterend', tableHTML);
                } else {
                    const noResultsHTML = `
                        <tr class="additional-table">
                            <td colspan="${currentRow.cells.length}">
                                <p>No additional results found.</p>
                            </td>
                        </tr>
                    `;
                    currentRow.insertAdjacentHTML('afterend', noResultsHTML);
                }

                button.textContent = '-';
                button.classList.remove('btn-outline-primary');
                button.classList.add('btn-success');
                activeButton = button;
            }
        });
    });
});