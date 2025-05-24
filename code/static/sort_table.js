document.addEventListener("DOMContentLoaded", () => {
    const table = document.querySelector("table");
    const headers = table.querySelectorAll("th");
    const tbody = table.querySelector("tbody");
    

    headers.forEach((header, index) => {
        header.addEventListener("click", () => {
            const rows = Array.from(tbody.querySelectorAll("tr"));
            const isAscending = header.classList.contains("asc");

            rows.sort((rowA, rowB) => {
                const cellA = rowA.children[index].textContent.trim();
                const cellB = rowB.children[index].textContent.trim();

                return compareMixed(cellA, cellB, isAscending);
            });

            headers.forEach((h) => h.classList.remove("asc", "desc"));
            header.classList.toggle("asc", !isAscending);
            header.classList.toggle("desc", isAscending);

            tbody.innerHTML = "";
            rows.forEach((row) => tbody.appendChild(row));
        });
    });

    function compareMixed(a, b, isAscending) {
        const regex = /(\D+|\d+)/g;
        const aParts = a.match(regex) || [];
        const bParts = b.match(regex) || [];

        for (let i = 0; i < Math.max(aParts.length, bParts.length); i++) {
            const aPart = aParts[i] || "";
            const bPart = bParts[i] || "";

            if (!isNaN(aPart) && !isNaN(bPart)) {
                const aNum = parseInt(aPart, 10);
                const bNum = parseInt(bPart, 10);
                if (aNum !== bNum) {
                    return isAscending ? aNum - bNum : bNum - aNum;
                }
            } else if (aPart !== bPart) {
                return isAscending
                    ? aPart.localeCompare(bPart)
                    : bPart.localeCompare(aPart);
            }
        }

        return 0;
    }
});

console.log("sort_table.js se ha cargado correctamente!");

