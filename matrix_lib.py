class Matrix:
    """
    Purpose: Object-Oriented Sparse Matrix library using Dictionary of Keys (DOK).
    Requirement: PHYSICAL storage of ONLY the upper triangle to minimize memory.
    Logic: Uses the property K_ij = K_ji; only keys where row <= col are stored.
    """
    def __init__(self, rows, cols, data=None):
        self.rows = rows
        self.cols = cols
        self.data = {} # Key: (row, col) where row <= col
        # Eğer başlangıçta bir liste (data) verilirse içeri aktar
        if data is not None:
            for r in range(len(data)):
                for c in range(len(data[0])):
                    self.add_value(r, c, data[r][c])

    def add_value(self, r, c, val):
        """Purpose: Assemble terms into the sparse matrix using symmetry mapping."""
        if abs(val) < 1e-18: return
        # Symmetry Filter: Always store in upper triangle (r <= c)
        if r > c: r, c = c, r
        self.data[(r, c)] = self.data.get((r, c), 0.0) + val

    def get_value(self, r, c):
        """Purpose: Retrieve value using symmetry if lower triangle (r > c) is requested."""
        if r > c: r, c = c, r
        return self.data.get((r, c), 0.0)

    def solve_sparse_system(self, F_vector):
        """
        Purpose: Solves [K]{D} = {F} using Gaussian Elimination.
        Logic: Reconstructs a working matrix from symmetry-aware storage for stable elimination.
        """
        n = self.rows
        A = [[self.get_value(i, j) for j in range(n)] for i in range(n)]
        B = F_vector[:]

        for i in range(n):
            pivot = A[i][i]
            if abs(pivot) < 1e-25: continue
            for k in range(i + 1, n):
                if abs(A[k][i]) < 1e-25: continue 
                factor = A[k][i] / pivot
                for j in range(i, n):
                    A[k][j] -= factor * A[i][j]
                B[k] -= factor * B[i]
        
        D = [0.0] * n
        for i in range(n - 1, -1, -1):
            if abs(A[i][i]) < 1e-25: continue
            sum_val = sum(A[i][j] * D[j] for j in range(i + 1, n))
            D[i] = (B[i] - sum_val) / A[i][i]
        return D

    def __repr__(self):
        """Purpose: Reconstruct full matrix for reporting and verification."""
        res = []
        for i in range(self.rows):
            row = [f"{self.get_value(i, j):12.4e}" for j in range(self.cols)]
            res.append(" ".join(row))
        return "\n".join(res)