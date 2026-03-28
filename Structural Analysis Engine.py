# =================================================================
# CE 4011 - Structural Analysis Software Development
# Part 1: Matrix Library
# =================================================================
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
# =================================================================
# CE 4011 - Structural Analysis Software Development
# Part 2: Analytical Engine (Member Computation & Global Assembly)
# =================================================================

def get_element_data(idx, XY, M, C):
    """
    Purpose: Computes local/global stiffness with precise theta transformation.
    Input: idx (int), XY (coordinates), M (material), C (connectivity).
    Output: k_loc, R, k_glob.
    """
    n1, n2 = int(C[idx][0]-1), int(C[idx][1]-1)
    mat_i = int(C[idx][2]-1)
    dx = XY[n2][0] - XY[n1][0]
    dy = XY[n2][1] - XY[n1][1]
    
    # Precise L calculation (Newton-Raphson)
    L_sq = dx**2 + dy**2
    L = L_sq / 2.0
    for _ in range(20):
        if L == 0: break
        L = (L + L_sq / L) / 2.0
    
    # Theta-based Cosine and Sine
    c = dx / L
    s = dy / L
    
    A_v, I_v, E_v = M[mat_i][0], M[mat_i][1], M[mat_i][2]
    EI = E_v * I_v
    EA_L = (E_v * A_v) / L

    # Local Stiffness [k_loc] (Euler-Bernoulli)
    k_loc = [[0.0 for _ in range(6)] for _ in range(6)]
    
    # Axial
    k_loc[0][0] = k_loc[3][3] = EA_L
    k_loc[0][3] = k_loc[3][0] = -EA_L
    
    # Shear/Flexure (Using L^3, L^2, L for precision)
    k_loc[1][1] = k_loc[4][4] = 12 * EI / (L**3)
    k_loc[1][4] = k_loc[4][1] = -12 * EI / (L**3)
    
    # 6EI/L^2 (The 10.67 value)
    val_6EI_L2 = 6 * EI / (L**2)
    k_loc[1][2] = k_loc[1][5] = k_loc[2][1] = k_loc[5][1] = val_6EI_L2
    k_loc[4][2] = k_loc[4][5] = k_loc[2][4] = k_loc[5][4] = -val_6EI_L2
    
    # 4EI/L (21.33) and 2EI/L (10.67)
    k_loc[2][2] = k_loc[5][5] = 4 * EI / L
    k_loc[2][5] = k_loc[5][2] = 2 * EI / L

    # Rotation Matrix [R]
    R = [[c, s, 0, 0, 0, 0], [-s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
         [0, 0, 0, c, s, 0], [0, 0, 0, -s, c, 0], [0, 0, 0, 0, 0, 1]]

    # k_glob = R^T * k_loc * R
    k_glob = [[0.0 for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(6):
            for m in range(6):
                for n in range(6):
                    k_glob[i][j] += R[m][i] * k_loc[m][n] * R[n][j]
                    
    return k_loc, R, k_glob

def main():
    try:
        with open('str_input.txt', 'r') as f:
            lines = [l.strip() for l in f if l.strip()]
        def get_v(s): return float(s.split('=')[1])
        
        ptr = 0
        nN = int(get_v(lines[ptr])); ptr += 1
        XY = [[get_v(p) for p in lines[ptr+i].split()[1:]] for i in range(nN)]; ptr += nN
        nE = int(get_v(lines[ptr])); ptr += 1
        Conn = [[get_v(p) for p in lines[ptr+i].split()[1:]] for i in range(nE)]; ptr += nE
        nM = int(get_v(lines[ptr])); ptr += 1
        Mat = [[get_v(p) for p in lines[ptr+i].split()[1:]] for i in range(nM)]; ptr += nM
        nS = int(get_v(lines[ptr])); ptr += 1
        Sup = [[int(get_v(p)) for p in lines[ptr+i].split()[1:]] for i in range(nS)]; ptr += nS
        nL = int(get_v(lines[ptr])); ptr += 1
        Load = [[int(get_v(lines[ptr+i].split()[1]))] + [float(lines[ptr+i].split()[j].split('=')[1]) for j in range(2, 5)] for i in range(nL)]

        E_raw = [[0, 0, 0] for _ in range(nN)]
        for s in Sup: E_raw[s[0]-1] = [s[1], s[2], s[3]]
        
        cnt = 1
        for i in range(nN):
            for j in range(3):
                if E_raw[i][j] == 0: E_raw[i][j] = cnt; cnt += 1
                else: E_raw[i][j] = 0
        NumEq = cnt - 1

        K_sparse = Matrix(NumEq, NumEq)
        elem_cache = []

        print("--- MEMBER LEVEL DATA ---")
        for i in range(nE):
            k_l, R, k_g = get_element_data(i, XY, Mat, Conn)
            n1, n2 = int(Conn[i][0]-1), int(Conn[i][1]-1)
            G = [E_raw[n1][0], E_raw[n1][1], E_raw[n1][2], E_raw[n2][0], E_raw[n2][1], E_raw[n2][2]]
            elem_cache.append((k_l, R, G))
            
            print(f"\nMember #{i+1}:")
            print(f"G-Vector: {G}")
            print("Local Stiffness Matrix [k']:")
            print(Matrix(6, 6, data=k_l))
            print("Global Stiffness Matrix [k_glob]:")
            print(Matrix(6, 6, data=k_g))
            
            for p in range(6):
                for q in range(6):
                    if G[p] != 0 and G[q] != 0 and G[q] >= G[p]:
                        K_sparse.add_value(G[p]-1, G[q]-1, k_g[p][q])

        F = [0.0] * NumEq
        for ld in Load:
            node_idx = int(ld[0]-1)
            for q in range(3):
                dof = E_raw[node_idx][q]
                if dof != 0: F[dof-1] += ld[q+1]

        D = K_sparse.solve_sparse_system(F)

        print("\n" + "="*60 + "\n--- STRUCTURAL LEVEL DATA ---")
        print("Global Sparse Stiffness Matrix [K]:\n", K_sparse)
        print("\nGlobal Force Vector {F}:", [f"{v:.4f}" for v in F])
        print("\nNodal Displacements {D}:", [f"{v:.6e}" for v in D])

        print("\n--- MEMBER END FORCES (Local) ---")
        for i in range(nE):
            k_l, R, G = elem_cache[i]
            d_glob = [D[idx-1] if idx != 0 else 0.0 for idx in G]
            d_loc = [sum(R[r][c] * d_glob[c] for c in range(6)) for r in range(6)]
            f_loc = [sum(k_l[r][c] * d_loc[c] for c in range(6)) for r in range(6)]
            print(f"Member #{i+1} f':", [f"{v:12.4e}" for v in f_loc])

    except Exception as e: print(f"Error: {e}")

if __name__ == "__main__": main()