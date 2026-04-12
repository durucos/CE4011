"""
PURPOSE: 
    A 2D Frame and Truss Structural Analysis Engine using the Direct Stiffness Method.
    Developed for CE 4011 - Assignment #3.

ASSUMPTIONS:
    1. Linear elastic material behavior (Hooke's Law).
    2. Small displacement theory (First-order theory).
    3. Euler-Bernoulli beam theory (Plane sections remain plane).
    4. Neglects shear deformation and axial-bending interaction (P-Delta).
    5. Units must be consistent (e.g., Meter, KiloNewton).

INPUT: 
    A structured .txt file defining nNode, nMaterial, nMember, and nLoads.
OUTPUT: 
    1. Independent Structure Count & Stability Status.
    2. Nodal Displacement Table [m, rad].
    3. Member-End Force Table (Local System) [kN, kNm].
    4. Support Reaction Table [kN, kNm].
"""

import numpy as np
import re
import os

# =================================================================
# 1. DATA CLASSES (UML: Material, Node)
# =================================================================

class Material:
    def __init__(self, m_id, E, area, inertia):
        self.id = m_id
        self.E = E          # Elastic Modulus
        self.A = area       # Cross-sectional Area
        self.I = inertia    # Moment of Inertia

class Node:
    def __init__(self, node_id, x, y, fixity):
        self.id = node_id
        self.x = x
        self.y = y
        self.fixity = fixity  # [trX, trY, rotZ] where 1 is fixed, 0 is free
        self.dof_indices = [-1, -1, -1] # Mapping to Global Matrix

# =================================================================
# 2. ELEMENT CLASS (UML: Member)
# =================================================================

class Member:
    def __init__(self, m_id, sn, en, mat, m_type="frame"):
        self.id = m_id
        self.sn = sn        # Start Node Object
        self.en = en        # End Node Object
        self.mat = mat      # Material Object
        self.type = m_type.lower()
        self.releases = [0, 0] # [rs, re] 1 for hinge, 0 for rigid
        self.udl = 0.0

    @property
    def length(self):
        return np.sqrt((self.en.x - self.sn.x)**2 + (self.en.y - self.sn.y)**2)

    def get_matrices(self):
        """Formulates Local k, Transformation R, and Global kg."""
        L = self.length
        c = (self.en.x - self.sn.x) / L
        s = (self.en.y - self.sn.y) / L
        E, A, I = self.mat.E, self.mat.A, self.mat.I
        
        ae_l = E * A / L
        k_loc = np.zeros((6, 6))
        
        if self.type == "truss":
            k_loc[0, 0] = k_loc[3, 3] = ae_l
            k_loc[0, 3] = k_loc[3, 0] = -ae_l
        else:
            ei_l, ei_l2, ei_l3 = E*I/L, E*I/L**2, E*I/L**3
            k_loc = np.array([
                [ae_l, 0, 0, -ae_l, 0, 0],
                [0, 12*ei_l3, 6*ei_l2, 0, -12*ei_l3, 6*ei_l2],
                [0, 6*ei_l2, 4*ei_l, 0, -6*ei_l2, 2*ei_l],
                [-ae_l, 0, 0, ae_l, 0, 0],
                [0, -12*ei_l3, -6*ei_l2, 0, 12*ei_l3, -6*ei_l2],
                [0, 6*ei_l2, 2*ei_l, 0, -6*ei_l2, 4*ei_l]
            ])
            # Static Condensation for Releases
            for rel_idx, dof_idx in [(0, 2), (1, 5)]:
                if self.releases[rel_idx] == 1:
                    for i in range(6):
                        if i == dof_idx: continue
                        k_loc[i,:] -= (k_loc[i,dof_idx]/k_loc[dof_idx,dof_idx]) * k_loc[dof_idx,:]
                    k_loc[dof_idx,:], k_loc[:,dof_idx] = 0, 0

        R = np.array([[c, s, 0, 0, 0, 0], [-s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0], [0, 0, 0, -s, c, 0], [0, 0, 0, 0, 0, 1]])
        
        return k_loc, R, R.T @ k_loc @ R

    def get_fixed_end_forces(self):
        """Equivalent nodal loads for UDL on frame elements."""
        L = self.length
        f_fe = np.zeros(6)
        if self.type == "frame":
            f_fe[1] = f_fe[4] = self.udl * L / 2
            f_fe[2], f_fe[5] = self.udl * L**2 / 12, -self.udl * L**2 / 12
        return f_fe

# =================================================================
# 3. ANALYSIS ENGINE (UML: StructuralModel)
# =================================================================

class StructuralModel:
    def __init__(self):
        self.nodes, self.materials, self.members, self.nodal_loads = {}, {}, [], []
        self.U_res = None

    def analyze_connectivity(self):
        """Graph-based BFS to detect Independent Structures."""
        visited, structures, adj = set(), [], {n_id: [] for n_id in self.nodes}
        for m in self.members:
            adj[m.sn.id].append(m.en.id); adj[m.en.id].append(m.sn.id)
        for n_id in self.nodes:
            if n_id not in visited:
                comp, q = [], [n_id]
                visited.add(n_id); 
                while q:
                    curr = q.pop(0); comp.append(curr)
                    for neighbor in adj[curr]:
                        if neighbor not in visited: visited.add(neighbor); q.append(neighbor)
                structures.append(comp)
        return structures

    def parse_file(self, filename):
        """Reads .txt file using RegEx patterns."""
        def extract(line, key):
            match = re.search(rf"{key}=([\d\.-]+|frame|truss)", line)
            return match.group(1) if match else None

        if not os.path.exists(filename):
            print(f"![ERROR] File '{filename}' not found."); return

        with open(filename, 'r') as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('#'): i += 1; continue
            
            if "nNode" in line:
                count = int(line.split('=')[1]); i += 1
                while len(self.nodes) < count:
                    curr = lines[i].strip()
                    if curr and not curr.startswith('#'):
                        n_id = int(extract(curr, 'n'))
                        fix = [int(extract(curr, 'trX')), int(extract(curr, 'trY')), int(extract(curr, 'rotZ'))]
                        self.nodes[n_id] = Node(n_id, float(extract(curr, 'x')), float(extract(curr, 'y')), fix)
                    i += 1
            elif "nMaterial" in line:
                count = int(line.split('=')[1]); i += 1
                while len(self.materials) < count:
                    curr = lines[i].strip()
                    if curr and not curr.startswith('#'):
                        m_id = int(extract(curr, 'n'))
                        self.materials[m_id] = Material(m_id, float(extract(curr, 'elasticmod')), float(extract(curr, 'area')), float(extract(curr, 'inertia')))
                    i += 1
            elif "nMember" in line:
                count = int(line.split('=')[1]); i += 1
                while len(self.members) < count:
                    curr = lines[i].strip()
                    if curr and not curr.startswith('#'):
                        sn_id, en_id = int(extract(curr, 'startnode')), int(extract(curr, 'endnode'))
                        m = Member(len(self.members)+1, self.nodes[sn_id], self.nodes[en_id], self.materials[int(extract(curr, 'matProp'))], extract(curr, 'type') or 'frame')
                        m.releases = [int(extract(curr, 'rs') or 0), int(extract(curr, 're') or 0)]
                        m.udl = float(extract(curr, 'udl') or 0)
                        self.members.append(m)
                    i += 1
            elif "nLoads" in line:
                count = int(line.split('=')[1]); i += 1
                while len(self.nodal_loads) < count:
                    curr = lines[i].strip()
                    if curr and not curr.startswith('#'):
                        self.nodal_loads.append([int(extract(curr, 'nodeId')), float(extract(curr, 'Fx')), float(extract(curr, 'Fy')), float(extract(curr, 'Mz'))])
                    i += 1
            i += 1

    def solve(self):
        """Assembles Global Matrix and Solves. Includes Stability Checks."""
        structures = self.analyze_connectivity()
        print(f"\n[SYSTEM] {len(structures)} Independent Structure(s) Detected.")
        
        # Mapping DOFs
        cnt = 0
        for n_id in sorted(self.nodes.keys()):
            for j in range(3):
                if self.nodes[n_id].fixity[j] == 0:
                    self.nodes[n_id].dof_indices[j] = cnt; cnt += 1
        
        K, F = np.zeros((cnt, cnt)), np.zeros(cnt)
        for m in self.members:
            kl, R, kg = m.get_matrices()
            G = m.sn.dof_indices + m.en.dof_indices
            fg = R.T @ m.get_fixed_end_forces()
            for p in range(6):
                if G[p] != -1:
                    F[G[p]] -= fg[p]
                    for q in range(6):
                        if G[q] != -1: K[G[p], G[q]] += kg[p, q]

        for ld in self.nodal_loads:
            for q in range(3):
                dof = self.nodes[ld[0]].dof_indices[q]
                if dof != -1: F[dof] += ld[q+1]

        if cnt > 0 and np.abs(np.linalg.det(K)) < 1e-12:
            print("![ERROR] UNSTABLE SYSTEM: Singular Matrix Detected.")
            return False
        
        self.U_res = np.linalg.solve(K, F)
        self.report()
        return True

    def report(self):
        """Generates Assignment-Required Output Tables."""
        print("\n" + "="*95)
        print(f"{'CE 4011 - STRUCTURAL ANALYSIS REPORT':^95}")
        print("="*95)

        print("\n1. NODAL DISPLACEMENTS")
        print("-" * 65)
        print(f"{'Node':<8} | {'DX (m)':>16} | {'DY (m)':>16} | {'RZ (rad)':>16}")
        for n_id in sorted(self.nodes.keys()):
            n = self.nodes[n_id]
            d = [self.U_res[idx] if idx != -1 else 0.0 for idx in n.dof_indices]
            print(f"{n_id:<8} | {d[0]:16.6e} | {d[1]:16.6e} | {d[2]:16.6e}")

        print("\n2. MEMBER-END FORCES (LOCAL SYSTEM)")
        print("-" * 115)
        print(f"{'Memb':<6} | {'N1_Axial':>15} | {'N1_Shear':>15} | {'N1_Mom':>15} | {'N2_Axial':>15} | {'N2_Shear':>15} | {'N2_Mom':>15}")
        
        reactions = np.zeros(len(self.nodes) * 3)
        for m in self.members:
            kl, R, kg = m.get_matrices()
            G = m.sn.dof_indices + m.en.dof_indices
            ug = np.array([self.U_res[idx] if idx != -1 else 0.0 for idx in G])
            f_loc = kl @ (R @ ug) + m.get_fixed_end_forces()
            f_glob = R.T @ f_loc
            for p in range(6):
                nid = m.sn.id if p < 3 else m.en.id
                reactions[(nid-1)*3 + (p%3)] += f_glob[p]
            print(f"{m.id:<6} " + "".join([f"| {v:15.3f} " for v in f_loc]))

        print("\n3. SUPPORT REACTIONS")
        print("-" * 45)
        for n_id in sorted(self.nodes.keys()):
            if any(self.nodes[n_id].fixity):
                for j, lbl in enumerate(["Fx", "Fy", "Mz"]):
                    if self.nodes[n_id].fixity[j] == 1:
                        val = reactions[(n_id-1)*3+j]
                        for ld in self.nodal_loads:
                            if ld[0] == n_id: val -= ld[j+1]
                        print(f"Node {n_id} {lbl}: {val:18.4f}")
                print("-" * 45)

if __name__ == "__main__":
    model = StructuralModel()
    model.parse_file('str_input.txt')
    model.solve()