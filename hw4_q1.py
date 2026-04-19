import numpy as np
import os
import re
from abc import ABC, abstractmethod

# --- 1. PHYSICAL PROPERTIES ---
class Material:
    def __init__(self, m_id, E, alpha_thermal):
        self.id, self.E, self.alpha_thermal = m_id, E, alpha_thermal
"""
    Purpose: Encapsulates linear elastic material properties.
    Assumptions: Material is homogeneous, isotropic, and follows Hooke's Law.
    Input: E (Elastic Modulus), alpha_thermal (Thermal expansion coefficient).
    Output: Material properties for stiffness and thermal force calculations.
    """
class Section:
    def __init__(self, s_id, area, inertia, depth):
        self.id, self.area, self.inertia, self.depth = s_id, area, inertia, depth
"""
    Purpose: Defines the cross-sectional geometry of structural members.
    Assumptions: Plane sections remain plane after deformation (Bernoulli-Euler).
    Input: Area (A), Moment of Inertia (I), Section Depth (h).
    Output: Geometric properties for axial and flexural stiffness components.
    """
# --- 2. NODES & SUPPORTS ---
class Support(ABC):
    def __init__(self, settlements=None):
        self.settlements = settlements if settlements else [0.0, 0.0, 0.0]
    @abstractmethod
    def get_fixity_vector(self): pass
    """
    Purpose: Defines boundary conditions and prescribed displacements.
    Assumptions: Supports are perfectly rigid in restrained directions.
    Input: Fixity vector (bool), Settlement vector (prescribed displacements).
    Output: Identification of Free (f) and Specified (s) Degrees of Freedom.
    """

class ManualSupport(Support):
    def __init__(self, fixity, settlements=None):
        super().__init__(settlements)
        self.fixity = fixity
    def get_fixity_vector(self): return self.fixity

class Node:
    def __init__(self, n_id, x, y, support=None):
        self.id, self.x, self.y = n_id, x, y
        self.support = support
        self.dof_numbers = [-1, -1, -1]
"""
    Purpose: Defines a point in 2D space where members meet or loads are applied.
    Assumptions: Displacement compatibility is maintained at the node.
    Input: X, Y coordinates, Support object (optional).
    Output: Global DOF mapping used for matrix assembly and solution retrieval.
    """
# --- 3. MEMBER HIERARCHY ---
class Member(ABC):
    def __init__(self, m_id, nodes, material, section, rs=0, re=0):
        self.id, self.nodes = m_id, nodes # nodes is a list [Node1, Node2]
        self.material, self.section = material, section
        self.rs, self.re = rs, re
    """
    Purpose: Abstract base class for 1D structural elements.
    Assumptions: Small displacement theory; local-global coordinate transformation is linear.
    Input: Start/End Nodes, Material, Section objects.
    Output: Transformation matrix (T) and length (L).
    """

    @property
    def length(self):
        return np.sqrt((self.nodes[1].x - self.nodes[0].x)**2 + (self.nodes[1].y - self.nodes[0].y)**2)

    @abstractmethod
    def calculate_local_k(self): pass

    def get_transformation_matrix(self):
        L, sn, en = self.length, self.nodes[0], self.nodes[1]
        c, s = (en.x - sn.x) / L, (en.y - sn.y) / L
        return np.array([[c, s, 0, 0, 0, 0], [-s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                         [0, 0, 0, c, s, 0], [0, 0, 0, -s, c, 0], [0, 0, 0, 0, 0, 1]])
    def get_matrices(self):
        """Helper to get local stiffness, transformation, and global stiffness in one go"""
        kl = self.calculate_local_k()
        R = self.get_transformation_matrix()
        kg = R.T @ kl @ R
        return kl, R, kg

class FrameMember(Member):
    def calculate_local_k(self):
        L, E, A, I = self.length, self.material.E, self.section.area, self.section.inertia
        k = np.array([[E*A/L, 0, 0, -E*A/L, 0, 0],
                      [0, 12*E*I/L**3, 6*E*I/L**2, 0, -12*E*I/L**3, 6*E*I/L**2],
                      [0, 6*E*I/L**2, 4*E*I/L, 0, -6*E*I/L**2, 2*E*I/L],
                      [-E*A/L, 0, 0, E*A/L, 0, 0],
                      [0, -12*E*I/L**3, -6*E*I/L**2, 0, 12*E*I/L**3, -6*E*I/L**2],
                      [0, 6*E*I/L**2, 2*E*I/L, 0, -6*E*I/L**2, 4*E*I/L]])
        if self.rs == 1: k[2,:] = k[:,2] = 0
        if self.re == 1: k[5,:] = k[:,5] = 0
        return k
    """
    Purpose: Models a 2D beam-column element with axial and flexural stiffness.
    Assumptions: Euler-Bernoulli beam theory (shear deformations are neglected).
    Input: Inherits from Member.
    Output: 6x6 local and global stiffness matrices (k_local, k_global).
    """

class TrussMember(Member):
    def calculate_local_k(self):
        L, E, A = self.length, self.material.E, self.section.area
        k = np.zeros((6,6))
        # Sadece eksenel serbestlikler (0 ve 3)
        k[0,0] = E*A/L
        k[3,3] = E*A/L
        k[0,3] = k[3,0] = -E*A/L
        return k
    """
    Purpose: Models a pin-jointed element carrying only axial loads.
    Assumptions: Members cannot carry moments; rotation at ends is not resisted.
    Input: Inherits from Member.
    Output: 6x6 global stiffness matrix with only axial DOF contributions.
    """

# --- 4. LOAD HIERARCHY ---
class Load(ABC):
    @abstractmethod
    def calc_equivalent_forces(self, target): pass

    """
    Purpose: Base class for external actions on the structure.
    Input: Target ID (Node or Member).
    Output: Equivalent nodal force vector in global coordinates.
    """

class NodalLoad(Load):
    def __init__(self, node_id, fx, fy, mz):
        self.node_id, self.forces = node_id, np.array([fx, fy, mz])
    def calc_equivalent_forces(self, node): return self.forces

class MemberLoad(Load, ABC):
    def __init__(self, member_id): self.member_id = member_id

class DistributedLoad(MemberLoad):
    def __init__(self, member_id, w):
        super().__init__(member_id); self.w = w
    def calc_equivalent_forces(self, member):
        L = member.length
        return np.array([0, self.w*L/2, self.w*L**2/12, 0, self.w*L/2, -self.w*L**2/12])

class ThermalLoad(MemberLoad, ABC):
    def __init__(self, member_id, tu, tb):
        super().__init__(member_id); self.tu, self.tb = tu, tb
    @abstractmethod
    def get_thermal_gradient(self): pass

class MemberThermalLoad(ThermalLoad):
    def get_thermal_gradient(self): return self.tb - self.tu
    def calc_equivalent_forces(self, member):
        E, A, I, h = member.material.E, member.section.area, member.section.inertia, member.section.depth
        alpha = member.material.alpha_thermal
        dT_avg = (self.tu + self.tb) / 2.0
        f_axial = dT_avg * alpha * E * A
        m_grad = (E * I * alpha * self.get_thermal_gradient()) / h if h != 0 else 0.0
        return np.array([f_axial, 0, m_grad, -f_axial, 0, -m_grad])

    """
    Purpose: Calculates fixed-end forces due to temperature changes.
    Assumptions: Linear temperature gradient across section depth; uniform expansion along length.
    Input: T_upper, T_bottom, Member object.
    Output: Global equivalent nodal force vector (Axial force and Moment).
    """

# --- 5. STRUCTURAL MODEL ---
    """
    Purpose: The main engine for structural assembly and solution.
    Assumptions: Static equilibrium is maintained; the global stiffness matrix is non-singular.
    Input: Collections of Nodes, Members, and Loads.
    Output: Global displacement vector (U), member internal forces, and reactions.
    Procedure: Implements Matrix Partitioning to solve [Kff]{Df} = {Pf} - [Kfs]{Ds}.
    """
class StructuralModel:
    def __init__(self):
        self.nodes, self.members, self.loads = {}, [], []
        self.materials, self.sections = {}, {}

    def parse_file(self, filename):
        if not os.path.exists(filename): 
            print("File not found!"); return
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip().replace(' ', '')
                if not line or any(x in line for x in ['#', 'nNode', 'nMaterial', 'nMember', 'nLoads']): continue
                if "res=" in line:
                    nid = int(re.search(r'n=(\d+)', line).group(1))
                    res = [int(i) for i in re.search(r'res=([\d,]+)', line).group(1).split(',')]
                    sett = [float(i) for i in re.search(r'settlements=([\d\.,-]+)', line).group(1).split(',')] if "settlements=" in line else [0,0,0]
                    self.nodes[nid] = Node(nid, float(re.search(r'x=([\d\.-]+)', line).group(1)), float(re.search(r'y=([\d\.-]+)', line).group(1)), ManualSupport(res, sett))
                elif "E=" in line:
                    mid = int(re.search(r'n=(\d+)', line).group(1))
                    self.materials[mid] = Material(mid, float(re.search(r'E=([\d\.-]+)', line).group(1)), float(re.search(r'alpha=([\d\.-]+)', line).group(1)))
                    self.sections[mid] = Section(mid, float(re.search(r'area=([\d\.-]+)', line).group(1)), float(re.search(r'inertia=([\d\.-]+)', line).group(1)), float(re.search(r'depth=([\d\.-]+)', line).group(1)))
                elif "sn=" in line:
                    mid = int(re.search(r'n=(\d+)', line).group(1))
                    sn, en, mat = int(re.search(r'sn=(\d+)', line).group(1)), int(re.search(r'en=(\d+)', line).group(1)), int(re.search(r'mat=(\d+)', line).group(1))
                    m_type = re.search(r'type=(\w+)', line).group(1) if "type=" in line else "frame"
                    if m_type == "truss": self.members.append(TrussMember(mid, [self.nodes[sn], self.nodes[en]], self.materials[mat], self.sections[mat]))
                    else: self.members.append(FrameMember(mid, [self.nodes[sn], self.nodes[en]], self.materials[mat], self.sections[mat], int(re.search(r'rs=(\d+)', line).group(1) if "rs=" in line else 0), int(re.search(r're=(\d+)', line).group(1) if "re=" in line else 0)))
                elif "loadType=" in line:
                    if "Thermal" in line: self.loads.append(MemberThermalLoad(int(re.search(r'membId=(\d+)', line).group(1)), float(re.search(r'tu=([\d\.-]+)', line).group(1)), float(re.search(r'tb=([\d\.-]+)', line).group(1))))
                    elif "UDL" in line: self.loads.append(DistributedLoad(int(re.search(r'membId=(\d+)', line).group(1)), float(re.search(r'w=([\d\.-]+)', line).group(1))))
                    elif "Nodal" in line: self.loads.append(NodalLoad(int(re.search(r'nodeId=(\d+)', line).group(1)), float(re.search(r'fx=([\d\.-]+)', line).group(1)), float(re.search(r'fy=([\d\.-]+)', line).group(1)), float(re.search(r'mz=([\d\.-]+)', line).group(1))))

   
    def solve(self):

            if not self.nodes: return
                
            # 1. DOF Numbering
            cnt = 0
            for nid in sorted(self.nodes.keys()):
                fix = self.nodes[nid].support.get_fixity_vector()
                for i in range(3):
                    if fix[i] == 0: 
                        self.nodes[nid].dof_numbers[i] = cnt
                        cnt += 1
                    else: self.nodes[nid].dof_numbers[i] = -1

            Kff = np.zeros((cnt, cnt))
            self.F_settlement = np.zeros(cnt)
            self.F_external = np.zeros(cnt)
                
            # 2. Kff ve F_settlement (-Kfr * Ur) Assembly
            for m in self.members:
                kl, R, kg = m.get_matrices()
                dofs = m.nodes[0].dof_numbers + m.nodes[1].dof_numbers
                m_sett = np.array(m.nodes[0].support.settlements + m.nodes[1].support.settlements)
                
                for i in range(6):
                    p_dof = dofs[i]
                    if p_dof != -1: 
                        for j in range(6):
                            q_dof = dofs[j]
                            if q_dof != -1: 
                                Kff[p_dof, q_dof] += kg[i, j]
                            else:
                                self.F_settlement[p_dof] -= kg[i, j] * m_sett[j]

            # 3. External Forces (Pf)
            for ld in self.loads:
                if isinstance(ld, NodalLoad):
                    dofs = self.nodes[ld.node_id].dof_numbers
                    for i in range(3):
                        if dofs[i] != -1: self.F_external[dofs[i]] += ld.forces[i]
                else:
                    m = next((obj for obj in self.members if obj.id == ld.member_id), None)
                    if m:
                        fef_glob = m.get_transformation_matrix().T @ ld.calc_equivalent_forces(m)
                        m_dofs = m.nodes[0].dof_numbers + m.nodes[1].dof_numbers
                        for i in range(6):
                            if m_dofs[i] != -1: self.F_external[m_dofs[i]] += fef_glob[i]

            # 4. Uf Solution
            F_total = self.F_external + self.F_settlement
            self.U_res = np.linalg.solve(Kff, F_total)

            # 5. Total Displacement (U_total = Uf + Ur)
            max_id = max(self.nodes.keys())
            U_total = np.zeros(max_id * 3)
            for nid, node in self.nodes.items():
                for i in range(3):
                    dof = node.dof_numbers[i]
                    U_total[(nid-1)*3 + i] = self.U_res[dof] if dof != -1 else node.support.settlements[i]

            # 6. Forces and Reactions
            self.member_forces = {}
            self.reactions = np.zeros(max_id * 3)
            self.member_settlement_fixed_forces = {} # f_sett = k_loc * u_sett
            self.support_settlement_fixed_reactions = np.zeros(max_id * 3) # R_sett = kg * u_sett

            for m in self.members:
                kl, R, kg = m.get_matrices()
                m_sett = np.array(m.nodes[0].support.settlements + m.nodes[1].support.settlements)
                u_glob = np.concatenate([U_total[(m.nodes[0].id-1)*3 : (m.nodes[0].id-1)*3 + 3], 
                                        U_total[(m.nodes[1].id-1)*3 : (m.nodes[1].id-1)*3 + 3]])
                
                # --- SETTLEMENT EFFECTS ---
                f_sett_local = kl @ (R @ m_sett)
                self.member_settlement_fixed_forces[m.id] = f_sett_local
                
                f_sett_glob = kg @ m_sett
                self.support_settlement_fixed_reactions[(m.nodes[0].id-1)*3 : (m.nodes[0].id-1)*3 + 3] += f_sett_glob[0:3]
                self.support_settlement_fixed_reactions[(m.nodes[1].id-1)*3 : (m.nodes[1].id-1)*3 + 3] += f_sett_glob[3:6]

                # ---FINAL FORCES (EQL) ---
                f_local = kl @ (R @ u_glob)
                f_fixed_glob = np.zeros(6)
                for ld in self.loads:
                    if not isinstance(ld, NodalLoad) and ld.member_id == m.id:
                        f_local -= ld.calc_equivalent_forces(m)
                        f_fixed_glob += m.get_transformation_matrix().T @ ld.calc_equivalent_forces(m)
                
                self.member_forces[m.id] = f_local
                f_glob_net = (kg @ u_glob) - f_fixed_glob
                self.reactions[(m.nodes[0].id-1)*3 : (m.nodes[0].id-1)*3 + 3] += f_glob_net[0:3]
                self.reactions[(m.nodes[1].id-1)*3 : (m.nodes[1].id-1)*3 + 3] += f_glob_net[3:6]

            self._print_results(U_total, F_total)

    def _print_results(self, U_total, F_total):
        # 1. Partitioning Data
        print("\n" + "="*30 + " PARTITIONING DETAILS " + "="*30)
        print(f"{'DOF':<6} {'F_total':<12} {'Uf (Result)':<12}")
        print("-" * 65)
        for i in range(len(F_total)):
            print(f"{i:<6}{F_total[i]:>12.2f} {self.U_res[i]:>12.6e}")


        # 3. Final Member End Forces
        print("\n" + "="*20 + " FINAL MEMBER END FORCES (LOCAL) " + "="*20)
        print(f"{'MID':<6} {'FX1':<9} {'VY1':<9} {'MZ1':<9} {'FX2':<9} {'VY2':<9} {'MZ2':<9}")
        for mid, f in self.member_forces.items():
            print(f"{mid:<6} " + " ".join([f"{val:>9.2f}" for val in f]))


        print("\n" + "="*20 + " FINAL SUPPORT REACTIONS (NET) " + "="*20)
        for nid in sorted(self.nodes.keys()):
            if any(f == 1 for f in self.nodes[nid].support.get_fixity_vector()):
                r = self.reactions[(nid-1)*3 : (nid-1)*3 + 3]
                print(f"Node {nid:<3}: RX={r[0]:>10.3f}, RY={r[1]:>10.3f}, MZ={r[2]:>10.3f}")


if __name__ == "__main__":
    model = StructuralModel() 
    model.parse_file("hw4_q1.txt")
    model.solve()