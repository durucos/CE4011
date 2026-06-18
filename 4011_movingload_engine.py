import os
import sys
import threading
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
matplotlib.use("Agg")

# ==============================================================================
# 1. MATERIAL & CROSS-SECTION PROPERTIES
# ==============================================================================
class Material:
    """
    Purpose: Represents material properties of the structural element.
    Assumptions: Linear-elastic and isotropic material.
    """
    def __init__(self, m_id, E, alpha_thermal):
        # Input: m_id (identifier), E (Elastic Modulus), alpha_thermal (Coefficient of Thermal Expansion)
        self.id = m_id
        self.E = E                                      
        self.alpha_thermal = alpha_thermal              

class Section:
    """
    Purpose: Defines geometric properties of the cross-section.
    """
    def __init__(self, s_id, area, inertia, depth):
        # Input: s_id (identifier), area (Cross-sectional area), inertia (Moment of Inertia), depth (Section depth)
        self.id = s_id
        self.area = area        
        self.inertia = inertia  
        self.depth = depth      

# ==============================================================================
# 2. NODES & BOUNDARY CONDITIONS
# ==============================================================================
class Support(ABC):
    """
    Purpose: Abstract base class for support conditions.
    """
    def __init__(self, settlements=None):
        self.settlements = settlements if settlements else [0.0, 0.0, 0.0]

    @abstractmethod
    def get_fixity_vector(self):
        pass

class ManualSupport(Support):
    """
    Purpose: Defines support conditions with explicit fixity vector.
    Input: fixity (list of 1s and 0s for fixed/free DOFs), settlements (list of presribed displacements)
    """
    def __init__(self, fixity, settlements=None):
        super().__init__(settlements)
        self.fixity = fixity  

    def get_fixity_vector(self):
        return self.fixity

class Node:
    """
    Purpose: Represents a specific coordinate point in the structural model.
    """
    def __init__(self, n_id, x, y, support=None):
        self.id = n_id
        self.x = x
        self.y = y
        self.support = support if support else ManualSupport([0, 0, 0])
        self.dof_numbers = [-1, -1, -1]  

# ==============================================================================
# 3. FINITE ELEMENT MEMBERS
# ==============================================================================
class Member(ABC):
    """
    Purpose: Abstract base class for structural members connecting two nodes.
    """
    def __init__(self, m_id, nodes, material, section, rs=0, re=0):
        self.id = m_id
        self.nodes = nodes
        self.material = material
        self.section = section
        self.rs = rs
        self.re = re

    @property
    def length(self):
        return np.sqrt((self.nodes[1].x - self.nodes[0].x)**2 + (self.nodes[1].y - self.nodes[0].y)**2)

    @abstractmethod
    def calculate_local_k(self):
        pass

    def get_transformation_matrix(self):
        """
        Purpose: Generates the transformation matrix for coordinate transformation.
        Output: 6x6 numpy array.
        """
        L = self.length
        sn, en = self.nodes[0], self.nodes[1]
        c = (en.x - sn.x) / L
        s = (en.y - sn.y) / L
        return np.array([
            [ c,  s, 0,  0,  0, 0],
            [-s,  c, 0,  0,  0, 0],
            [ 0,  0, 1,  0,  0, 0],
            [ 0,  0, 0,  c,  s, 0],
            [ 0,  0, 0, -s,  c, 0],
            [ 0,  0, 0,  0,  0, 1]
        ])

    def get_matrices(self):
        """
        Purpose: Retrieves local and global stiffness matrices along with transformation matrix.
        Output: Tuple of (local stiffness, transformation matrix, global stiffness)
        """
        kl = self.calculate_local_k()
        R  = self.get_transformation_matrix()
        kg = R.T @ kl @ R
        return kl, R, kg

class FrameMember(Member):
    """
    Purpose: 2D Euler-Bernoulli frame element (3 DOF per node).
    """
    def calculate_local_k(self):
        """
        Purpose: Calculates the local stiffness matrix for a 2D frame element.
        """
        L, E, A, I = (self.length, self.material.E, self.section.area, self.section.inertia)
        k = np.array([
            [ E*A/L,   0,            0,           -E*A/L,  0,            0            ],
            [ 0,       12*E*I/L**3,  6*E*I/L**2,  0,      -12*E*I/L**3,  6*E*I/L**2   ],
            [ 0,       6*E*I/L**2,   4*E*I/L,     0,      -6*E*I/L**2,   2*E*I/L      ],
            [-E*A/L,   0,            0,           E*A/L,   0,            0            ],
            [ 0,      -12*E*I/L**3, -6*E*I/L**2,  0,       12*E*I/L**3, -6*E*I/L**2   ],
            [ 0,       6*E*I/L**2,   2*E*I/L,     0,      -6*E*I/L**2,   4*E*I/L      ]
        ])
        if self.rs == 1: k[2, :] = k[:, 2] = 0
        if self.re == 1: k[5, :] = k[:, 5] = 0
        return k

# ==============================================================================
# 4. LOAD HIERARCHY
# ==============================================================================
class Load(ABC):
    @abstractmethod
    def calc_equivalent_forces(self, target):
        pass

class NodalLoad(Load):
    """
    Purpose: Represents a concentrated force/moment applied directly to a node.
    """
    def __init__(self, node_id, fx, fy, mz):
        self.node_id = node_id
        self.forces = np.array([fx, fy, mz])
    def calc_equivalent_forces(self, node): return self.forces

class MemberLoad(Load, ABC):
    def __init__(self, member_id): self.member_id = member_id

class DistributedLoad(MemberLoad):
    """
    Purpose: Represents a uniformly distributed vertical load on a member.
    """
    def __init__(self, member_id, w):
        super().__init__(member_id)
        self.w = w
    def calc_equivalent_forces(self, member):
        L = member.length
        return np.array([0, self.w * L / 2, self.w * L**2 / 12, 0, self.w * L / 2, -self.w * L**2 / 12])

class MemberThermalLoad(MemberLoad):
    """
    Purpose: Computes fixed-end forces due to uniform temperature and thermal gradient.
    Assumptions: Linear temperature gradient through the depth.
    """
    def __init__(self, member_id, tu, tb):
        super().__init__(member_id)
        self.tu = tu 
        self.tb = tb 
    def calc_equivalent_forces(self, member):
        E, A, I = (member.material.E, member.section.area, member.section.inertia)
        h = member.section.depth
        alpha = member.material.alpha_thermal
        dT_avg = (self.tu + self.tb) / 2.0
        m_grad = (E * I * alpha * (self.tb - self.tu)) / h if h != 0 else 0.0
        f_axial = dT_avg * alpha * E * A
        return np.array([-f_axial, 0, -m_grad, f_axial, 0, m_grad])

class MemberPointLoad(MemberLoad):
    """
    Purpose: Computes fixed-end forces for a concentrated load on a member span.
    """
    def __init__(self, member_id, p, a, h_force=0.0):
        super().__init__(member_id)
        self.p = p             
        self.a = a             
        self.h_force = h_force 
    def calc_equivalent_forces(self, member):
        L = member.length
        a = self.a
        b = L - a
        V1 = self.p * (b**2) * (3*a + b) / (L**3)
        M1 = self.p * a * (b**2) / (L**2)
        V2 = self.p * (a**2) * (a + 3*b) / (L**3)
        M2 = -self.p * (a**2) * b / (L**2)
        Fx1 = -self.h_force * (b / L)
        Fx2 = -self.h_force * (a / L)
        return np.array([Fx1, V1, M1, Fx2, V2, M2])

# ==============================================================================
# 5. CORE DSM SOLVER (OPTIMIZED)
# ==============================================================================
class StructuralModel:
    """
    Purpose: Assembles the global stiffness matrix and solves [K]{U} = {F}
    using the Direct Stiffness Method (DSM).
    Assumptions: Linear structural behavior, small deformations.
    """
    def __init__(self):
        self.nodes = {}
        self.members = []
        self.loads = []
        self.springs = {}
        self.member_forces = {}
        self.last_displacements = None
        self._K_inv = None

    def clear_cache(self):
        self._K_inv = None

    def reset_model(self):
        self.loads = []
        self.member_forces = {}
        self.last_displacements = None

    def _node_at(self, x_coord):
        for nid, node in self.nodes.items():
            if abs(node.x - x_coord) < 1e-4: return node
        return None

    def solve(self, apply_settlements=True):
        """
        Purpose: Assembles load vector, applies boundary conditions, solves for displacements.
        Input: apply_settlements (bool) indicating whether to include prescribed settlements.
        Output: Dictionary mapping member IDs to their 6 internal forces.
        """
        node_ids = sorted(self.nodes.keys())
        if not node_ids: return None

        cnt = 0
        for nid in node_ids:
            fix = self.nodes[nid].support.get_fixity_vector()
            for i in range(3):
                if fix[i] == 1: self.nodes[nid].dof_numbers[i] = -1
                else: self.nodes[nid].dof_numbers[i] = cnt; cnt += 1
        if cnt == 0: return None

        if self._K_inv is None:
            Kff = np.zeros((cnt, cnt))
            for m in self.members:
                _, _, kg = m.get_matrices()
                dofs = self.nodes[int(m.nodes[0].id)].dof_numbers + self.nodes[int(m.nodes[1].id)].dof_numbers
                for i in range(6):
                    if dofs[i] != -1:
                        for j in range(6):
                            if dofs[j] != -1: Kff[dofs[i], dofs[j]] += kg[i, j]

            for x_coord, s_vals in self.springs.items():
                node = self._node_at(float(x_coord))
                if node is None: continue
                dofs = node.dof_numbers
                for i in range(3):
                    if dofs[i] != -1 and s_vals[i] > 0.0:
                        Kff[dofs[i], dofs[i]] += s_vals[i]

            Kff += np.eye(cnt) * 1e-9
            try:
                self._K_inv = np.linalg.inv(Kff)
            except: return None

        F_ext = np.zeros(cnt)

        for load in self.loads:
            if isinstance(load, NodalLoad):
                dofs = self.nodes[int(load.node_id)].dof_numbers
                for i in range(3):
                    if dofs[i] != -1: F_ext[dofs[i]] += load.forces[i]

        for load in self.loads:
            if isinstance(load, (DistributedLoad, MemberThermalLoad, MemberPointLoad)):
                m = next((obj for obj in self.members if obj.id == load.member_id), None)
                if m is None: continue
                fef_local = load.calc_equivalent_forces(m)
                _, R, _ = m.get_matrices()
                fef_global = R.T @ fef_local
                m_dofs = self.nodes[int(m.nodes[0].id)].dof_numbers + self.nodes[int(m.nodes[1].id)].dof_numbers
                for i in range(6):
                    if m_dofs[i] != -1: F_ext[m_dofs[i]] -= fef_global[i]

        if apply_settlements:
            for x_coord, s_vals in self.springs.items():
                node = self._node_at(float(x_coord))
                if node is None: continue
                dofs = node.dof_numbers
                for i in range(3):
                    if dofs[i] != -1 and s_vals[i] > 0.0:
                        if abs(node.support.settlements[i]) > 0.0:
                            F_ext[dofs[i]] -= s_vals[i] * node.support.settlements[i]

        U_res = self._K_inv @ F_ext

        max_nid = max(node_ids)
        U_total = np.zeros(max_nid * 3)
        for nid, node in self.nodes.items():
            base = (int(nid) - 1) * 3
            for i in range(3):
                if node.dof_numbers[i] != -1: U_total[base + i] = U_res[node.dof_numbers[i]]
                else: U_total[base + i] = node.support.settlements[i]

        self.last_displacements = U_total

        for m in self.members:
            kl, R, _ = m.get_matrices()
            b0 = (int(m.nodes[0].id) - 1) * 3
            b1 = (int(m.nodes[1].id) - 1) * 3
            u_el = np.concatenate([U_total[b0:b0+3], U_total[b1:b1+3]])
            f_local = kl @ (R @ u_el)
            
            for load in self.loads:
                if isinstance(load, (DistributedLoad, MemberThermalLoad, MemberPointLoad)) and load.member_id == m.id:
                    f_local += load.calc_equivalent_forces(m)
            self.member_forces[m.id] = f_local

        return self.member_forces

class EnvelopeTracker:
    """
    Purpose: Tracks the maximum and minimum envelope forces for multiple load steps.
    """
    def __init__(self, member_ids):
        self.max_envelope = {mid: np.full(6, -1e12) for mid in member_ids}
        self.min_envelope = {mid: np.full(6, 1e12) for mid in member_ids}

    def update(self, current_forces):
        for mid, forces in current_forces.items():
            if mid in self.max_envelope:
                self.max_envelope[mid] = np.maximum(self.max_envelope[mid], forces)
                self.min_envelope[mid] = np.minimum(self.min_envelope[mid], forces)

# ==============================================================================
# 9. TIME-STEP ORCHESTRATOR
# ==============================================================================
class TimeStepOrchestrator:
    """
    Purpose: Manages sequential load cases including dead, thermal, settlement, 
             and moving train loads.
    """
    def __init__(self, config_dict, progress_callback=None):
        self.cfg = config_dict
        self.progress_callback = progress_callback
        self.model = StructuralModel()
        self._assemble_mesh()
        
        self.dead_load_displacements = None
        self.thermal_displacements = None
        self.settlement_displacements = None
        
        self.live_uy_max = None
        self.live_uy_min = None

    def _assemble_mesh(self):
        """
        Purpose: Generates members and nodes based on the geometry dictionary.
        """
        self.model.clear_cache()  
        geom = self.cfg["geom"]
        sec = self.cfg["section"]
        self.model.springs = self.cfg["springs"]
        mat_obj = Material(1, sec["E"], sec["alpha"])
        sec_obj = Section(1, sec["A"], sec["I"], sec["h"])
        
        node_id = 1
        cumulative_x = 0.0
        for span_idx, span_L in enumerate(geom["span_lengths"]):
            n_elements = int(max(1, round(span_L * geom["elems_per_m"])))
            dx = span_L / n_elements
            for k in range(0 if span_idx == 0 else 1, n_elements + 1):
                node_x = cumulative_x + k * dx
                supp = ManualSupport([0, 0, 0]) 
                for s_x, s_vals in self.cfg["settlements"].items():
                    if abs(node_x - float(s_x)) < 1e-4: supp.settlements = s_vals
                self.model.nodes[node_id] = Node(node_id, node_x, 0.0, supp)
                node_id += 1
            cumulative_x += span_L
            
        for i in range(len(self.model.nodes) - 1):
            n1 = self.model.nodes[i + 1]
            n2 = self.model.nodes[i + 2]
            self.model.members.append(FrameMember(i + 1, [n1, n2], mat_obj, sec_obj))

    def _inject_permanent_loads(self):
        """
        Purpose: Adds self-weight and externally defined distributed/point static loads.
        """
        if self.cfg["section"]["sw_flag"] == 1:
            sw_mag = self.cfg["section"]["A"] * self.cfg["section"]["gamma"]
            for m in self.model.members: self.model.loads.append(DistributedLoad(m.id, sw_mag))
        
        for el in self.cfg["extra_loads"]:
            if el["type"] == "DIST":
                for m in self.model.members:
                    if min(m.nodes[1].x, el["end_x"]) > max(m.nodes[0].x, el["start_x"]):
                        self.model.loads.append(DistributedLoad(m.id, el["mag"]))
            elif el["type"] == "POINT":
                for m in self.model.members:
                    if m.nodes[0].x <= el["start_x"] <= m.nodes[1].x:
                        a = el["start_x"] - m.nodes[0].x
                        if a < 1e-5: a = 1e-5
                        if a > m.length - 1e-5: a = m.length - 1e-5
                        self.model.loads.append(MemberPointLoad(m.id, el["mag"], a, h_force=0.0))
                        break

    def solve_isolated_dead(self):
        self.model.reset_model()
        self._inject_permanent_loads()
        res = self.model.solve(apply_settlements=False) 
        self.dead_load_displacements = np.copy(self.model.last_displacements)
        return res

    def solve_isolated_thermal(self):
        self.model.reset_model()
        t_cfg = self.cfg["thermal"]
        for m in self.model.members:
            self.model.loads.append(MemberThermalLoad(m.id, t_cfg["t_top"], t_cfg["t_bot"]))
        res = self.model.solve(apply_settlements=False)
        self.thermal_displacements = np.copy(self.model.last_displacements)
        return res

    def solve_isolated_settlement(self):
        self.model.reset_model()
        res = self.model.solve(apply_settlements=True)
        self.settlement_displacements = np.copy(self.model.last_displacements)
        return res

    def execute_simulation_march(self):
        """
        Purpose: Advances moving train loads step-by-step and extracts force envelopes.
        """
        total_L = self.cfg["geom"]["total_L"]
        dx_step = self.cfg["sim"]["dx"]
        axles = self.cfg["sim"]["loads"]
        spacing = self.cfg["sim"]["spacing"]
        acc = self.cfg["sim"]["acc"] 
        
        member_ids = [m.id for m in self.model.members]
        self.envelope = EnvelopeTracker(member_ids)
        
        num_dofs = len(self.model.nodes) * 3
        self.live_uy_max = np.zeros(num_dofs)
        self.live_uy_min = np.zeros(num_dofs)
        
        current_nose = 0.0
        max_dist = total_L + 20.0
        
        last_percent = -1  
        
        while current_nose <= max_dist:
            self.model.reset_model()
            for i, P in enumerate(axles):
                x_axle = current_nose - i * spacing
                if 0.0 <= x_axle <= total_L:
                    H = (P / 9.81) * acc 
                    for m in self.model.members:
                        if m.nodes[0].x <= x_axle <= m.nodes[1].x:
                            a = x_axle - m.nodes[0].x
                            if a < 1e-5: a = 1e-5
                            if a > m.length - 1e-5: a = m.length - 1e-5
                            self.model.loads.append(MemberPointLoad(m.id, P, a, h_force=H))
                            break
            
            step_results = self.model.solve(apply_settlements=False)
            
            if step_results:
                self.envelope.update(step_results)
                if self.model.last_displacements is not None:
                    self.live_uy_max = np.maximum(self.live_uy_max, self.model.last_displacements)
                    self.live_uy_min = np.minimum(self.live_uy_min, self.model.last_displacements)
            
            if self.progress_callback:
                percent = int((current_nose / max_dist) * 100)
                if percent > last_percent:
                    self.progress_callback(percent)
                    last_percent = percent
                    
            current_nose += dx_step
            
        if self.progress_callback: self.progress_callback(100)
        return self.envelope

# ==============================================================================
# VISUALIZER
# ==============================================================================
class PostProcessVisualizer:
    """
    Purpose: Provides static methods for generating analysis result plots.
    """
    @staticmethod
    def plot_schematic_model(data_deck):
        spans = data_deck["geom"]["span_lengths"]
        supports = data_deck["springs"].keys()
        total_L = data_deck["geom"]["total_L"]
        extra_loads = data_deck.get("extra_loads", [])

        fig, ax = plt.subplots(figsize=(10, 3))
        
        cumulative = 0
        colors = ["#7f8c8d", "#95a5a6"]
        
        if total_L > 0:
            for idx, span_L in enumerate(spans):
                ax.plot([cumulative, cumulative + span_L], [0, 0], color=colors[idx % 2], linewidth=8, solid_capstyle='butt')
                ax.text(cumulative + span_L/2, -0.3, f"Span {idx+1}\nL={span_L}m", ha='center', fontsize=9, fontweight='bold', color="#2c3e50")
                cumulative += span_L
        else:
            ax.text(5, 0, "Please enter span length...", ha='center', fontsize=12, color='gray')
            total_L = 10 
        
        for sup_x in supports:
            ax.plot(sup_x, 0, marker='^', markersize=14, color='#c0392b') 
            ax.text(sup_x, -0.7, f"Support\nx={sup_x}m", ha='center', fontsize=8, fontweight='bold')

        hw = total_L * 0.015 if total_L > 0 else 0.5  
        
        for load in extra_loads:
            l_type = load["type"]
            mag = load["mag"]
            st_x = load["start_x"]
            en_x = load["end_x"]
            
            if l_type == "DIST" and mag != 0:
                if en_x < st_x: st_x, en_x = en_x, st_x
                if en_x == st_x: continue 
                
                ax.fill_between([st_x, en_x], 0.15, 1.0, color='dodgerblue', alpha=0.2)
                ax.plot([st_x, en_x], [1.0, 1.0], color='dodgerblue', lw=2)
                
                num_arrows = max(3, int((en_x - st_x) / 4) + 1)
                x_arr = np.linspace(st_x, en_x, num_arrows)
                for xa in x_arr:
                    ax.arrow(xa, 1.0, 0, -0.65, head_width=hw*0.8, head_length=0.2, fc='dodgerblue', ec='dodgerblue')
                
                ax.text((st_x + en_x)/2, 1.15, f"w = {mag} kN/m", ha='center', color='dodgerblue', fontweight='bold', fontsize=9)

            elif l_type == "POINT" and mag != 0:
                ax.arrow(st_x, 1.6, 0, -1.2, head_width=hw, head_length=0.3, fc='purple', ec='purple', lw=2.5)
                ax.text(st_x, 1.75, f"P = {mag} kN", ha='center', color='purple', fontweight='bold', fontsize=9)

        max_x = max([total_L] + list(supports)) if supports else total_L
        
        ax.set_ylim(-1.5, 2.5) 
        ax.set_xlim(-5, max_x + 5)
        ax.set_title("Live Bridge Schematic Preview", fontweight='bold', color="#2C3E50")
        ax.set_xlabel("Bridge Length (m)")
        ax.set_yticks([])
        ax.grid(True, linestyle="--", alpha=0.5)
        plt.tight_layout()
        return fig

    @staticmethod
    def plot_isolated_live_forces(model, envelope, support_xs, save_path=None):
        sorted_members = sorted(model.members, key=lambda m: m.nodes[0].x)
        x_pts, max_m, min_m, max_v, min_v, max_n, min_n = [], [], [], [], [], [], []
        for m in sorted_members:
            f_max, f_min = envelope.max_envelope[m.id], envelope.min_envelope[m.id]
            x_pts.extend([m.nodes[0].x, m.nodes[1].x])
            max_m.extend([f_max[2], -f_min[5]]); min_m.extend([f_min[2], -f_max[5]])
            max_v.extend([f_max[1], -f_min[4]]); min_v.extend([f_min[1], -f_max[4]])
            max_n.extend([-f_min[0], f_max[3]]); min_n.extend([-f_max[0], f_min[3]])
            
        return PostProcessVisualizer._plot_core_diagrams(model, x_pts, max_m, min_m, max_v, min_v, max_n, min_n, support_xs, "ISOLATED MOVING LOAD ENVELOPE", save_path)

    @staticmethod
    def _plot_core_diagrams(model, x_pts, max_m, min_m, max_v, min_v, max_n, min_n, support_xs, title, save_path=None):
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        
        ax1.plot(x_pts, max_m, color='darkred', label='Max Mz', lw=1.5)
        ax1.plot(x_pts, min_m, color='navy', label='Min Mz', lw=1.5)
        ax1.set_title(f"{title} - Bending Moment (Mz)", fontsize=11, fontweight='bold')
        ax1.set_ylabel("Moment [kNm]", fontsize=10)
        ax1.grid(True, linestyle=':', alpha=0.6)
        PostProcessVisualizer._annotate_peaks(ax1, x_pts, max_m, "MAX")
        PostProcessVisualizer._annotate_peaks(ax1, x_pts, min_m, "MIN")
        PostProcessVisualizer._annotate_supports(ax1, x_pts, max_m, support_xs)
        PostProcessVisualizer._annotate_supports(ax1, x_pts, min_m, support_xs)

        ax2.plot(x_pts, max_v, color='crimson', label='Max Vy', lw=1.5)
        ax2.plot(x_pts, min_v, color='royalblue', label='Min Vy', lw=1.5)
        ax2.set_title(f"{title} - Shear Force (Vy)", fontsize=11, fontweight='bold')
        ax2.set_ylabel("Shear [kN]", fontsize=10)
        ax2.grid(True, linestyle=':', alpha=0.6)
        PostProcessVisualizer._annotate_peaks(ax2, x_pts, max_v, "MAX")
        PostProcessVisualizer._annotate_peaks(ax2, x_pts, min_v, "MIN")
        PostProcessVisualizer._annotate_supports(ax2, x_pts, max_v, support_xs)
        PostProcessVisualizer._annotate_supports(ax2, x_pts, min_v, support_xs)
        
        ax3.plot(x_pts, max_n, color='darkgreen', label='Max Fx (Tension)', lw=1.5)
        ax3.plot(x_pts, min_n, color='darkgoldenrod', label='Min Fx (Compression)', lw=1.5)
        ax3.set_title(f"{title} - Axial Force (Fx)", fontsize=11, fontweight='bold')
        ax3.set_xlabel("Bridge Length [m]", fontsize=10)
        ax3.set_ylabel("Axial [kN]", fontsize=10)
        ax3.grid(True, linestyle=':', alpha=0.6)
        PostProcessVisualizer._annotate_peaks(ax3, x_pts, max_n, "MAX")
        PostProcessVisualizer._annotate_peaks(ax3, x_pts, min_n, "MIN")
        PostProcessVisualizer._annotate_supports(ax3, x_pts, max_n, support_xs)
        PostProcessVisualizer._annotate_supports(ax3, x_pts, min_n, support_xs)
        
        plt.tight_layout()
        if save_path: plt.savefig(save_path, bbox_inches='tight', dpi=150)
        return fig 

    @staticmethod
    def _annotate_peaks(ax, x, y, label):
        if len(y) == 0: return
        if label == "MAX": 
            idx = np.argmax(y); color = "red"; offset = 15
        else: 
            idx = np.argmin(y); color = "blue"; offset = -20
        ax.plot(x[idx], y[idx], 'o', color=color, markersize=6)
        ax.annotate(f"{y[idx]:.1f}", xy=(x[idx], y[idx]), 
                    xytext=(0, offset), textcoords='offset points', 
                    ha='center', fontsize=9, fontweight='bold', color=color,
                    bbox=dict(boxstyle="round", fc="white", ec=color, alpha=0.8))

    @staticmethod
    def _annotate_supports(ax, x, y, support_xs):
        if len(y) == 0: return
        for sx in support_xs:
            idx = np.argmin(np.abs(np.array(x) - sx))
            ax.plot(x[idx], y[idx], 'ko', markersize=4)
            ax.annotate(f"{y[idx]:.1f}", xy=(x[idx], y[idx]), 
                        xytext=(0, -15), textcoords='offset points', 
                        ha='center', fontsize=8, color='black',
                        bbox=dict(boxstyle="round", fc="white", ec="black", alpha=0.5))

    @staticmethod
    def plot_single_case_forces(model, forces, support_xs, title, save_path=None):
        if not forces: return None
        sorted_members = sorted(model.members, key=lambda m: m.nodes[0].x)
        x_pts, M, V, N = [], [], [], []
        for m in sorted_members:
            f = forces[m.id]
            x_pts.extend([m.nodes[0].x, m.nodes[1].x])
            M.extend([f[2], -f[5]])
            V.extend([f[1], -f[4]])
            N.extend([-f[0], f[3]])
        return PostProcessVisualizer._plot_core_diagrams(model, x_pts, M, M, V, V, N, N, support_xs, f"{title} FORCE DIAGRAMS", save_path)

    @staticmethod
    def plot_single_case_deflection(model, displacements, support_xs, title, line_color="green", save_path=None, direction="y"):
        if displacements is None: return None
        x_pts, u_pts = [], []
        for nid in sorted(model.nodes.keys()):
            node = model.nodes[nid]
            base = (int(node.id) - 1) * 3
            x_pts.append(node.x)
            if direction == "x": u_pts.append(displacements[base] * 1000)      
            else: u_pts.append(displacements[base + 1] * 1000)  
        
        fig = plt.figure(figsize=(10, 4))
        plt.plot(x_pts, u_pts, color=line_color, lw=2, label=title)
        
        idx_max = np.argmax(np.array(u_pts))
        idx_min = np.argmin(np.array(u_pts))
        plt.annotate(f"Max: {u_pts[idx_max]:.2f} mm", xy=(x_pts[idx_max], u_pts[idx_max]),
                     xytext=(0, 15), textcoords='offset points', ha='center', fontsize=9, fontweight='bold', color='black')
        plt.annotate(f"Min: {u_pts[idx_min]:.2f} mm", xy=(x_pts[idx_min], u_pts[idx_min]),
                     xytext=(0, -20), textcoords='offset points', ha='center', fontsize=9, fontweight='bold', color='black')
        
        for sx in support_xs:
            idx = np.argmin(np.abs(np.array(x_pts) - sx))
            plt.plot(x_pts[idx], u_pts[idx], 'ko', markersize=4)

        dir_text = "Horizontal Displacement (Ux)" if direction == "x" else "Vertical Deflection (Uy)"
        plt.title(f"{title} - {dir_text}", fontsize=12, fontweight='bold')
        plt.xlabel("Bridge Length [m]", fontsize=10)
        plt.ylabel(f"{dir_text} [mm]", fontsize=10)
        plt.grid(True, linestyle='--')
        plt.legend()
        plt.tight_layout()
        if save_path: plt.savefig(save_path, bbox_inches='tight', dpi=150)
        return fig

    @staticmethod
    def plot_isolated_live_deflection(model, live_uy_max, live_uy_min, support_xs, save_path=None):
        if live_uy_max is None or live_uy_min is None: return None
        x_pts, uy_max, uy_min = [], [], []
        for nid in sorted(model.nodes.keys()):
            node = model.nodes[nid]
            base = (int(node.id) - 1) * 3
            x_pts.append(node.x)
            uy_max.append(live_uy_max[base + 1] * 1000)
            uy_min.append(live_uy_min[base + 1] * 1000)

        fig = plt.figure(figsize=(10, 4))
        plt.plot(x_pts, uy_max, color='royalblue', lw=1.5, label='Max Envelope (Upward)')
        plt.plot(x_pts, uy_min, color='red', lw=2.5, label='Min Envelope (Downward)')
        plt.plot(x_pts, np.zeros_like(x_pts), color='black', lw=1, linestyle='--') 
        
        idx_max = np.argmax(uy_max)
        idx_min = np.argmin(uy_min)
        plt.annotate(f"Max: {uy_max[idx_max]:.2f} mm", xy=(x_pts[idx_max], uy_max[idx_max]),
                     xytext=(0, 15), textcoords='offset points', ha='center', fontsize=9, fontweight='bold', color='royalblue')
        plt.annotate(f"Min: {uy_min[idx_min]:.2f} mm", xy=(x_pts[idx_min], uy_min[idx_min]),
                     xytext=(0, -20), textcoords='offset points', ha='center', fontsize=9, fontweight='bold', color='red')
                     
        for sx in support_xs:
            idx = np.argmin(np.abs(np.array(x_pts) - sx))
            plt.plot(x_pts[idx], uy_min[idx], 'ko', markersize=4)

        plt.title("ISOLATED MOVING LOAD DEFLECTION ENVELOPE", fontsize=12, fontweight='bold')
        plt.xlabel("Bridge Length [m]", fontsize=10)
        plt.ylabel("Vertical Deflection (Uy) [mm]", fontsize=10)
        plt.legend()
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        if save_path: plt.savefig(save_path, bbox_inches='tight', dpi=150)
        return fig

    @staticmethod
    def plot_combined_forces(model, combined_max, combined_min, support_xs, save_path=None):
        sorted_members = sorted(model.members, key=lambda m: m.nodes[0].x)
        x_pts, max_m, min_m, max_v, min_v, max_n, min_n = [], [], [], [], [], [], []
        for m in sorted_members:
            f_max, f_min = combined_max[m.id], combined_min[m.id]
            x_pts.extend([m.nodes[0].x, m.nodes[1].x])
            max_m.extend([f_max[2], -f_min[5]]); min_m.extend([f_min[2], -f_max[5]])
            max_v.extend([f_max[1], -f_min[4]]); min_v.extend([f_min[1], -f_max[4]])
            max_n.extend([-f_min[0], f_max[3]]); min_n.extend([-f_max[0], f_min[3]])
            
        return PostProcessVisualizer._plot_core_diagrams(model, x_pts, max_m, min_m, max_v, min_v, max_n, min_n, support_xs, "COMBINED FORCE ENVELOPE (G+T+S+Q)", save_path)

    @staticmethod
    def plot_combined_deflection(orch, support_xs, save_path=None):
        model = orch.model
        x_pts, uy_base, uy_max_env, uy_min_env = [], [], [], []
        
        base_disp = np.zeros(len(model.nodes)*3)
        if orch.dead_load_displacements is not None: base_disp += orch.dead_load_displacements
        if orch.thermal_displacements is not None: base_disp += orch.thermal_displacements
        if orch.settlement_displacements is not None: base_disp += orch.settlement_displacements

        for nid in sorted(model.nodes.keys()):
            node = model.nodes[nid]
            base = (int(node.id) - 1) * 3
            x_pts.append(node.x)
            
            d_base = base_disp[base + 1] * 1000
            d_live_max = (orch.live_uy_max[base + 1] * 1000) if orch.live_uy_max is not None else 0
            d_live_min = (orch.live_uy_min[base + 1] * 1000) if orch.live_uy_min is not None else 0
            
            uy_base.append(d_base)
            uy_max_env.append(d_base + d_live_max)
            uy_min_env.append(d_base + d_live_min)

        fig = plt.figure(figsize=(10, 4))
        plt.plot(x_pts, uy_base, color='black', lw=1.5, linestyle='--', label='Permanent Load Baseline (G+T+S)')
        plt.plot(x_pts, uy_max_env, color='royalblue', lw=1.5, label='Max Combined Deflection')
        plt.plot(x_pts, uy_min_env, color='red', lw=2.5, label='Min Combined Deflection')
        
        idx_max = np.argmax(uy_max_env)
        idx_min = np.argmin(uy_min_env)
        plt.annotate(f"Max: {uy_max_env[idx_max]:.2f} mm", xy=(x_pts[idx_max], uy_max_env[idx_max]),
                     xytext=(0, 15), textcoords='offset points', ha='center', fontsize=9, fontweight='bold', color='royalblue')
        plt.annotate(f"Min: {uy_min_env[idx_min]:.2f} mm", xy=(x_pts[idx_min], uy_min_env[idx_min]),
                     xytext=(0, -20), textcoords='offset points', ha='center', fontsize=9, fontweight='bold', color='red')
                     
        for sx in support_xs:
            idx = np.argmin(np.abs(np.array(x_pts) - sx))
            plt.plot(x_pts[idx], uy_min_env[idx], 'ko', markersize=4)
                     
        plt.title("COMBINED DEFLECTION ENVELOPE (G+T+S+Q)", fontsize=12, fontweight='bold')
        plt.xlabel("Bridge Length [m]", fontsize=10)
        plt.ylabel("Combined Deflection (Uy) [mm]", fontsize=10)
        plt.legend()
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        if save_path: plt.savefig(save_path, bbox_inches='tight', dpi=150)
        return fig

def create_disp_dataframe(model, disp_array):
    if disp_array is None: return None
    data = []
    for nid in sorted(model.nodes.keys()):
        node = model.nodes[nid]
        base = (int(node.id) - 1) * 3
        data.append({
            "Node ID": node.id,
            "X Coord (m)": round(node.x, 3),
            "Ux (mm)": round(disp_array[base] * 1000, 4),      
            "Uy (mm)": round(disp_array[base + 1] * 1000, 4),  
            "Rz (rad)": round(disp_array[base + 2], 6)         
        })
    df = pd.DataFrame(data)
    df.set_index("Node ID", inplace=True)
    return df

# ==============================================================================
# GUI APPLICATION (TKINTER)
# ==============================================================================
class BridgeGUI:
    """
    Purpose: Provides the graphical user interface for the analysis engine.
    Input: User-defined parameters via Tkinter widgets.
    Output: Modifies and initiates the backend process to yield structural responses.
    """
    def __init__(self, root):
        self.root = root
        self.root.title("CE 4011 - Structural Computational Engine")
        self.root.geometry("1500x900") 
        
        self.preview_timer = None 
        
        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.style.configure("TButton", font=("Segoe UI", 10, "bold"), padding=5)
        self.style.configure("Accent.TButton", foreground="white", background="#2980B9", font=("Segoe UI", 10, "bold"))
        self.style.configure("TLabel", font=("Segoe UI", 9))
        self.style.configure("Header.TLabel", font=("Segoe UI", 10, "bold"), foreground="#2C3E50")
        
        self.cfg = {}
        self._setup_ui()
        self._set_default_values() 
        
    def _display_df_to_frame(self, frame, df):
        """
        Purpose: Displays the DataFrame content in a Treeview table.
        Input: Tkinter Frame and a Pandas DataFrame.
        """
        for widget in frame.winfo_children(): widget.destroy()
        if df is None: return

        # Add scrollbars so the table is scrollable if large
        x_scroll = ttk.Scrollbar(frame, orient=tk.HORIZONTAL)
        y_scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL)
        
        # Create Treeview
        tree = ttk.Treeview(frame, show='headings', 
                            xscrollcommand=x_scroll.set, yscrollcommand=y_scroll.set)
        
        # Configure columns
        df_reset = df.reset_index()
        tree["columns"] = list(df_reset.columns)
        
        for col in df_reset.columns:
            tree.heading(col, text=col)
            tree.column(col, width=100)
            
        # Populate data
        for _, row in df_reset.iterrows():
            tree.insert("", "end", values=list(row))
            
        x_scroll.config(command=tree.xview)
        y_scroll.config(command=tree.yview)
        x_scroll.pack(side=tk.BOTTOM, fill=tk.X)
        y_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        tree.pack(fill=tk.BOTH, expand=True)

    def _setup_ui(self):
        control_frame = ttk.Frame(self.root, padding=10)
        control_frame.pack(side=tk.TOP, fill=tk.X)
        
        self.btn_refresh = ttk.Button(control_frame, text="👁️ Preview Bridge Model", command=self.refresh_model_view)
        self.btn_refresh.pack(side=tk.LEFT, padx=10)
        self.btn_new = ttk.Button(control_frame, text="🆕 New Geometry", command=self.define_new_geometry)
        self.btn_new.pack(side=tk.LEFT, padx=5)

        self.btn_run = ttk.Button(control_frame, text="▶ RUN ANALYSIS", command=self.start_analysis_thread, style="Accent.TButton")
        self.btn_run.pack(side=tk.LEFT, padx=10)
        self.btn_defaults = ttk.Button(control_frame, text="🔄 Load Defaults", command=self.load_defaults_command)
        self.btn_defaults.pack(side=tk.LEFT, padx=5)
        
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(control_frame, orient="horizontal", length=300, mode="determinate", variable=self.progress_var)
        self.progress_bar.pack(side=tk.LEFT, padx=20, pady=5)
        self.status_label = ttk.Label(control_frame, text="Status: Ready", font=("Segoe UI", 9, "italic"))
        self.status_label.pack(side=tk.LEFT, padx=5, pady=5)

        
        self.paned_window = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        self.paned_window.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.input_frame = ttk.Frame(self.paned_window, width=380, relief=tk.SUNKEN)
        self.paned_window.add(self.input_frame, weight=1)
        
        input_notebook = ttk.Notebook(self.input_frame)
        input_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.tab_geom = ttk.Frame(input_notebook)
        self.tab_bound = ttk.Frame(input_notebook)
        self.tab_loads = ttk.Frame(input_notebook)
        self.tab_train = ttk.Frame(input_notebook)
        
        input_notebook.add(self.tab_geom, text="Geometry & Section")
        input_notebook.add(self.tab_bound, text="Supports (Springs)")
        input_notebook.add(self.tab_loads, text="Static & Thermal")
        input_notebook.add(self.tab_train, text="Moving Load")
        
        self._build_geom_tab()
        self._build_bound_tab()
        self._build_loads_tab()
        self._build_train_tab()

        self.results_frame = ttk.Frame(self.paned_window)
        self.paned_window.add(self.results_frame, weight=4)
        
        self.res_notebook = ttk.Notebook(self.results_frame)
        self.res_notebook.pack(fill=tk.BOTH, expand=True)
        
        self.tab_model_view = ttk.Frame(self.res_notebook)
        self.res_notebook.add(self.tab_model_view, text="Model View")

        self.frames = {
            "DEAD": self._create_res_sub_tabs(self.res_notebook, "Dead Load"),
            "THERMAL": self._create_res_sub_tabs(self.res_notebook, "Thermal"),
            "SETTLEMENT": self._create_res_sub_tabs(self.res_notebook, "Settlement"),
            "LIVE": self._create_res_sub_tabs(self.res_notebook, "Moving Load"),
            "COMBINED": self._create_res_sub_tabs(self.res_notebook, "Combined Envelope")
        }

    def _create_res_sub_tabs(self, parent_notebook, title):
        main_tab = ttk.Notebook(parent_notebook)
        parent_notebook.add(main_tab, text=title)
        
        frame_forces_plot = ttk.Frame(main_tab)
        frame_forces_table = ttk.Frame(main_tab)
        frame_deflections_plot = ttk.Frame(main_tab)
        frame_deflections_table = ttk.Frame(main_tab)
        
        main_tab.add(frame_forces_plot, text="Forces (Plot)")
        main_tab.add(frame_forces_table, text="Forces (Table)")
        main_tab.add(frame_deflections_plot, text="Deflections (Plot)")
        main_tab.add(frame_deflections_table, text="Deflections (Table)")
        
        return {
            "forces": {"plot": frame_forces_plot, "table": frame_forces_table},
            "deflections": {"plot": frame_deflections_plot, "table": frame_deflections_table}
        }

    def _build_geom_tab(self):
        ttk.Label(self.tab_geom, text="Span Lengths (m) - 10 Spans Max:", style="Header.TLabel").grid(row=0, column=0, columnspan=2, pady=5, sticky="w")
        self.span_vars = []
        for i in range(10):
            ttk.Label(self.tab_geom, text=f"Span {i+1}:").grid(row=i+1, column=0, sticky="e")
            var = tk.StringVar()
            var.trace_add("write", lambda *args: self._on_span_change())
            ttk.Entry(self.tab_geom, textvariable=var, width=12).grid(row=i+1, column=1, pady=2, sticky="w")
            self.span_vars.append(var)
            
        ttk.Label(self.tab_geom, text="Section Properties:", style="Header.TLabel").grid(row=11, column=0, columnspan=2, pady=10, sticky="w")
        
        props = [("Area (m2):", "A"), ("Inertia (m4):", "I"), ("Depth h (m):", "h"),
                 ("Elastic Modulus E (kN/m2):", "E"), ("Thermal Alpha (1/°C):", "alpha"),
                 ("Unit Weight (kN/m3):", "gamma"), ("Elements per meter:", "epm")]
        self.sec_vars = {}
        for idx, (lbl, key) in enumerate(props):
            ttk.Label(self.tab_geom, text=lbl).grid(row=12+idx, column=0, sticky="e")
            var = tk.StringVar()
            ttk.Entry(self.tab_geom, textvariable=var, width=15).grid(row=12+idx, column=1, pady=2, sticky="w")
            self.sec_vars[key] = var

        self.sw_var = tk.IntVar()
        ttk.Checkbutton(self.tab_geom, text="Include Self-Weight", variable=self.sw_var).grid(row=20, column=0, columnspan=2, pady=5)

    def _on_span_change(self):
        cumulative_x = 0.0
        support_idx = 0
        
        if len(self.support_rows) > 0:
            self.support_rows[support_idx]["x"].set("0.0")
            support_idx += 1
            
        for var in self.span_vars:
            try:
                val = var.get().strip()
                if val:
                    span_L = float(val)
                    if span_L > 0:
                        cumulative_x += span_L
                        if support_idx < len(self.support_rows):
                            self.support_rows[support_idx]["x"].set(str(round(cumulative_x, 3)))
                            support_idx += 1
            except ValueError:
                pass
                
        for i in range(support_idx, len(self.support_rows)):
            self.support_rows[i]["x"].set("")
            
        self.schedule_preview_update()

    def _build_bound_tab(self):
        ttk.Label(self.tab_bound, text="Supports (Location X, Kx, Ky, Kz, Settlement Y)", style="Header.TLabel").grid(row=0, column=0, columnspan=5, pady=5, sticky="w")
        headers = ["X (m)", "Kx (kN/m)", "Ky (kN/m)", "Kz (rad)", "Settlement"]
        for col, h in enumerate(headers):
            ttk.Label(self.tab_bound, text=h).grid(row=1, column=col, padx=2)
            
        self.support_rows = []
        for i in range(12): 
            row_vars = {}
            for col, key in enumerate(["x", "kx", "ky", "kz", "sy"]):
                var = tk.StringVar()
                var.trace_add("write", lambda *args: self.schedule_preview_update())
                ttk.Entry(self.tab_bound, textvariable=var, width=10).grid(row=2+i, column=col, padx=2, pady=2)
                row_vars[key] = var
            self.support_rows.append(row_vars)

    def _build_loads_tab(self):
        ttk.Label(self.tab_loads, text="Thermal Gradient (°C)", style="Header.TLabel").grid(row=0, column=0, columnspan=2, pady=5, sticky="w")
        self.t_top_var = tk.StringVar()
        self.t_bot_var = tk.StringVar()
        ttk.Label(self.tab_loads, text="Top Temp:").grid(row=1, column=0, sticky="e")
        ttk.Entry(self.tab_loads, textvariable=self.t_top_var, width=10).grid(row=1, column=1, sticky="w")
        ttk.Label(self.tab_loads, text="Bottom Temp:").grid(row=2, column=0, sticky="e")
        ttk.Entry(self.tab_loads, textvariable=self.t_bot_var, width=10).grid(row=2, column=1, sticky="w")
        
        ttk.Label(self.tab_loads, text="Extra Static Loads", style="Header.TLabel").grid(row=3, column=0, columnspan=4, pady=10, sticky="w")
        headers = ["Type", "Mag (kN or kN/m)", "Start X (m)", "End X (m)"]
        for col, h in enumerate(headers):
            ttk.Label(self.tab_loads, text=h).grid(row=4, column=col, padx=2)
            
        self.extra_load_rows = []
        for i in range(4): 
            row_vars = {}
            for col, key in enumerate(["type", "mag", "start", "end"]):
                var = tk.StringVar()
                var.trace_add("write", lambda *args: self.schedule_preview_update())
                
                if key == "type":
                    cb = ttk.Combobox(self.tab_loads, textvariable=var, values=["", "DIST", "POINT"], width=8, state="readonly")
                    cb.grid(row=5+i, column=col, padx=2, pady=2)
                else:
                    ttk.Entry(self.tab_loads, textvariable=var, width=10).grid(row=5+i, column=col, padx=2, pady=2)
                
                row_vars[key] = var
            self.extra_load_rows.append(row_vars)
            
    def _build_train_tab(self):
        ttk.Label(self.tab_train, text="Moving Train Parameters", style="Header.TLabel").grid(row=0, column=0, columnspan=2, pady=5, sticky="w")
        
        self.train_vars = {}
        props = [("Train Name:", "name"), ("Axle Spacing (m):", "spacing"), 
                 ("Acceleration (m/s2):", "acc"), ("Step dx (m):", "dx")]
        for idx, (lbl, key) in enumerate(props):
            ttk.Label(self.tab_train, text=lbl).grid(row=1+idx, column=0, sticky="e")
            var = tk.StringVar()
            ttk.Entry(self.tab_train, textvariable=var, width=15).grid(row=1+idx, column=1, pady=2, sticky="w")
            self.train_vars[key] = var
            
        ttk.Label(self.tab_train, text="Axle Loads (kN) - Max 6 axles:", style="Header.TLabel").grid(row=6, column=0, columnspan=2, pady=10, sticky="w")
        self.axle_vars = []
        for i in range(6):
            ttk.Label(self.tab_train, text=f"Axle {i+1} P:").grid(row=7+i, column=0, sticky="e")
            var = tk.StringVar()
            ttk.Entry(self.tab_train, textvariable=var, width=10).grid(row=7+i, column=1, pady=2, sticky="w")
            self.axle_vars.append(var)

    def _set_default_values(self):
        defaults_spans = ["30.0", "50.0", "30.0", "", "", "", "", "", "", ""]
        for var, val in zip(self.span_vars, defaults_spans): var.set(val)
        
        self.sec_vars["A"].set("1.5")
        self.sec_vars["I"].set("0.125")
        self.sec_vars["h"].set("1.5")
        self.sec_vars["E"].set("36000000.0")  
        self.sec_vars["alpha"].set("1e-05")
        self.sec_vars["gamma"].set("25.0")
        self.sec_vars["epm"].set("1.0")
        self.sw_var.set(1) 

        sup_defs = [
            ["0.0", "2000000.0", "15000000.0", "100000000.0", "0.0"],  
            ["30.0", "600000.0", "8000000.0", "30000000.0", "0.5"],    
            ["80.0", "600000.0", "8000000.0", "30000000.0", "0.0"],    
            ["110.0", "2000000.0", "15000000.0", "100000000.0", "0.0"] 
        ]
        for i, data in enumerate(sup_defs):
            self.support_rows[i]["x"].set(data[0])
            self.support_rows[i]["kx"].set(data[1])
            self.support_rows[i]["ky"].set(data[2])
            self.support_rows[i]["kz"].set(data[3])
            self.support_rows[i]["sy"].set(data[4])

        self.t_top_var.set("45.0")
        self.t_bot_var.set("10.0")
        
        self.extra_load_rows[0]["type"].set("DIST")
        self.extra_load_rows[0]["mag"].set("100.0")
        self.extra_load_rows[0]["start"].set("30.0")
        self.extra_load_rows[0]["end"].set("80.0")

        self.extra_load_rows[1]["type"].set("POINT")
        self.extra_load_rows[1]["mag"].set("100.0")
        self.extra_load_rows[1]["start"].set("50.0")
        self.extra_load_rows[1]["end"].set("0.0") 

        # Default parameters set per request
        self.train_vars["name"].set("Train1")
        self.train_vars["spacing"].set("1.6")
        self.train_vars["acc"].set("1.2")  # Requested default
        self.train_vars["dx"].set("0.25")  # Requested default
        
        axle_defs = ["250.0", "250.0", "250.0", "250.0", "", ""]
        for var, val in zip(self.axle_vars, axle_defs): var.set(val)

    def _compile_config_from_gui(self):
        spans = []
        for var in self.span_vars:
            try:
                if var.get().strip(): spans.append(float(var.get()))
            except ValueError:
                pass 
            
        sec = {
            "A": float(self.sec_vars["A"].get() or 1), "I": float(self.sec_vars["I"].get() or 1),
            "h": float(self.sec_vars["h"].get() or 1), "E": float(self.sec_vars["E"].get() or 1),
            "alpha": float(self.sec_vars["alpha"].get() or 0), "gamma": float(self.sec_vars["gamma"].get() or 0),
            "sw_flag": self.sw_var.get()
        }
        
        springs, settlements = {}, {}
        for row in self.support_rows:
            try:
                if row["x"].get().strip():
                    x = float(row["x"].get())
                    springs[x] = [float(row["kx"].get() or 0), float(row["ky"].get() or 0), float(row["kz"].get() or 0)]
                    settlements[x] = [0.0, float(row["sy"].get() or 0), 0.0]
            except ValueError:
                pass

        extra = []
        for row in self.extra_load_rows:
            if row["type"].get().strip().upper() in ["DIST", "POINT"]:
                extra.append({
                    "type": row["type"].get().strip().upper(),
                    "mag": float(row["mag"].get() or 0),
                    "start_x": float(row["start"].get() or 0),
                    "end_x": float(row["end"].get() or 0)
                })

        axles = []
        for var in self.axle_vars:
            try:
                if var.get().strip(): axles.append(float(var.get()))
            except ValueError: pass

        return {
            "geom": {"span_lengths": spans, "total_L": sum(spans), "elems_per_m": float(self.sec_vars["epm"].get() or 1)},
            "section": sec,
            "extra_loads": extra,
            "springs": springs,
            "settlements": settlements,
            "thermal": {"t_top": float(self.t_top_var.get() or 0), "t_bot": float(self.t_bot_var.get() or 0)},
            "sim": {
                "train_name": self.train_vars["name"].get(), 
                "acc": float(self.train_vars["acc"].get() or 0), 
                "dx": float(self.train_vars["dx"].get() or 1), 
                "loads": axles, 
                "spacing": float(self.train_vars["spacing"].get() or 0)
            }
        }

    def embed_plot(self, fig, parent_frame):
        for widget in parent_frame.winfo_children(): widget.destroy()
        if fig is None: return
        canvas = FigureCanvasTkAgg(fig, master=parent_frame)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, parent_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        plt.close(fig) 

    def schedule_preview_update(self):
        if self.preview_timer:
            self.root.after_cancel(self.preview_timer)
        self.preview_timer = self.root.after(500, self.refresh_model_view_silent)

    def refresh_model_view_silent(self):
        try:
            data_deck = self._compile_config_from_gui()
            fig = PostProcessVisualizer.plot_schematic_model(data_deck)
            self.embed_plot(fig, self.tab_model_view)
        except Exception:
            pass 

    def refresh_model_view(self):
        try:
            data_deck = self._compile_config_from_gui()
            fig = PostProcessVisualizer.plot_schematic_model(data_deck)
            self.embed_plot(fig, self.tab_model_view)
            self.res_notebook.select(self.tab_model_view)
        except Exception as e:
            messagebox.showerror("Error", f"Invalid input parameters:\n{e}")

    def update_progress(self, percent):
        def _update():
            self.progress_var.set(percent)
            self.status_label.config(text=f"Status: Simulating Moving Load... {percent}%")
        self.root.after(0, _update)

    def set_status(self, msg):
        self.root.after(0, lambda: self.status_label.config(text=msg))

    def start_analysis_thread(self):
        try:
            data_deck = self._compile_config_from_gui()
            if data_deck["geom"]["total_L"] == 0:
                messagebox.showwarning("Warning", "Bridge span length cannot be 0.")
                return
        except Exception as e:
            messagebox.showerror("Input Error", f"Invalid input parameters:\n{e}")
            return

        self.btn_run.config(state=tk.DISABLED)
        self.progress_var.set(0)
        self.status_label.config(text="Status: Preparing matrices...")
        
        analysis_thread = threading.Thread(target=self.run_analysis, args=(data_deck,))
        analysis_thread.start()

    def run_analysis(self, data_deck):
        """Runs structural analysis, collects all results, and updates UI tables/plots."""
        try:
            # 1. Start Engine
            orch = TimeStepOrchestrator(data_deck, progress_callback=self.update_progress)
            
            self.set_status("Status: Analyzing...")
            
            # 2. Solve Isolated Load Cases
            res_dead = orch.solve_isolated_dead()
            res_thermal = orch.solve_isolated_thermal()
            res_settlement = orch.solve_isolated_settlement()
            res_live = orch.execute_simulation_march()
            
            # 3. Superposition (Calculation of Envelopes)
            combined_max_f, combined_min_f = {}, {}
            for m in orch.model.members:
                mid = m.id
                f_d = res_dead.get(mid, np.zeros(6)) if res_dead else np.zeros(6)
                f_t = res_thermal.get(mid, np.zeros(6)) if res_thermal else np.zeros(6)
                f_s = res_settlement.get(mid, np.zeros(6)) if res_settlement else np.zeros(6)
                
                f_base = f_d + f_t + f_s
                f_lmax = res_live.max_envelope.get(mid, np.zeros(6)) if res_live else np.zeros(6)
                f_lmin = res_live.min_envelope.get(mid, np.zeros(6)) if res_live else np.zeros(6)
                
                combined_max_f[mid] = f_base + f_lmax
                combined_min_f[mid] = f_base + f_lmin
                
            base_disp = np.zeros(len(orch.model.nodes) * 3)
            if orch.dead_load_displacements is not None: base_disp += orch.dead_load_displacements
            if orch.thermal_displacements is not None: base_disp += orch.thermal_displacements
            if orch.settlement_displacements is not None: base_disp += orch.settlement_displacements
            
            comb_disp_max = base_disp + (orch.live_uy_max if orch.live_uy_max is not None else 0)
            comb_disp_min = base_disp + (orch.live_uy_min if orch.live_uy_min is not None else 0)

            # 4. Prepare UI DataFrames (Rounded to 4 digits)
            # This dictionary maps the DataFrames to their respective UI Frame categories
            results_to_ui = {
                "FORCE_DEAD": (pd.DataFrame.from_dict(res_dead, orient='index', columns=['Fx_i', 'Vy_i', 'Mz_i', 'Fx_j', 'Vy_j', 'Mz_j']).round(4), "DEAD", "forces"),
                "FORCE_THERMAL": (pd.DataFrame.from_dict(res_thermal, orient='index', columns=['Fx_i', 'Vy_i', 'Mz_i', 'Fx_j', 'Vy_j', 'Mz_j']).round(4), "THERMAL", "forces"),
                "FORCE_SETTLE": (pd.DataFrame.from_dict(res_settlement, orient='index', columns=['Fx_i', 'Vy_i', 'Mz_i', 'Fx_j', 'Vy_j', 'Mz_j']).round(4), "SETTLEMENT", "forces"),
                "FORCE_LIVE": (pd.DataFrame(res_live.max_envelope).T.round(4), "LIVE", "forces"),
                "FORCE_COMB": (pd.DataFrame(combined_max_f).T.round(4), "COMBINED", "forces"),
                
                "DISP_DEAD": (create_disp_dataframe(orch.model, orch.dead_load_displacements).round(4), "DEAD", "deflections"),
                "DISP_THERMAL": (create_disp_dataframe(orch.model, orch.thermal_displacements).round(4), "THERMAL", "deflections"),
                "DISP_SETTLE": (create_disp_dataframe(orch.model, orch.settlement_displacements).round(4), "SETTLEMENT", "deflections"),
                "DISP_LIVE": (create_disp_dataframe(orch.model, orch.live_uy_max).round(4), "LIVE", "deflections"),
                "DISP_COMB": (create_disp_dataframe(orch.model, comb_disp_max).round(4), "COMBINED", "deflections")
            }

            # 5. Push data to UI Tables
            for key, (df, cat, sub) in results_to_ui.items():
                if df is not None:
                    self.root.after(0, self._display_df_to_frame, self.frames[cat][sub]["table"], df)

            # 6. Update Plots
            support_locs = sorted(list(data_deck["springs"].keys()))
            
            # Embed all necessary plots
            self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_single_case_forces(orch.model, res_dead, support_locs, "DEAD LOAD"), self.frames["DEAD"]["forces"]["plot"])
            self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_single_case_deflection(orch.model, orch.dead_load_displacements, support_locs, "DEAD LOAD"), self.frames["DEAD"]["deflections"]["plot"])
            
            self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_single_case_forces(orch.model, res_thermal, support_locs, "THERMAL LOAD"), self.frames["THERMAL"]["forces"]["plot"])
            # We add the direction="x" parameter when plotting the displacement graph of the thermal load
            self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_single_case_deflection(
                orch.model, 
                orch.thermal_displacements, 
                support_locs, 
                "THERMAL LOAD", 
                line_color="orange", 
                direction="x"  # THIS PART WAS CORRECTED
            ), self.frames["THERMAL"]["deflections"]["plot"])
            
            self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_single_case_forces(orch.model, res_settlement, support_locs, "SETTLEMENT"), self.frames["SETTLEMENT"]["forces"]["plot"])
            self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_single_case_deflection(orch.model, orch.settlement_displacements, support_locs, "SETTLEMENT", line_color="purple"), self.frames["SETTLEMENT"]["deflections"]["plot"])

            if res_live:
                self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_isolated_live_forces(orch.model, res_live, support_locs), self.frames["LIVE"]["forces"]["plot"])
                self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_isolated_live_deflection(orch.model, orch.live_uy_max, orch.live_uy_min, support_locs), self.frames["LIVE"]["deflections"]["plot"])
            
            self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_combined_forces(orch.model, combined_max_f, combined_min_f, support_locs), self.frames["COMBINED"]["forces"]["plot"])
            self.root.after(0, self.embed_plot, PostProcessVisualizer.plot_combined_deflection(orch, support_locs), self.frames["COMBINED"]["deflections"]["plot"])

            self.root.after(0, lambda: [self.btn_run.config(state=tk.NORMAL), 
                                        self.status_label.config(text="Status: Analysis Complete!")])
            messagebox.showinfo("Success", "Analysis finished successfully.")
        except Exception as e:
            self.root.after(0, lambda: [self.btn_run.config(state=tk.NORMAL), 
                                        self.status_label.config(text="Status: Error encountered.")])
            messagebox.showerror("Execution Error", f"An error occurred:\n{str(e)}")
    
    def reset_all_inputs(self):
        """Clears all input fields to default states."""
        # 1. Clear Spans
        for var in self.span_vars: var.set("")
        
        # 2. Clear Section Properties
        for var in self.sec_vars.values(): var.set("")
        
        # 3. Clear Supports
        for row in self.support_rows:
            for var in row.values(): var.set("")
            
        # 4. Clear Extra Loads
        for row in self.extra_load_rows:
            for var in row.values(): var.set("")
            
        # 5. Clear Train Params
        for var in self.train_vars.values(): var.set("")
        for var in self.axle_vars: var.set("")
        
        self.status_label.config(text="Status: Inputs cleared. Define new geometry.")
    
    def load_defaults_command(self):
        """Resets all fields to preset engineering defaults."""
        if messagebox.askyesno("Load Defaults", "Are you sure you want to reset all parameters to defaults?"):
            self._set_default_values() # Loads current defaults
            self.refresh_model_view_silent() # Update the interface
            self.status_label.config(text="Status: Default parameters loaded.")
    
    def define_new_geometry(self):
        """Resets inputs, clears results/plots, and reactivates the Run button."""
        if messagebox.askyesno("New Design", "Are you sure you want to clear everything and start a new design?"):
            # 1. Clear Inputs
            self.reset_all_inputs()
            # 3. Activate Run Analysis Button
            self.btn_run.config(state=tk.NORMAL)
            self.status_label.config(text="Status: Ready for new design.")

if __name__ == "__main__":
    root = tk.Tk()
    app = BridgeGUI(root)
    root.after(500, app.refresh_model_view)
    root.mainloop()
