import unittest
import numpy as np
import sys
import os
import re
from str_engine_hw3 import Material, Node, Member, StructuralModel

class SuppressOutput:
    def __enter__(self):
        self._original_stdout, self._original_stderr = sys.stdout, sys.stderr
        self._fnull = open(os.devnull, 'w')
        sys.stdout = sys.stderr = self._fnull
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout, sys.stderr = self._original_stdout, self._original_stderr
        self._fnull.close()

class CleanAuditResult(unittest.TextTestResult):
    def startTest(self, test):
        super().startTest(test)
        name = test._testMethodName
        cat = "UNIT" if "unit" in name else "INTERFACE" if "interface" in name else "REGRESSION"
        sys.stdout.write(f"[{cat} TEST] {test._testMethodDoc or name}\n")

    def addSuccess(self, test):
        sys.stdout.write(f"   Status: PASSED\n\n")

    def addFailure(self, test, err):
        sys.stdout.write(f"   Status: FAILED\n   Reason: {str(err[1])}\n\n")

    def addError(self, test, err):
        sys.stdout.write(f"   Status: ERROR\n   Reason: Engine Crash -> {str(err[1]).splitlines()[-1]}\n\n")

class StructuralAuditor(unittest.TestCase):
    def setUp(self):
        self.model = StructuralModel()
        self.active_file = self._get_active_file()

    def _get_active_file(self):
        try:
            with open('str_engine_hw3.py', 'r', encoding='utf-8') as f:
                content = f.read()
                match = re.search(r"parse_file\(['\"](.+?\.txt)['\"]\)", content)
                return match.group(1) if match else None
        except: return None

    # --- UNIT TESTS ---
    def test_unit_01_node_integrity(self):
        """Verify Node object correctly stores fixity and coordinates."""
        n = Node(1, 10.0, 5.0, [1, 1, 0])
        self.assertEqual(n.x, 10.0)
        self.assertEqual(n.fixity, [1, 1, 0])

    def test_unit_02_member_geometry(self):
        """Verify Member.length calculation."""
        n1, n2 = Node(1, 0, 0, [1, 1, 1]), Node(2, 3, 4, [0, 0, 0])
        m = Member(1, n1, n2, Material(1, 1, 1, 1))
        self.assertAlmostEqual(m.length, 5.0)

    def test_unit_03_rotation_matrix(self):
        """Verify Rotation Matrix (R) for a vertical member."""
        n1, n2 = Node(1, 0, 0, [1, 1, 1]), Node(2, 0, 2, [0, 0, 0])
        m = Member(1, n1, n2, Material(1, 1, 1, 1))
        _, R, _ = m.get_matrices()
        self.assertAlmostEqual(R[0,1], 1.0) 

    def test_unit_04_axial_stiffness(self):
        """Verify local axial stiffness component (EA/L)."""
        mat = Material(1, 200, 2, 1) 
        n1, n2 = Node(1, 0, 0, [1, 1, 1]), Node(2, 1, 0, [0, 0, 0])
        m = Member(1, n1, n2, mat)
        k_loc, _, _ = m.get_matrices()
        self.assertAlmostEqual(k_loc[0,0], 400.0)

    def test_unit_05_dof_mapping(self):
        """Verify DOF indices correctly map free/fixed DOFs."""
        self.model.nodes = {1: Node(1,0,0,[1,1,1]), 2: Node(2,1,0,[0,0,0])}
        with SuppressOutput(): self.model.solve()
        self.assertEqual(self.model.nodes[1].dof_indices, [-1, -1, -1])

    # --- INTERFACE TESTS ---
    def test_interface_06_load_mapping(self):
        """Verify nodal loads are mapped to global results."""
        if self.active_file:
            with SuppressOutput(): 
                self.model.parse_file(self.active_file)
                self.model.solve()
            self.assertIsNotNone(self.model.U_res)
        else: self.skipTest("No file.")

    def test_interface_07_kg_symmetry(self):
        """Verify Global Stiffness Matrix (kg) from members is symmetric."""
        m = Member(1, Node(1,0,0,[0,0,0]), Node(2,1,1,[0,0,0]), Material(1,2e8,0.1,0.1))
        _, _, kg = m.get_matrices()
        self.assertTrue(np.allclose(kg, kg.T))

    def test_interface_08_vector_dimensions(self):
        """Verify Result vector size matches active DOFs in file."""
        if self.active_file:
            with SuppressOutput():
                self.model.parse_file(self.active_file)
                self.model.solve()
            active_dofs = sum(1 for n in self.model.nodes.values() for f in n.fixity if f == 0)
            if self.model.U_res is not None:
                self.assertEqual(len(self.model.U_res), active_dofs)
        else: self.skipTest("No file.")

    # --- REGRESSION TESTS  ---
    def test_reg_09_instability_audit(self):
        """Audit Instability: Detect mechanisms or singular matrices in the model."""
        if self.active_file:
            with SuppressOutput():
                self.model.parse_file(self.active_file)
                status = self.model.solve()
            
            # If the status false:
            if status is False:
                self.fail("INSTABILITY: Global matrix is singular or solver failed.")
            
            # If the mechanism exist:
            if self.model.U_res is not None and len(self.model.U_res) > 0:
                if np.any(np.abs(self.model.U_res) > 1e10):
                    self.fail("INSTABILITY: Large displacements detected (Rigid Body Motion).")
        else: self.skipTest("No file.")

    def test_reg_10_connectivity_audit(self):
        """Audit Connectivity: Detect number of independent structures."""
        if self.active_file:
            with SuppressOutput(): self.model.parse_file(self.active_file)
            islands = self.model.analyze_connectivity()
            if len(islands) > 1:
                self.fail(f"MULTI-STRUCTURE: {len(islands)} independent structures detected.")
        else: self.skipTest("No file.")

    def test_reg_11_flying_object_audit(self):
        """Audit Flying Objects: Identify islands with zero boundary conditions."""
        if self.active_file:
            with SuppressOutput(): self.model.parse_file(self.active_file)
            islands = self.model.analyze_connectivity()
            for i, island in enumerate(islands):
                # At least one node of a seperate strcuture should be restrained.
                has_support = any(any(self.model.nodes[nid].fixity) for nid in island)
                if not has_support:
                    self.fail(f"FLYING OBJECT: Structure Island #{i+1} has no supports!")
        else: self.skipTest("No file.")

    def test_reg_12_equilibrium(self):
        """Verify solve produces results for loaded systems."""
        if self.active_file:
            with SuppressOutput(): 
                self.model.parse_file(self.active_file)
                self.model.solve()
            if self.model.U_res is not None:
                self.assertTrue(np.any(np.abs(self.model.U_res) > 1e-12))
        else: self.skipTest("No file.")

    def test_reg_13_island_stability_audit(self):
        """Audit Island Stability: Independent check for each structural part."""
        if self.active_file:
            islands = self.model.analyze_connectivity()
            for i, island in enumerate(islands):
                # Each independent structure should have at least one horizontal restrain.
                # If the independent structure is stable: no failure.
                if not any(self.model.nodes[nid].fixity[0] for nid in island):
                    self.fail(f"LOCAL INSTABILITY: Island #{i+1} lacks horizontal restraint (trX=0).")
        else: self.skipTest("No file.")

if __name__ == "__main__":
    print("\n" + "="*65)
    print(f"{'CE 4011 FINAL AUDIT REPORT':^65}")
    print("="*65 + "\n")
    suite = unittest.TestLoader().loadTestsFromTestCase(StructuralAuditor)
    runner = unittest.TextTestRunner(resultclass=CleanAuditResult, verbosity=0, stream=sys.stdout)
    runner.run(suite)
  
    
    suite = unittest.TestLoader().loadTestsFromTestCase(StructuralAuditor)
    #Run
    result = unittest.TestResult()
    suite.run(result)
    
    # Failure Count
    total_tests = result.testsRun
    failed_count = len(result.failures)
    error_count = len(result.errors)
    passed_count = total_tests - failed_count - error_count
    
    print("="*65)
    print(f"AUDIT SUMMARY:")
    print(f"  - Total Tests Run: {total_tests}")
    print(f"  - PASSED: {passed_count}")
    print(f"  - FAILED: {failed_count}")
    print(f"  - ERRORS: {error_count}")
    print("-" * 65)

   
    if failed_count > 0 or error_count > 0:
        print(f"!!! ATTENTION: {failed_count + error_count} AUDIT ISSUE(S) DETECTED !!!")
        print("WARNING: THE MODEL MAY BE UNSTABLE OR IMPROPERLY DEFINED.")
        print("CALCULATIONS MAY BE INCORRECT OR PHYSICALLY INVALID.")
        print("DO NOT USE THESE RESULTS FOR FINAL DESIGN.")
    else:
        print("AUDIT COMPLETE: All structural checks passed.")
        print("The model appears stable and properly connected.")
    
    print("="*65 + "\n")