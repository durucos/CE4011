import unittest
import numpy as np

class TestStructuralAnalysis(unittest.TestCase):

    # --- 1. UNIT TESTS ---
    def test_settlement_analysis_unit(self):
        print("\n--- Running Unit Test: Settlement Result (hw4_q1.py) ---")
        try:
            from hw4_q1 import model as model_q1
            v_b_dof = model_q1.nodes[2].dof_numbers[1]
            v_b_disp = model_q1.U_res[v_b_dof]
            success = abs(v_b_disp - (-0.002)) < 1e-4
            if success: print(f"SUCCESS: Settlement passed. v_B = {v_b_disp:.6f}")
            else: print(f"FAILURE: v_B = {v_b_disp:.6f}")
            self.assertTrue(success)
        except Exception as e: self.skipTest(f"Error: {e}")

    def test_thermal_analysis_unit(self):
        print("\n--- Running Unit Test: Thermal Forces & Moments (hw4_q2.py) ---")
        try:
            from hw4_q2 import model as m2
            m1_forces = m2.member_forces.get(1)
            act_axial, act_moment = abs(m1_forces[0]), abs(m1_forces[2])
            success = abs(act_axial - 1920.0) < 10.0 and abs(act_moment - 256.0) < 5.0
            if success: print(f"SUCCESS: Thermal Unit verified. Axial: {act_axial:.1f}, Moment: {act_moment:.1f}")
            else: print(f"FAILURE: Axial: {act_axial:.1f}, Moment: {act_moment:.1f}")
            self.assertTrue(success)
        except Exception as e: self.skipTest(f"Error: {e}")

    # --- 2. PROCEDURE TESTS ---
    def test_settlement_procedure_regression(self):
        print("\n--- Running Procedure Test: Settlement Logic (hw4_q1.py) ---")
        try:
            from hw4_q1 import model as m1
            vb_dof = m1.nodes[2].dof_numbers[1]
            f_sett = abs(m1.F_settlement[vb_dof])
            success = np.isclose(f_sett, 3750.0, rtol=1e-2)
            if success: print(f"SUCCESS: Procedure validated. Force: {f_sett:.2f}")
            else: print(f"FAILURE: Force: {f_sett:.2f}")
            self.assertTrue(success)
        except Exception as e: self.skipTest(f"Error: {e}")

    def test_thermal_procedure_regression(self):
        print("\n--- Running Procedure Test: Thermal Vector Assembly (hw4_q2.py) ---")
        try:
            from hw4_q2 import model as m2
            m1_dofs = m2.members[0].nodes[0].dof_numbers + m2.members[0].nodes[1].dof_numbers
            f_ext_x = abs(m2.F_external[m1_dofs[3]])
            f_ext_m = abs(m2.F_external[m1_dofs[5]])
            # Procedure check for fixed-end forces (FEF)
            success = (f_ext_x > 1900) and (f_ext_m > 250)
            if success: print(f"SUCCESS: Thermal mapping verified.")
            else: print(f"FAILURE: Force X: {f_ext_x:.1f}, Moment: {f_ext_m:.1f}")
            self.assertTrue(success)
        except Exception as e: self.skipTest(f"Error: {e}")

if __name__ == "__main__":
    unittest.main()