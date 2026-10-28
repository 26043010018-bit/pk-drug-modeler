import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class PharmacokineticModeler:
    """
    Bioengineering tool to simulate and analyze one-compartment 
    pharmacokinetic (PK) drug concentration-time profiles following IV bolus administration.
    """

    def __init__(self, dose_mg: float, volume_of_distribution_L: float, elimination_rate_constant: float):
        self.dose = float(dose_mg)
        self.vd = float(volume_of_distribution_L)
        self.ke = float(elimination_rate_constant)
        self.c0 = self.dose / self.vd  # Initial peak concentration (C0)

    def calculate_parameters(self) -> dict:
        """Derives fundamental pharmacokinetic parameters from the compartmental model."""
        half_life = np.log(2) / self.ke
        clearance = self.ke * self.vd
        auc = self.dose / clearance  # Area Under the Curve (0 to infinity)

        return {
            'initial_concentration_mg_L': self.c0,
            'half_life_hr': half_life,
            'clearance_L_hr': clearance,
            'auc_mg_hr_L': auc
        }

    def simulate_profile(self, time_hours: np.ndarray) -> np.ndarray:
        """Calculates plasma drug concentration at given time points using C(t) = C0 * e^(-ke * t)."""
        return self.c0 * np.exp(-self.ke * time_hours)

    def generate_pk_plot(self, max_time: float = 24.0, output_filename: str = 'pk_profile.png'):
        """Generates and saves a pharmacokinetic plasma concentration vs. time plot."""
        time = np.linspace(0, max_time, 200)
        concentration = self.simulate_profile(time)
        params = self.calculate_parameters()

        plt.figure(figsize=(8, 5))
        
        # Plot concentration curve
        plt.plot(time, concentration, color='#1b9e77', linewidth=2, 
                 label=f'Plasma Concentration\n(t1/2 = {params["half_life_hr"]:.2f} hrs)')
        
        # Highlight half-life point on curve
        hl = params['half_life_hr']
        c_hl = self.c0 / 2.0
        plt.scatter([hl], [c_hl], color='#d95f02', s=70, zorder=5, 
                    label=f'Half-life Marker ({hl:.1f}h, {c_hl:.1f} mg/L)')
        
        plt.title('One-Compartment Pharmacokinetic (PK) Profile', fontsize=12, fontweight='bold')
        plt.xlabel('Time Post-Dose (Hours)')
        plt.ylabel('Plasma Concentration (mg/L)')
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.legend(loc='upper right')
        
        plt.tight_layout()
        plt.savefig(output_filename, dpi=300)
        plt.close()
        print(f"[+] PK plot successfully saved to {output_filename}")


if __name__ == '__main__':
    # Simulation parameters: 
    # Dose = 500 mg, Volume of Distribution (Vd) = 30 Liters, Elimination rate constant (ke) = 0.15 hr^-1
    dose = 500.0
    vd = 30.0
    ke = 0.15

    print("=== PHARMACOKINETIC MODELER EXECUTING ===")
    pk_tool = PharmacokineticModeler(dose_mg=dose, volume_of_distribution_L=vd, elimination_rate_constant=ke)

    # Compute parameters
    metrics = pk_tool.calculate_parameters()

    print("\n[Derived Pharmacokinetic Metrics]")
    print(f"  Initial Peak Concentration (C0) : {metrics['initial_concentration_mg_L']:.2f} mg/L")
    print(f"  Biological Half-Life (t1/2)     : {metrics['half_life_hr']:.2f} hours")
    print(f"  Total Clearance (Cl)            : {metrics['clearance_L_hr']:.2f} L/hr")
    print(f"  Area Under Curve (AUC)          : {metrics['auc_mg_hr_L']:.2f} mg·hr/L")

    # Generate simulation plot
    pk_tool.generate_pk_plot(max_time=24.0)
    print("\n[✓] Pharmacokinetic simulation complete.")