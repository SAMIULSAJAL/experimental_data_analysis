import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from scipy import optimize, special
from scipy.stats import linregress 

# --- Configuration ---
BASE_LAB_DIR = "/Users/ssajal/Desktop/academic/Atomic_and_molecular_lab/M_1.8_Mikrowellen_Spektroskopie/our_lab/"

# Section 1 Paths
SECTION1_DATA_FOLDER = os.path.join(BASE_LAB_DIR, "data", "20250610")
SECTION1_OUTPUT_FOLDER = os.path.join(BASE_LAB_DIR, "output", "section_1_mode_plots")

# Section 2 Paths
SECTION2_DATA_FOLDER = os.path.join(BASE_LAB_DIR, "data", "20250610", "2f_fm_transition")
SECTION2_OUTPUT_FOLDER = os.path.join(BASE_LAB_DIR, "output", "section_2_analysis_results")

# Section 3 Paths
SECTION3_DATA_FOLDER = os.path.join(BASE_LAB_DIR, "data", "20250610", "2f_fm_integration time(6,1)")
SECTION3_OUTPUT_FOLDER = os.path.join(BASE_LAB_DIR, "output", "section_3")

# Section 4 Paths
SECTION4_DATA_FOLDER = os.path.join(BASE_LAB_DIR, "data", "20250610", "2f_fm_amplitude(3,3)")
SECTION4_OUTPUT_FOLDER = os.path.join(BASE_LAB_DIR, "output", "section_4")
SECTION4_TABLES_FILE = os.path.join(SECTION4_OUTPUT_FOLDER, 'overleaf_section4_tables.txt')

# Section 5 Paths - NEW
SECTION5_DATA_FOLDER = os.path.join(BASE_LAB_DIR, "data", "20250610", "2f_fm_pressure(3,3)")
SECTION5_OUTPUT_FOLDER = os.path.join(BASE_LAB_DIR, "output", "section_5")
SECTION5_TABLES_FILE = os.path.join(SECTION5_OUTPUT_FOLDER, 'overleaf_section5_tables.txt')


# Ensure output directories exist
os.makedirs(SECTION1_OUTPUT_FOLDER, exist_ok=True)
os.makedirs(SECTION2_OUTPUT_FOLDER, exist_ok=True)
os.makedirs(SECTION3_OUTPUT_FOLDER, exist_ok=True)
os.makedirs(SECTION4_OUTPUT_FOLDER, exist_ok=True)
os.makedirs(SECTION5_OUTPUT_FOLDER, exist_ok=True) # New output directory

# Theoretical values for comparison
NU_THEORY_LOOKUP = {
    (3,3): 23870.11,
    (4,3): 22688.24,
    (4,4): 24139.39,
    (5,3): 21285.3,
    (5,5): 24532.94, # Theoretical value for (5,5) still listed for reference
    (6,6): 25056.04,
    (7,6): 26000.00 # Placeholder: Update with exact theoretical value from tutor!
}
NU0_THEORY = 23786
A_THEORY = 151.5
B_THEORY = 59.9

# --- Helper Functions ---

def load_data(filepath, use_columns=(0, 1)):
    """Loads frequency and signal from a specified two-column, whitespace-delimited file."""
    try:
        data = np.genfromtxt(filepath, comments='#', delimiter=None, usecols=use_columns, encoding='utf-8')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] < 2:
            raise ValueError("File does not contain at least two data columns.")
        return data.astype(float)
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def load_pressure_data_from_header(filepath):
    """
    Loads frequency and signal from a file, and extracts pressure voltages from header comments.
    Returns (xs, ys, pressure_mV). pressure_mV is the average if 'before' and 'after' are present.
    """
    pressure_before_mV = None
    pressure_after_mV = None
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                if line.strip().startswith('#'):
                    if 'Pressure before:' in line:
                        pressure_before_mV = float(line.split(':')[1].strip().replace('mV', ''))
                    elif 'Pressure after:' in line:
                        pressure_after_mV = float(line.split(':')[1].strip().replace('mV', ''))
                else:
                    break # Stop reading comments once data lines start
        
        # Load the actual numerical data
        data = np.genfromtxt(filepath, comments='#', delimiter=None, usecols=(0, 1), encoding='utf-8')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] < 2:
            raise ValueError("File does not contain at least two data columns.")
        
        # Determine the pressure_mV to return
        if pressure_before_mV is not None and pressure_after_mV is not None:
            avg_pressure_mV = (pressure_before_mV + pressure_after_mV) / 2
        elif pressure_before_mV is not None:
            avg_pressure_mV = pressure_before_mV
        elif pressure_after_mV is not None: # Less likely but for completeness
            avg_pressure_mV = pressure_after_mV
        else:
            avg_pressure_mV = np.nan # No pressure info found in header

        return data[:, 0], data[:, 1], avg_pressure_mV

    except Exception as e:
        print(f"Error loading {filepath} or parsing header: {e}")
        return None, None, np.nan


def plot_data(xs, ys, title, xlabel, ylabel, save_path, fit_xs=None, fit_ys=None, legend_labels=None, show_legend=True):
    """Plots data, optionally with a fit, and saves the figure."""
    plt.figure(figsize=(10, 6))
    plt.plot(xs, ys, color='black', label=legend_labels[0] if legend_labels else 'Experimental Data', linewidth=0.8)
    if fit_xs is not None and fit_ys is not None:
        plt.plot(fit_xs, fit_ys, color='red', label=legend_labels[1] if legend_labels else 'Fit', linewidth=1.5)
    
    # Increased font sizes for better visibility in figures
    plt.title(title, fontsize=24)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    if show_legend:
        plt.legend(fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

# --- Tutor's Fitting Functions ---

def Voigt(derivative, x, x0, amp, fwhm_gauss, fwhm_lorentz):
    sigma = fwhm_gauss / (2 * np.sqrt(2 * np.log(2)))
    gamma = fwhm_lorentz / 2

    if np.isclose(sigma, 0) and np.isclose(gamma, 0):
        if derivative == 0:
            return np.array([amp if np.isclose(val, x0, atol=1e-9) else 0 for val in x])
        else:
            return np.zeros_like(x)

    z = (x - x0 + 1j * gamma) / (sigma * np.sqrt(2))
    wz = special.wofz(z)
    
    z0 = (1j * gamma) / (sigma * np.sqrt(2))
    w0 = special.wofz(z0)

    if derivative == 0:
        tmp = lambda x_val, x0_val, wz_val, sigma_val, gamma_val: np.real(wz_val) / (sigma_val * np.sqrt(2 * np.pi))            
    elif derivative == 1:
        tmp = lambda x_val, x0_val, wz_val, sigma_val, gamma_val: 1 / (sigma_val**3 * np.sqrt(2 * np.pi)) * (gamma_val * np.imag(wz_val) - (x_val - x0_val) * np.real(wz_val))
    elif derivative == 2:
        tmp = lambda x_val, x0_val, wz_val, sigma_val, gamma_val: 1 / (sigma_val**5 * np.sqrt(2 * np.pi)) * (gamma_val * (2 * (x_val - x0_val) * np.imag(wz_val) - sigma_val * np.sqrt(2 / np.pi)) + (gamma_val**2 + sigma_val**2 - (x_val - x0_val)**2) * np.real(wz_val))
    else:
        raise NotImplementedError('Only the zeroth, first, and second derivatives of a Voigt profile are implemented.')
    
    ys = tmp(x, x0, wz, sigma, gamma)
    
    ymax_ref = tmp(0, 0, w0, sigma, gamma)
    if np.isclose(ymax_ref, 0):
        ymax_ref = 1e-9 if amp >= 0 else -1e-9
    
    ys *= amp / ymax_ref
    return ys

def lineshape(xs, x0, y0, gauss, lorentz, a, b, c):
    """Combines 2nd derivative of Voigt with a quadratic baseline."""
    ys_voigt_deriv2 = Voigt(2, xs, x0, y0, gauss, lorentz)
    ys_baseline = a + xs * b + xs ** 2 * c
    return ys_voigt_deriv2 + ys_baseline

# --- Linear Model Function for the Ansatz (Equation 42) ---
def inversion_frequency_ansatz_linear(X_terms, nu0, a_const, b_const):
    """Linear model for ammonia inversion frequencies: nu = nu0 - a(J(J+1) - K^2) + bK^2."""
    J_K_term = X_terms[:, 0]
    K2_term = X_terms[:, 1]
    return nu0 - a_const * J_K_term + b_const * K2_term

def convert_voltage_to_pressure_mbar(voltage_mV):
    """Converts voltage in mV to pressure in mbar using the tutor's formula."""
    voltage_Volt = float(voltage_mV) / 1000.0 # Convert mV to V
    pressure_bar = 10 ** (voltage_Volt - 5.5)

    # Convert to mbar based on the likely intended scale for these measurements.
    # The tutor's script converts to 'ubar' if power_of_ten < 0, but since the reference
    # uses mbar for similar low values, we explicitly convert to mbar here for consistency.
    pressure_mbar = pressure_bar * 1000 

    return pressure_mbar


# --- Analysis Sections ---

def run_section1_plots():
    print("\n--- Section 1: Plotting 1F-FM, 2F-FM, and AM modes ---")
    
    file_1ffm_path = os.path.join(SECTION1_DATA_FOLDER, '1ffm.csv')
    file_2ffm_path = os.path.join(SECTION1_DATA_FOLDER, '2ffm.csv')
    file_am_path = os.path.join(SECTION1_DATA_FOLDER, 'am.csv')

    data_1ffm = load_data(file_1ffm_path)
    data_2ffm = load_data(file_2ffm_path)
    data_am = load_data(file_am_path)

    if data_1ffm is None or data_2ffm is None or data_am is None:
        print("Skipping Section 1 plots due to data loading errors.")
        return

    plot_data(data_2ffm[:, 0], data_2ffm[:, 1], '2F FM', # Title changed
              'Frequency (MHz)', 'Intensity', # Y-label changed
              os.path.join(SECTION1_OUTPUT_FOLDER, '3_3_NH3_2FFM.png'),
              legend_labels=['2F FM Mode'], show_legend=False) # No legend for single data series

    plot_data(data_1ffm[:, 0], data_1ffm[:, 1], '1F FM', # Title changed
              'Frequency (MHz)', 'Intensity', # Y-label changed
              os.path.join(SECTION1_OUTPUT_FOLDER, '3_3_NH3_1FFM.png'),
              legend_labels=['1F FM Mode'], show_legend=False) # No legend for single data series

    plot_data(data_am[:, 0], data_am[:, 1], 'AM', # Title changed
              'Frequency (MHz)', 'Intensity', # Y-label changed
              os.path.join(SECTION1_OUTPUT_FOLDER, '3_3_NH3_FAM.png'),
              legend_labels=['AM Mode'], show_legend=False) # No legend for single data series

    print(f"Section 1 plots saved to: {SECTION1_OUTPUT_FOLDER}")

def run_section2_analysis():
    print("\n--- Section 2: Peak Fitting for Multiple Transitions & Constant Determination ---")
    
    TRANSITION_FILES_INFO = [
        {'filename': '3,3.csv', 'J': 3, 'K': 3},
        {'filename': '4,3_500.csv', 'J': 4, 'K': 3},
        {'filename': '4,4_500.csv', 'J': 4, 'K': 4},
        {'filename': '5,3_500.csv', 'J': 5, 'K': 3},
        {'filename': '7,6_500.csv', 'J': 7, 'K': 6},
        {'filename': '6,6_500.csv', 'J': 6, 'K': 6},
    ]

    measured_peak_centers = []

    print("\n--- Fitting individual transition peaks ---")
    for transition_info in TRANSITION_FILES_INFO:
        fname = transition_info['filename']
        J, K = transition_info['J'], transition_info['K']
        filepath = os.path.join(SECTION2_DATA_FOLDER, fname)

        if not os.path.exists(filepath):
            print(f"Skipping {fname} (J={J}, K={K}): File not found at {filepath}")
            continue

        data = load_data(filepath)
        if data is None:
            print(f"Skipping fitting for {fname} due to loading error.")
            continue

        xs, ys = data[:, 0], data[:, 1]

        xmin, xmax = xs.min(), xs.max()
        x0_guess_initial = (xmax + xmin) / 2
        ymin, ymax = ys.min(), ys.max()
        yptp = ymax - ymin

        xs_weight_factor = 4
        ys_weighted = (ys - ymin) * np.exp(- np.abs(np.abs(xs - x0_guess_initial) / (xmax - xmin)) * xs_weight_factor)
        x0_guess = xs[np.argmax(ys_weighted)]

        # Initial guesses for parameters: x0, y0, gauss, lorentz, a, b, c
        p0 = [x0_guess, yptp, 1, 1, 0, 0, 0]

        bounds = ([xmin, -np.inf, 0.01, 0.01, -np.inf, -np.inf, -np.inf],
                  [xmax, np.inf, xmax - xmin, xmax - xmin, np.inf, np.inf, np.inf])

        try:
            popt, pcov = optimize.curve_fit(lineshape, xs, ys, p0=p0, bounds=bounds, maxfev=50000)
            perr = np.sqrt(np.diag(pcov))

            nu_exp = popt[0]
            nu_exp_error = perr[0]

            y_fit = lineshape(xs, *popt)
            chi_squared = np.sum(((ys - y_fit) / np.std(ys))**2)

            # Store all necessary parameters for the LaTeX tables
            measured_peak_centers.append({
                'J': J,
                'K': K,
                'nu_exp': nu_exp,
                'nu_exp_error': nu_exp_error,
                'gauss_fwhm': popt[2],
                'lorentz_fwhm': popt[3],
                'amplitude_y0': popt[1],
                'chi_squared': chi_squared
            })

            print(f"  Fitted ({J},{K}): nu_exp = {nu_exp:.3f} +/- {nu_exp_error:.3f} MHz")
            
            fit_xs = np.linspace(xs.min(), xs.max(), 1000)
            fit_ys = lineshape(fit_xs, *popt)
            plot_data(xs, ys, f'{J},{K}', # Title changed
                      'Frequency (MHz)', 'Signal Intensity (V)',
                      os.path.join(SECTION2_OUTPUT_FOLDER, f'Fit_Transition_{J}_{K}.png'),
                      fit_xs, fit_ys, legend_labels=['Data', 'Fit'], show_legend=True) # Legend keys changed

        except Exception as e:
            print(f"  ERROR fitting ({J},{K}) transition: {e}. Skipping.")


    # --- Linear Regression for Constants ---
    if not measured_peak_centers:
        print("\nNo successful peak fits obtained for linear regression.")
        return

    measured_peak_centers.sort(key=lambda x: (x['J'], x['K']))

    nu_exp_array = np.array([mpc['nu_exp'] for mpc in measured_peak_centers])
    nu_exp_errors_array = np.array([mpc['nu_exp_error'] for mpc in measured_peak_centers])

    X_reg_terms = np.array([[(mpc['J'] * (mpc['J'] + 1)) - mpc['K']**2, mpc['K']**2]
                            for mpc in measured_peak_centers])

    print("\n--- Determining Constants (nu0, a_mol, b_mol) using Weighted Least Squares ---")
    try:
        popt_ansatz, pcov_ansatz = optimize.curve_fit(inversion_frequency_ansatz_linear,
                                                       X_reg_terms, nu_exp_array,
                                                       sigma=nu_exp_errors_array, absolute_sigma=True,
                                                       p0=[NU0_THEORY, A_THEORY, B_THEORY])

        nu0_exp, a_exp, b_exp = popt_ansatz
        perr_ansatz = np.sqrt(np.diag(pcov_ansatz))
        nu0_error, a_error, b_error = perr_ansatz

        print(f"Experimental nu0: {nu0_exp:.2f} +/- {nu0_error:.2f} MHz")
        print(f"Experimental a: {a_exp:.2f} +/- {a_error:.2f} MHz")
        print(f"Experimental b: {b_exp:.2f} +/- {b_error:.2f} MHz")

        # --- GENERATE OVERLEAF LATEX TABLES ---
        output_overleaf_file = os.path.join(SECTION2_OUTPUT_FOLDER, 'overleaf.txt')
        with open(output_overleaf_file, 'w') as f:
            f.write("% --- LaTeX Code for Tables (Requires \\usepackage{amsmath} and optionally \\usepackage{multirow} for complex layouts) ---\n\n")

            # Table 1: Individual Transition Frequencies
            f.write("\\begin{table}[h]\n")
            f.write("    \\centering\n")
            f.write("    \\begin{tabular}{|c|c|c|}\n")
            f.write("        \\hline\n")
            f.write("        Transition & $\\nu_{\\text{exp}}$ [MHz] & $\\nu_{\\text{theory}}$ [MHz] \\\\\n") # Headers changed
            f.write("        \\hline\n")
            for mpc in measured_peak_centers:
                J, K = mpc['J'], mpc['K']
                nu_exp_val = mpc['nu_exp']
                nu_exp_err = mpc['nu_exp_error']
                nu_theory_val = NU_THEORY_LOOKUP.get((J, K), "N/A")
                
                nu_theory_str = f"{nu_theory_val:.2f}" if isinstance(nu_theory_val, (int, float)) else str(nu_theory_val)

                f.write(f"        {J},{K} & {nu_exp_val:.3f}$\\pm${nu_exp_err:.3f} & {nu_theory_str} \\\\\n") # Data format changed
            f.write("        \\hline\n")
            f.write("    \\end{tabular}\n")
            f.write("    \\caption{Inversion frequencies for 6 different transitions}\n") # Caption changed and moved
            f.write("    \\label{tab:individual_frequencies}\n")
            f.write("\\end{table}\n\n")

            # Table 2: Values obtained for the constants nu0, a and b
            f.write("\\begin{table}[h]\n")
            f.write("    \\centering\n")
            f.write("    \\begin{tabular}{|c|c|c|}\n")
            f.write("        \\hline\n")
            f.write("        Constant & Experimental value & Theoretical value \\\\\n") # Headers changed
            f.write("        \\hline\n")
            f.write(f"        $\\nu_0$ & {nu0_exp:.2f}$\\pm${nu0_error:.2f} & {NU0_THEORY:.2f} \\\\\n") # Data format changed
            f.write(f"        a & {a_exp:.2f}$\\pm${a_error:.2f} & {A_THEORY:.2f} \\\\\n") # Data format changed
            f.write(f"        b & {b_exp:.2f}$\\pm${b_error:.2f} & {B_THEORY:.2f} \\\\\n") # Data format changed
            f.write("        \\hline\n")
            f.write("    \\end{tabular}\n")
            f.write("    \\caption{Values obtained for the constants $\\nu_0$, a and b}\n") # Caption changed and moved
            f.write("    \\label{tab:molecular_constants}\n")
            f.write("\\end{table}\n")
            print(f"\nLaTeX tables saved to '{output_overleaf_file}' successfully.")

    except Exception as e:
        print(f"Error during linear regression or LaTeX generation: {e}")

    print(f"Section 2 analysis plots and results saved to: {SECTION2_OUTPUT_FOLDER}")

def run_section3_analysis():
    print("\n--- Section 3: Variation of Line Intensity with Integration Time (Transition 6,1) ---")

    integration_times = {
        '100ms': '100ms.csv',
        '200ms': '200ms.csv',
        '500ms': '500ms.csv'
    }
    
    J_val, K_val = 6, 1 # Specific transition for this section

    for label, filename in integration_times.items():
        filepath = os.path.join(SECTION3_DATA_FOLDER, filename)
        data = load_data(filepath)

        if data is None:
            print(f"Skipping plot for {label} due to data loading error.")
            continue

        xs, ys = data[:, 0], data[:, 1]

        # --- Fitting the data ---
        xmin, xmax = xs.min(), xs.max()
        x0_guess_initial = (xmax + xmin) / 2
        ymin, ymax = ys.min(), ys.max()
        yptp = ymax - ymin

        xs_weight_factor = 4
        # Apply a weighted guess to focus on the peak region
        ys_weighted = (ys - ymin) * np.exp(- np.abs(np.abs(xs - x0_guess_initial) / (xmax - xmin)) * xs_weight_factor)
        x0_guess = xs[np.argmax(ys_weighted)]

        # Initial guesses for parameters: x0, y0, gauss, lorentz, a, b, c
        p0 = [x0_guess, yptp, 1, 1, 0, 0, 0]

        bounds = ([xmin, -np.inf, 0.01, 0.01, -np.inf, -np.inf, -np.inf],
                  [xmax, np.inf, xmax - xmin, xmax - xmin, np.inf, np.inf, np.inf])

        try:
            popt, pcov = optimize.curve_fit(lineshape, xs, ys, p0=p0, bounds=bounds, maxfev=50000)
            fit_xs = np.linspace(xs.min(), xs.max(), 1000)
            fit_ys = lineshape(fit_xs, *popt)

            # Plotting the data with the fit
            plot_data(xs, ys, 
                      title=f'{J_val},{K_val} - {label}', # Title as "6,1 - 100ms"
                      xlabel='Frequency (MHz)', 
                      ylabel='Intensity', # Y-label changed
                      save_path=os.path.join(SECTION3_OUTPUT_FOLDER, f'fit_{J_val}_{K_val}_{label}.png'),
                      fit_xs=fit_xs, 
                      fit_ys=fit_ys, 
                      legend_labels=['Data', 'Fit'], # Legend keys
                      show_legend=True)

            print(f"  Plot for ({J_val},{K_val}) at {label} with fit saved.")

        except Exception as e:
            print(f"  ERROR fitting and plotting ({J_val},{K_val}) at {label}: {e}. Skipping.")

    print(f"Section 3 plots saved to: {SECTION3_OUTPUT_FOLDER}")

def run_section4_analysis():
    print("\n--- Section 4: Dependence of Line Shape and Intensity on FM Amplitude (Transition 3,3) ---")

    # Mapping of file names (e.g., '100') to FM Amplitude in MHz
    # We will use the first 4 as per the reference, assuming 100 = 0.1 MHz, etc.
    fm_amplitudes_settings = {
        '100': 0.1,
        '200': 0.2,
        '300': 0.3,
        '400': 0.4
    }
    
    J_val, K_val = 3, 3 # Specific transition for this section

    results_fm_amplitude = [] # To store data for tables

    for setting_str, fm_amp_mhz in fm_amplitudes_settings.items():
        filename = f'{setting_str}.csv'
        filepath = os.path.join(SECTION4_DATA_FOLDER, filename)
        
        # Load data, specifically using columns 0 (frequency) and 1 (signal)
        data = load_data(filepath, use_columns=(0, 1)) 

        if data is None:
            print(f"Skipping analysis for FM Amplitude {fm_amp_mhz} MHz due to data loading error.")
            continue

        xs, ys = data[:, 0], data[:, 1]

        # --- Fitting the data (similar to Section 2/3) ---
        xmin, xmax = xs.min(), xs.max()
        x0_guess_initial = (xmax + xmin) / 2
        ymin, ymax = ys.min(), ys.max()
        yptp = ymax - ymin

        xs_weight_factor = 4
        # Apply a weighted guess to focus on the peak region
        ys_weighted = (ys - ymin) * np.exp(- np.abs(np.abs(xs - x0_guess_initial) / (xmax - xmin)) * xs_weight_factor)
        x0_guess = xs[np.argmax(ys_weighted)] # x0_guess is now correctly defined before p0

        # Increased the range for FWHM guesses. This might help if initial guesses were too restrictive.
        max_fwhm_guess = xmax - xmin 
        
        # Initial guesses for parameters: x0, y0, gauss, lorentz, a, b, c
        p0 = [x0_guess, yptp, 0.5, 0.5, 0, 0, 0] # Slightly adjusted FWHM initial guesses

        # Bounds: fwhm_gauss and fwhm_lorentz must be positive and within reasonable range
        # CORRECTED: xmax - x min changed to max_fwhm_guess (variable name)
        bounds = ([xmin, -np.inf, 0.01, 0.01, -np.inf, -np.inf, -np.inf],
                  [xmax, np.inf, max_fwhm_guess, max_fwhm_guess, np.inf, np.inf, np.inf])

        try:
            popt, pcov = optimize.curve_fit(lineshape, xs, ys, p0=p0, bounds=bounds, maxfev=50000)
            perr = np.sqrt(np.diag(pcov))

            nu_exp_val = popt[0]
            nu_exp_err = perr[0]
            peak_intensity_val = popt[1] # y0 is the amplitude from the fit
            peak_intensity_err = perr[1]
            fwhm_gauss_val = popt[2] 
            fwhm_gauss_err = perr[2]
            fwhm_lorentz_val = popt[3]
            fwhm_lorentz_err = perr[3]

            # --- Calculate Total Voigt FWHM and its error ---
            # F_V = 0.5346 * F_L + sqrt(0.2166 * F_L^2 + F_G^2)
            F_G = fwhm_gauss_val
            F_L = fwhm_lorentz_val

            # Prevent negative or zero under sqrt (shouldn't happen with positive F_L, F_G)
            sqrt_term = 0.2166 * F_L**2 + F_G**2
            if sqrt_term < 0: # Should not happen with positive F_L, F_G
                print(f"Warning: Negative term under sqrt for FM {fm_amp_mhz}. Setting Voigt FWHM to NaN.")
                fwhm_voigt_val = np.nan
                fwhm_voigt_err = np.nan
            else:
                fwhm_voigt_val = 0.5346 * F_L + np.sqrt(sqrt_term)

                # Error propagation for F_V using partial derivatives
                # Partial derivative d(FV)/d(FG)
                if np.isclose(np.sqrt(sqrt_term), 0): # Avoid division by zero if sqrt_term is exactly zero
                    dfv_dfg = 0 
                else:
                    dfv_dfg = F_G / np.sqrt(sqrt_term)

                # Partial derivative d(FV)/d(FL)
                if np.isclose(np.sqrt(sqrt_term), 0): # Avoid division by zero
                    dfv_dfl = 0.5346 # In this case, it's just the Lorentzian term
                else:
                    dfv_dfl = 0.5346 + (0.2166 * F_L) / np.sqrt(sqrt_term)

                fwhm_voigt_err = np.sqrt(
                    (dfv_dfl * fwhm_lorentz_err)**2 + 
                    (dfv_dfg * fwhm_gauss_err)**2
                )


            # --- DIAGNOSTIC PRINT ---
            print(f"\n--- Fit Results for FM Amplitude {fm_amp_mhz:.1f} MHz ---")
            print(f"  Peak Freq (x0): {nu_exp_val:.3f} +/- {nu_exp_err:.3f} MHz")
            print(f"  Peak Intensity (y0): {peak_intensity_val:.3f} +/- {peak_intensity_err:.3f}")
            print(f"  Gaussian FWHM: {fwhm_gauss_val:.3f} +/- {fwhm_gauss_err:.3f} MHz")
            print(f"  Lorentzian FWHM: {fwhm_lorentz_val:.3f} +/- {fwhm_lorentz_err:.3f} MHz")
            print(f"  Calculated Voigt FWHM: {fwhm_voigt_val:.3f} +/- {fwhm_voigt_err:.3f} MHz")
            # --- END DIAGNOSTIC PRINT ---

            results_fm_amplitude.append({
                'fm_mhz': fm_amp_mhz,
                'nu_exp': nu_exp_val,
                'nu_exp_error': nu_exp_err,
                'peak_intensity': peak_intensity_val,
                'peak_intensity_error': peak_intensity_err,
                'fwhm_gauss': fwhm_gauss_val, 
                'fwhm_gauss_error': fwhm_gauss_err,
                'fwhm_lorentz': fwhm_lorentz_val,
                'fwhm_lorentz_error': fwhm_lorentz_err,
                'fwhm_voigt': fwhm_voigt_val,
                'fwhm_voigt_error': fwhm_voigt_err
            })

            fit_xs = np.linspace(xs.min(), xs.max(), 1000)
            fit_ys = lineshape(fit_xs, *popt)

            # Plotting the data with the fit
            plot_data(xs, ys, 
                      title=f'{J_val},{K_val} - FM Amplitude: {fm_amp_mhz:.1f} MHz',
                      xlabel='Frequency (MHz)', 
                      ylabel='Intensity',
                      save_path=os.path.join(SECTION4_OUTPUT_FOLDER, f'fit_{J_val}_{K_val}_FM_Amp_{setting_str}.png'),
                      fit_xs=fit_xs, 
                      fit_ys=fit_ys, 
                      legend_labels=['Data', 'Fit'], 
                      show_legend=True)

            print(f"  Analysis and plot for FM Amplitude {fm_amp_mhz:.1f} MHz saved.")

        except Exception as e:
            print(f"  ERROR fitting and plotting for FM Amplitude {fm_amp_mhz:.1f} MHz: {e}. Skipping.")

    # --- Generate LaTeX Tables for Section 4 ---
    if results_fm_amplitude:
        with open(SECTION4_TABLES_FILE, 'w') as f:
            f.write("% --- LaTeX Code for Section 4 Tables (Requires \\usepackage{amsmath}) ---\n\n")

            # Table 3: Peak intensity and peak position at 4 different FM amplitudes.
            f.write("\\begin{table}[h]\n")
            f.write("    \\centering\n")
            f.write("    \\begin{tabular}{|c|c|c|}\n")
            f.write("        \\hline\n")
            f.write("        FM (MHz) & $\\nu$ (MHz) & Peak Intensity (Amplitude) \\\\\n")
            f.write("        \\hline\n")
            for res in results_fm_amplitude:
                f.write(f"        {res['fm_mhz']:.1f} & {res['nu_exp']:.3f}$\\pm${res['nu_exp_error']:.3f} & {res['peak_intensity']:.3f}$\\pm${res['peak_intensity_error']:.3f} \\\\\n")
            f.write("        \\hline\n")
            f.write("    \\end{tabular}\n")
            f.write("    \\caption{Peak intensity and peak position at 4 different FM amplitudes.}\n")
            f.write("    \\label{tab:fm_amplitude_peak_data}\n")
            f.write("\\end{table}\n\n")

            # Table 4: FWHM for different frequency modulations (now Voigt FWHM)
            f.write("\\begin{table}[h]\n")
            f.write("    \\centering\n")
            f.write("    \\begin{tabular}{|c|c|}\n") # Only one FWHM column
            f.write("        \\hline\n")
            f.write("        FM (MHz) & Voigt FWHM (MHz) \\\\\n") # Updated header for total Voigt FWHM
            f.write("        \\hline\n")
            for res in results_fm_amplitude:
                f.write(f"        {res['fm_mhz']:.1f} & {res['fwhm_voigt']:.3f}$\\pm${res['fwhm_voigt_error']:.3f} \\\\\n")
            f.write("        \\hline\n")
            f.write("    \\end{tabular}\n")
            f.write("    \\caption{Voigt FWHM for different frequency modulations.}\n") # Updated caption
            f.write("    \\label{tab:fm_amplitude_fwhm}\n")
            f.write("\\end{table}\n")
        
        print(f"\nLaTeX tables saved to '{SECTION4_TABLES_FILE}' successfully.")
    else:
        print("\nNo successful fits for Section 4, so no tables were generated.")

    print(f"Section 4 plots saved to: {SECTION4_OUTPUT_FOLDER}")


def run_section5_analysis():
    print("\n--- Section 5: Dependence of Line Shape and Intensity on Pressure (Transition 3,3) ---")

    # Reverted to analyzing only these three specific files as requested
    file_labels = ['30', '65', '95']
    
    J_val, K_val = 3, 3 # Specific transition for this section

    results_pressure_dependence = [] # To store data for tables and plotting FWHM vs Pressure

    for label in file_labels: # Iterate through the specified labels
        filename = f'{label}.csv'
        filepath = os.path.join(SECTION5_DATA_FOLDER, filename)
        
        xs, ys, pressure_mV_from_header = load_pressure_data_from_header(filepath)

        if xs is None or ys is None or np.isnan(pressure_mV_from_header):
            print(f"Skipping analysis for file {filename} due to data loading or pressure extraction error.")
            continue

        pressure_mbar = convert_voltage_to_pressure_mbar(pressure_mV_from_header)
        print(f"\nAnalyzing file {filename} (Pressure Sensor Reading: {pressure_mV_from_header:.0f}mV, Calculated Pressure: {pressure_mbar:.4f} mbar)")

        # --- Fitting the data ---
        xmin, xmax = xs.min(), xs.max()
        x0_guess_initial = (xmax + xmin) / 2
        ymin, ymax = ys.min(), ys.max()
        yptp = ymax - ymin

        xs_weight_factor = 4
        ys_weighted = (ys - ymin) * np.exp(- np.abs(np.abs(xs - x0_guess_initial) / (xmax - xmin)) * xs_weight_factor)
        x0_guess = xs[np.argmax(ys_weighted)]

        max_fwhm_guess = xmax - xmin 
        
        # Initial guesses for parameters: x0, y0, gauss, lorentz, a, b, c
        p0 = [x0_guess, yptp, 0.5, 0.5, 0, 0, 0] 

        bounds = ([xmin, -np.inf, 0.01, 0.01, -np.inf, -np.inf, -np.inf],
                  [xmax, np.inf, max_fwhm_guess, max_fwhm_guess, np.inf, np.inf, np.inf])

        try:
            popt, pcov = optimize.curve_fit(lineshape, xs, ys, p0=p0, bounds=bounds, maxfev=50000)
            perr = np.sqrt(np.diag(pcov))

            nu_exp_val = popt[0]
            nu_exp_err = perr[0]
            peak_intensity_val = popt[1]
            peak_intensity_err = perr[1]
            fwhm_gauss_val = popt[2] 
            fwhm_gauss_err = perr[2]
            fwhm_lorentz_val = popt[3]
            fwhm_lorentz_err = perr[3]

            # --- Calculate Total Voigt FWHM and its error ---
            F_G = fwhm_gauss_val
            F_L = fwhm_lorentz_val

            sqrt_term = 0.2166 * F_L**2 + F_G**2
            if sqrt_term < 0:
                print(f"Warning: Negative term under sqrt for Pressure {pressure_mbar:.4f} mbar. Setting FWHM to NaN.")
                fwhm_voigt_val = np.nan
                fwhm_voigt_err = np.nan
            else:
                fwhm_voigt_val = 0.5346 * F_L + np.sqrt(sqrt_term)

                if np.isclose(np.sqrt(sqrt_term), 0):
                    dfv_dfg = 0 
                    dfv_dfl = 0.5346 
                else:
                    dfv_dfg = F_G / np.sqrt(sqrt_term)
                    dfv_dfl = 0.5346 + (0.2166 * F_L) / np.sqrt(sqrt_term)

                fwhm_voigt_err = np.sqrt(
                    (dfv_dfl * fwhm_lorentz_err)**2 + 
                    (dfv_dfg * fwhm_gauss_err)**2
                )
            # --- End Voigt FWHM Calculation ---

            results_pressure_dependence.append({
                'filename': filename, 
                'pressure_mbar': pressure_mbar,
                'nu_exp': nu_exp_val,
                'nu_exp_error': nu_exp_err,
                'peak_intensity': peak_intensity_val,
                'peak_intensity_error': peak_intensity_err,
                'fwhm_voigt': fwhm_voigt_val, # Still storing as fwhm_voigt internally
                'fwhm_voigt_error': fwhm_voigt_err
            })

            fit_xs = np.linspace(xs.min(), xs.max(), 1000)
            fit_ys = lineshape(fit_xs, *popt)

            plot_data(xs, ys, 
                      title=f'{J_val},{K_val} - Pressure: {pressure_mbar:.4f} mbar',
                      xlabel='Frequency (MHz)', 
                      ylabel='Intensity',
                      save_path=os.path.join(SECTION5_OUTPUT_FOLDER, f'fit_{J_val}_{K_val}_Pressure_{label}.png'),
                      fit_xs=fit_xs, 
                      fit_ys=fit_ys, 
                      legend_labels=['Data', 'Fit'], 
                      show_legend=True)

            print(f"  Analysis and plot for Pressure {pressure_mbar:.4f} mbar saved.")

        except Exception as e:
            print(f"  ERROR fitting and plotting for Pressure from file {filename}: {e}. Skipping.")

    # --- Generate LaTeX Tables for Section 5 ---
    if results_pressure_dependence:
        # Sort results by pressure for clean tables and plots
        results_pressure_dependence.sort(key=lambda x: x['pressure_mbar'])

        with open(SECTION5_TABLES_FILE, 'w') as f:
            f.write("% --- LaTeX Code for Section 5 Tables (Requires \\usepackage{amsmath}) ---\n\n")

            # Table 5: Peak centres and amplitudes for various pressures
            f.write("\\begin{table}[h]\n")
            f.write("    \\centering\n")
            f.write("    \\begin{tabular}{|c|c|c|}\n")
            f.write("        \\hline\n")
            f.write("        Pressure (mbar) & $\\nu$ (MHz) & Peak Intensity (Amplitude) \\\\\n")
            f.write("        \\hline\n")
            for res in results_pressure_dependence:
                f.write(f"        {res['pressure_mbar']:.4f} & {res['nu_exp']:.3f}$\\pm${res['nu_exp_error']:.3f} & {res['peak_intensity']:.3f}$\\pm${res['peak_intensity_error']:.3f} \\\\\n")
            f.write("        \\hline\n")
            f.write("    \\end{tabular}\n")
            f.write("    \\caption{Peak centres and amplitudes for various pressures.}\n")
            f.write("    \\label{tab:pressure_peak_data}\n")
            f.write("\\end{table}\n\n")

            # Table 6: FWHM for different pressures - Changed "Voigt FWHM" to "FWHM"
            f.write("\\begin{table}[h]\n")
            f.write("    \\centering\n")
            f.write("    \\begin{tabular}{|c|c|}\n")
            f.write("        \\hline\n")
            f.write("        Pressure (mbar) & FWHM (MHz) \\\\\n") # Changed header here
            f.write("        \\hline\n")
            for res in results_pressure_dependence:
                f.write(f"        {res['pressure_mbar']:.4f} & {res['fwhm_voigt']:.3f}$\\pm${res['fwhm_voigt_error']:.3f} \\\\\n")
            f.write("        \\hline\n")
            f.write("    \\end{tabular}\n")
            f.write("    \\caption{FWHM for different pressures.}\n") # Changed caption here
            f.write("    \\label{tab:pressure_fwhm}\n")
            f.write("\\end{table}\n")
        
        print(f"\nSection 5 LaTeX tables saved to '{SECTION5_TABLES_FILE}' successfully.")

        # --- Plot FWHM vs. Pressure ---
        pressures = np.array([res['pressure_mbar'] for res in results_pressure_dependence])
        fwhms = np.array([res['fwhm_voigt'] for res in results_pressure_dependence])
        fwhm_errors = np.array([res['fwhm_voigt_error'] for res in results_pressure_dependence])

        # Attempt linear fit for FWHM vs Pressure
        if len(pressures) >= 2 and not np.any(np.isnan(fwhms)) and not np.all(fwhm_errors == 0): # Ensure at least 2 points for regression
            # Filter out NaNs or infs that could result from failed fits for plot/regression
            valid_indices = ~np.isnan(pressures) & ~np.isnan(fwhms) & ~np.isinf(fwhms) & ~np.isinf(pressures)
            pressures_valid = pressures[valid_indices]
            fwhms_valid = fwhms[valid_indices]
            fwhm_errors_valid = fwhm_errors[valid_indices]
            
            if len(pressures_valid) >= 2:
                # Weighted linear regression if errors are available
                if np.any(fwhm_errors_valid > 0): # Check if errors are non-zero for weighting
                    weights = 1 / fwhm_errors_valid**2
                    p_fit, cov_fit = np.polyfit(pressures_valid, fwhms_valid, 1, w=weights, cov=True)
                    slope, intercept = p_fit
                    slope_err, intercept_err = np.sqrt(np.diag(cov_fit))
                    
                else: # Unweighted if errors are all zero or not provided
                    slope, intercept, r_value, p_value, std_err = linregress(pressures_valid, fwhms_valid)
                    slope_err, intercept_err = std_err, float('nan') 
                
                fit_x = np.linspace(min(pressures_valid) * 0.9, max(pressures_valid) * 1.1, 100)
                fit_y = slope * fit_x + intercept

                plt.figure(figsize=(10, 6))
                # Changed from plt.errorbar to plt.plot to remove error bars, keeping circular markers
                plt.plot(pressures_valid, fwhms_valid, 'o', color='black', label='Experimental Data') 
                plt.plot(fit_x, fit_y, color='red', linestyle='--', label=f'Linear Fit: y = ({slope:.2f}$\\pm${slope_err:.2f})x + ({intercept:.3f}$\\pm${intercept_err:.3f})')
                
                plt.title(f'{J_val},{K_val} Transition: FWHM vs. Pressure', fontsize=24) # Changed title here
                plt.xlabel('Pressure (mbar)', fontsize=20)
                plt.ylabel('FWHM (MHz)', fontsize=20) # Changed ylabel here
                plt.legend(fontsize=16) 
                plt.tick_params(axis='both', which='major', labelsize=20)
                plt.grid(True, linestyle='--', alpha=0.6)
                plt.tight_layout()
                fwhm_plot_path = os.path.join(SECTION5_OUTPUT_FOLDER, f'FWHM_vs_Pressure_{J_val}_{K_val}.png')
                plt.savefig(fwhm_plot_path)
                plt.close()
                print(f"  FWHM vs. Pressure plot saved to: {fwhm_plot_path}")
                print(f"  Linear fit parameters: Slope = {slope:.2f} +/- {slope_err:.2f}, Intercept = {intercept:.3f} +/- {intercept_err:.3f}")
            else:
                print("  Not enough valid data points after filtering for plotting FWHM vs. Pressure or performing linear regression.")
        else:
            print("  Not enough valid data points to plot FWHM vs. Pressure or perform linear regression.")

    else:
        print("\nNo successful fits for Section 5, so no tables or plots were generated.")

    print(f"Section 5 plots and results saved to: {SECTION5_OUTPUT_FOLDER}")


# --- Main Program Entry Point ---
if __name__ == "__main__":
    print("--- Starting Full Ammonia Spectroscopy Analysis ---")
    run_section1_plots()
    run_section2_analysis()
    run_section3_analysis()
    run_section4_analysis()
    run_section5_analysis() # Call the new section
    print("\n--- Full Analysis Script Completed ---")
    print("\nREMINDER: Update the theoretical value for (7,6) in NU_THEORY_LOOKUP once you get the exact value from your tutor!")