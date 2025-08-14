// EhlersTASCWrapperSelfContained.mqh
// Self-contained TASC Dominant Cycle for EA use

#ifndef EHLERS_TASC_WRAPPER_SELFCONTAINED_MQH
#define EHLERS_TASC_WRAPPER_SELFCONTAINED_MQH

#include <Math\Stat\Math.mqh>
#include <Arrays\ArrayDouble.mqh>

// Define constants (can be adjusted or made inputs to Init if needed)
#define TASC_SC_NUM_PERIODS 51      // Max period for internal arrays
#define TASC_SC_MINN_CALC 8
#define TASC_SC_MAXN_CALC 50
#define TASC_SC_MINN_CG 10
#define TASC_SC_MEDIAN_LOOKBACK 10
#define TASC_SC_HPF_PERIOD_FOR_ALPHA 40 // Corresponds to 9 degrees for HPF alpha

class EhlersTASCWrapperSelfContained
{
private:
   // --- Configuration ---
   bool   m_initialized;
   double m_alpha1_hpf; // HPF alpha
   double m_log10_val;

   // --- State Variables (Internal History) ---
   // For HP and SmoothHP calculations
   double m_price_prev1;
   double m_hp_prev1, m_hp_prev2, m_hp_prev3, m_hp_prev4, m_hp_prev5;
   double m_smoothHP_prev1;

   // For Spectral Analysis (per period N)
   double m_EhlersI_Prev1[TASC_SC_NUM_PERIODS];
   double m_EhlersI_Prev2[TASC_SC_NUM_PERIODS];
   double m_Q_Prev1[TASC_SC_NUM_PERIODS];
   double m_Q_Prev2[TASC_SC_NUM_PERIODS];
   double m_Real_Prev1[TASC_SC_NUM_PERIODS];
   double m_Real_Prev2[TASC_SC_NUM_PERIODS];
   double m_Imag_Prev1[TASC_SC_NUM_PERIODS];
   double m_Imag_Prev2[TASC_SC_NUM_PERIODS];

   // For DC and Median DC
   double m_DC_hist[TASC_SC_MEDIAN_LOOKBACK]; // History of calculated raw DCs
   int    m_dc_hist_idx;                     // Circular buffer index for DC_hist
   double m_medianDC_prev1;                  // Previous bar's final median filtered DC

   // --- Private Median Helper (using standard indexing) ---
   double CalculateMedian(const double &arr[], int count); // Takes full array, uses first 'count' elements

public:
   // --- Constructor ---
                     EhlersTASCWrapperSelfContained(void);
   // --- Destructor ---
                    ~EhlersTASCWrapperSelfContained(void);

   // --- Initialization ---
   bool              Init(void); // Could add HPF period as param if desired

   // --- Calculation Method ---
   double            CalculateDominantCycle(const double &price_current, // ONLY current bar's price
                                            const int bar_number       // Current absolute bar number (for Delta calc)
                                           );
};

// --- Constructor Implementation ---
EhlersTASCWrapperSelfContained::EhlersTASCWrapperSelfContained(void) :
   m_initialized(false), m_alpha1_hpf(0.0), m_log10_val(0.0),
   m_price_prev1(0.0),
   m_hp_prev1(0.0), m_hp_prev2(0.0), m_hp_prev3(0.0), m_hp_prev4(0.0), m_hp_prev5(0.0),
   m_smoothHP_prev1(0.0),
   m_dc_hist_idx(0),
   m_medianDC_prev1(20.0) // Default starting period
{
   ArrayInitialize(m_EhlersI_Prev1, 0.0); ArrayInitialize(m_EhlersI_Prev2, 0.0);
   ArrayInitialize(m_Q_Prev1, 0.0); ArrayInitialize(m_Q_Prev2, 0.0);
   ArrayInitialize(m_Real_Prev1, 0.0); ArrayInitialize(m_Real_Prev2, 0.0);
   ArrayInitialize(m_Imag_Prev1, 0.0); ArrayInitialize(m_Imag_Prev2, 0.0);
   ArrayInitialize(m_DC_hist, 20.0); // Initialize DC history with default
}

// --- Destructor Implementation ---
EhlersTASCWrapperSelfContained::~EhlersTASCWrapperSelfContained(void) {}

// --- Initialization Implementation ---
bool EhlersTASCWrapperSelfContained::Init(void)
{
   if(m_initialized) return true;

   // Calculate HPF alpha (fixed 40 period equivalent / 9 degrees)
   double angle_rad_hpf = (2.0 * M_PI / TASC_SC_HPF_PERIOD_FOR_ALPHA); // For a 40-period HPF
   double cos_val = MathCos(angle_rad_hpf);
   if(MathAbs(cos_val) > 1e-10) {
      m_alpha1_hpf = (1.0 - MathSin(angle_rad_hpf)) / cos_val;
   } else {
      m_alpha1_hpf = 0.0; // Default if cosine is zero
      Print("Warning: EhlersTASCWrapperSelfContained Cosine for HPF Alpha is zero in Init.");
   }
   // Ensure alpha is reasonable (e.g., for 40 period, it's around 0.7 to 0.8)
   // If using Ehlers' (1-sin(deg))/(cos(deg)) for 9 degrees (equiv to ~40 period HPF):
   double deg9_rad = 9.0 * M_PI / 180.0;
   double cos9 = MathCos(deg9_rad);
   if(MathAbs(cos9) > 1e-10) m_alpha1_hpf = (1.0 - MathSin(deg9_rad)) / cos9;
   else m_alpha1_hpf = 0.7; // Fallback


   m_log10_val = MathLog(10.0);
   if (MathAbs(m_log10_val) < 1e-10) {
       Print("Error: EhlersTASCWrapperSelfContained MathLog(10.0) is zero!");
       m_log10_val = 2.302585092994046; // Fallback
   }

   // Reset state arrays
   m_price_prev1 = 0.0;
   m_hp_prev1 = 0.0; m_hp_prev2 = 0.0; m_hp_prev3 = 0.0; m_hp_prev4 = 0.0; m_hp_prev5 = 0.0;
   m_smoothHP_prev1 = 0.0;
   ArrayInitialize(m_EhlersI_Prev1, 0.0); ArrayInitialize(m_EhlersI_Prev2, 0.0);
   ArrayInitialize(m_Q_Prev1, 0.0); ArrayInitialize(m_Q_Prev2, 0.0);
   ArrayInitialize(m_Real_Prev1, 0.0); ArrayInitialize(m_Real_Prev2, 0.0);
   ArrayInitialize(m_Imag_Prev1, 0.0); ArrayInitialize(m_Imag_Prev2, 0.0);
   ArrayInitialize(m_DC_hist, 20.0);
   m_dc_hist_idx = 0;
   m_medianDC_prev1 = 20.0;


   m_initialized = true;
   Print("EhlersTASCWrapperSelfContained Initialized. HPF Alpha=", m_alpha1_hpf);
   return true;
}

// --- Private Median Helper Implementation ---
double EhlersTASCWrapperSelfContained::CalculateMedian(const double &arr[], int count)
{
    if(count <= 0 || count > ArraySize(arr)) { return 0.0; } // Basic check
    double temp_arr[];
    if(ArrayResize(temp_arr, count) < 0) return 0.0;
    if(ArrayCopy(temp_arr, arr, 0, 0, count) != count) return 0.0;
    ArraySort(temp_arr);
    if(count % 2 == 1) return temp_arr[count / 2];
    else if(count >= 2) return (temp_arr[count / 2 - 1] + temp_arr[count / 2]) / 2.0;
    return 0.0;
}

// --- Calculation Method Implementation ---
double EhlersTASCWrapperSelfContained::CalculateDominantCycle(
   const double &price_current, // Only current bar's price (e.g., priceFinal)
   const int bar_number        // Absolute bar number (e.g., Bars(_Symbol,_Period) - 1 for current)
   )
{
   if(!m_initialized) { Print("EhlersTASCWrapperSelfContained Error: Not initialized!"); return m_medianDC_prev1; }
   // Min history is implicitly handled by how many previous values are stored (m_hp_prev5 etc.)

   // --- Calculate HP for current bar ---
   double hp_curr = 0.0;
   if (bar_number > 0) { // Need at least one previous price
       hp_curr = 0.5 * (m_alpha1_hpf + 1.0) * (price_current - m_price_prev1) + m_alpha1_hpf * m_hp_prev1;
   }
   if (!MathIsValidNumber(hp_curr)) hp_curr = m_hp_prev1;


   // --- Calculate SmoothHP for current bar ---
   double smoothHP_curr = 0.0;
   if (bar_number >= 5) { // Need 5 previous HP values
        smoothHP_curr = (hp_curr + 2*m_hp_prev1 + 3*m_hp_prev2 + 3*m_hp_prev3 + 2*m_hp_prev4 + m_hp_prev5) / 12.0;
   } else { smoothHP_curr = hp_curr; } // Fallback for early bars
   if (!MathIsValidNumber(smoothHP_curr)) smoothHP_curr = m_smoothHP_prev1;


   // --- Delta (Bandwidth Factor) ---
   // bar_number should be the 0-based index of the current bar being calculated
   double EhlersDelta = -0.015 * (bar_number + 1.0) + 0.5; // +1 because EL CurrentBar is 1-based
   EhlersDelta = MathMax(0.15, EhlersDelta);

   // --- Spectral Analysis Loop ---
   double MaxAmpl = 0.0;
   double TempDB[TASC_SC_NUM_PERIODS]; ArrayInitialize(TempDB, 60.0);

   for(int N = TASC_SC_MINN_CALC; N <= TASC_SC_MAXN_CALC; N++)
   {
      if (N < 0 || N >= TASC_SC_NUM_PERIODS) continue;

      double angle_rad_N = (N != 0) ? (2.0*M_PI / N) : 0.0; double EhlersBeta = MathCos(angle_rad_N);
      double angle_rad_720_delta_N = (N != 0) ? (720.0*EhlersDelta*M_PI / (180.0*N)) : 0.0; double Cos720Delta = MathCos(angle_rad_720_delta_N);
      double EhlersGamma = 1.0; if(MathAbs(Cos720Delta) > 1e-10) EhlersGamma = 1.0 / Cos720Delta; else EhlersGamma = 1e9; // Avoid div by zero
      double gamma_sq = EhlersGamma*EhlersGamma; double alpha_bp = 0.0;
      if(gamma_sq >= 1.0) alpha_bp = EhlersGamma - MathSqrt(gamma_sq - 1.0);
      else alpha_bp = EhlersGamma - MathSqrt(1.0 - gamma_sq); // Should not happen if Cos720Delta is not > 1
      alpha_bp = MathMax(0.0, MathMin(alpha_bp, 1.0)); double OneMinusAlpha = 1.0-alpha_bp; double OnePlusAlpha = 1.0+alpha_bp;

      double q_val = (N != 0) ? (N/(2.0*M_PI))*(smoothHP_curr - m_smoothHP_prev1) : 0.0; double i_val = smoothHP_curr;

      double i_val_bar2 = m_EhlersI_Prev2[N]; double q_val_bar2 = m_Q_Prev2[N];
      double real_prev1 = m_Real_Prev1[N]; double real_prev2 = m_Real_Prev2[N];
      double imag_prev1 = m_Imag_Prev1[N]; double imag_prev2 = m_Imag_Prev2[N];
      double real_curr = 0.5*OneMinusAlpha*(i_val - i_val_bar2) + EhlersBeta*OnePlusAlpha*real_prev1 - alpha_bp*real_prev2;
      double imag_curr = 0.5*OneMinusAlpha*(q_val - q_val_bar2) + EhlersBeta*OnePlusAlpha*imag_prev1 - alpha_bp*imag_prev2;

      m_EhlersI_Prev2[N] = m_EhlersI_Prev1[N]; m_EhlersI_Prev1[N] = i_val;
      m_Q_Prev2[N] = m_Q_Prev1[N]; m_Q_Prev1[N] = q_val;
      m_Real_Prev2[N] = m_Real_Prev1[N]; m_Real_Prev1[N] = real_curr;
      m_Imag_Prev2[N] = m_Imag_Prev1[N]; m_Imag_Prev1[N] = imag_curr;

      double ampl_sq = real_curr*real_curr + imag_curr*imag_curr; if(ampl_sq > MaxAmpl) MaxAmpl = ampl_sq;
      double db_val = 60.0;
      if(MaxAmpl > 1e-12) {
          double amplRatio = (MaxAmpl > 0) ? (ampl_sq / MaxAmpl) : 0.0;
          if(amplRatio > 1e-10 && amplRatio < 1.0) {
              double logArg = 0.01/(1.0 - 0.99*amplRatio);
              if(logArg > 1e-10 && MathAbs(m_log10_val) > 1e-10) { db_val = -10.0*MathLog(logArg)/m_log10_val; }
          } else if (amplRatio >= 1.0) { db_val = 0.0; }
      }
      db_val = fmax(0.0, fmin(db_val, 20.0)); // Clamp dB
      TempDB[N] = db_val;
   }

   // --- Center of Gravity ---
   double Num_CG = 0.0; double Denom_CG = 0.0;
   for(int N = TASC_SC_MINN_CG; N <= TASC_SC_MAXN_CALC; N++) {
       if (N < 0 || N >= TASC_SC_NUM_PERIODS) continue;
       if(TempDB[N] <= 3.0) { Num_CG += N * (20.0 - TempDB[N]); Denom_CG += (20.0 - TempDB[N]); }
   }

   // --- Calculate Raw DC for current bar ---
   double dc_raw_curr = m_medianDC_prev1; // Default to previous final DC
   if (Denom_CG > 1e-10) { dc_raw_curr = Num_CG / Denom_CG; }
   if (!MathIsValidNumber(dc_raw_curr)) dc_raw_curr = m_medianDC_prev1;

   // Store raw DC in history for median calculation
   m_DC_hist[m_dc_hist_idx] = dc_raw_curr;
   m_dc_hist_idx = (m_dc_hist_idx + 1) % TASC_SC_MEDIAN_LOOKBACK;

   // --- Apply median filter ---
   double finalDomCyc = CalculateMedian(m_DC_hist, TASC_SC_MEDIAN_LOOKBACK);

   // --- Final Clamp ---
   finalDomCyc = fmax(8.0, fmin(finalDomCyc, 50.0));
   if (!MathIsValidNumber(finalDomCyc)) finalDomCyc = m_medianDC_prev1; // Ultimate fallback

   // --- Update State Variables for *NEXT* bar ---
   m_price_prev1 = price_current;
   m_hp_prev5 = m_hp_prev4; m_hp_prev4 = m_hp_prev3; m_hp_prev3 = m_hp_prev2;
   m_hp_prev2 = m_hp_prev1; m_hp_prev1 = hp_curr;
   m_smoothHP_prev1 = smoothHP_curr;
   m_medianDC_prev1 = finalDomCyc; // Store the final calculated cycle

   return finalDomCyc;
}

#endif // EHLERS_TASC_WRAPPER_SELFCONTAINED_MQH