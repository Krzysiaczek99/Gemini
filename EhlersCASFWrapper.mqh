// EhlersCASFWrapper.mqh
// Wrapper for John Ehlers' CASF Dominant Cycle calculation (v1.1 - Fixes)

#ifndef EHLERS_CASF_WRAPPER_MQH
#define EHLERS_CASF_WRAPPER_MQH

#include <Math\Stat\Math.mqh>
#include <Arrays/Array.mqh> // For ArrayInitialize
//#include <Arrays/ArraySort.mqh> // For ArraySort in Median

// --- Constants ---
#define CASF_MEDIAN_LOOKBACK 5 // Lookback for MedianDelta calculation

//+------------------------------------------------------------------+
//| Median Helper Function (Standalone)                              |
//| Needs to be accessible by the wrapper                            |
//+------------------------------------------------------------------+
double CASF_MedianHelper(const double &arr[], int start_idx, int count)
{
   // ... (Median function remains the same as previous correct version) ...
   if(count <= 0 || start_idx < 0 || start_idx + count > ArraySize(arr)) { return 0.0; }
   double temp_arr[];
   if(ArrayResize(temp_arr, count) < 0) { Print("CASF_MedianHelper Error: Failed to resize temp array"); return 0.0; }
   bool copy_ok = true;
   for(int i=0; i<count; i++) {
       int source_idx = start_idx + i;
       if(source_idx >= ArraySize(arr)) { copy_ok=false; break; }
       double val = arr[source_idx];
       if(!MathIsValidNumber(val)) { PrintFormat("CASF_MedianHelper Warning: Invalid value arr[%d]=%.4g at source index %d. Returning 0.", i, val, source_idx); return 0.0; }
       temp_arr[i] = val;
   }
   if (!copy_ok) { Print("CASF_MedianHelper Error: Array index out of bounds during copy."); return 0.0; }
   ArraySort(temp_arr);
   if(count % 2 == 1) return temp_arr[count / 2];
   else if(count >= 2) return (temp_arr[count / 2 - 1] + temp_arr[count / 2]) / 2.0;
   else return 0.0;
}
//+------------------------------------------------------------------+
//| EhlersCASFWrapper Class                                          |
//+------------------------------------------------------------------+
class EhlersCASFWrapper
{
private:
   // --- Configuration ---
   double m_alpha;
   bool   m_initialized;

   // --- State Variables (History) ---
   double m_Smooth_prev1, m_Smooth_prev2;
   double m_Cycle_prev1, m_Cycle_prev2, m_Cycle_prev3, m_Cycle_prev4, m_Cycle_prev6;
   double m_Q1_prev1;
   double m_I1_prev1;
   double m_DeltaPhase_hist[CASF_MEDIAN_LOOKBACK];
   int    m_dp_hist_idx;
   double m_InstPeriod_prev1;
   double m_Period_prev1;

   // --- Constants (Declare static const members inside, define outside) ---
   static const int    MinHistory;
   static const double DefaultPeriod;
   static const double DefaultDeltaPhase;


public:
   // --- Constructor ---
                     EhlersCASFWrapper(void);
   // --- Destructor ---
                    ~EhlersCASFWrapper(void);

   // --- Initialization ---
   bool              Init(double alpha);

   // --- Calculation Method ---
   double            CalculateDominantCycle(const double &price[],
                                            const int bar,
                                            const int rates_total);
};

// --- Define and Initialize static const members OUTSIDE the class ---
const int    EhlersCASFWrapper::MinHistory       = 6;
const double EhlersCASFWrapper::DefaultPeriod    = 15.0;
const double EhlersCASFWrapper::DefaultDeltaPhase = 0.1;
// --- Constructor Implementation ---
EhlersCASFWrapper::EhlersCASFWrapper(void) :
   m_alpha(0.07), // Default alpha
   m_initialized(false),
   m_Smooth_prev1(0.0), m_Smooth_prev2(0.0),
   m_Cycle_prev1(0.0), m_Cycle_prev2(0.0), m_Cycle_prev3(0.0), m_Cycle_prev4(0.0), m_Cycle_prev6(0.0),
   m_Q1_prev1(0.0),
   m_I1_prev1(0.0),
   m_dp_hist_idx(0),
   m_InstPeriod_prev1(EhlersCASFWrapper::DefaultPeriod), // Use static const for init
   m_Period_prev1(EhlersCASFWrapper::DefaultPeriod)      // Use static const for init
{
   ArrayInitialize(m_DeltaPhase_hist, EhlersCASFWrapper::DefaultDeltaPhase); // Initialize history with default
}

// --- Destructor Implementation ---
EhlersCASFWrapper::~EhlersCASFWrapper(void) {}

// --- Initialization Implementation ---
bool EhlersCASFWrapper::Init(double alpha)
{
   if (alpha <= 0 || alpha >= 1.0) {
       Print("EhlersCASFWrapper Error: Alpha must be between 0 and 1 (exclusive). Received: ", alpha);
       m_initialized = false;
       return false;
   }
   m_alpha = alpha;

   // Reset state variables to defaults using static const
   m_Smooth_prev1 = 0.0; m_Smooth_prev2 = 0.0;
   m_Cycle_prev1 = 0.0; m_Cycle_prev2 = 0.0; m_Cycle_prev3 = 0.0; m_Cycle_prev4 = 0.0; m_Cycle_prev6 = 0.0;
   m_Q1_prev1 = 0.0;
   m_I1_prev1 = 0.0;
   ArrayInitialize(m_DeltaPhase_hist, EhlersCASFWrapper::DefaultDeltaPhase);
   m_dp_hist_idx = 0;
   m_InstPeriod_prev1 = EhlersCASFWrapper::DefaultPeriod;
   m_Period_prev1 = EhlersCASFWrapper::DefaultPeriod;

   m_initialized = true;
   Print("EhlersCASFWrapper Initialized. Alpha=", m_alpha);
   return true;
}
//+------------------------------------------------------------------+
//| CalculateDominantCycle Method                                    |
//+------------------------------------------------------------------+
double EhlersCASFWrapper::CalculateDominantCycle(
   const double &price[],
   const int bar,
   const int rates_total)
{
    // Use static const members defined outside
    bool print_casf_debug = (bar % 100 == 0 || bar >= rates_total - 5);

    // --- Check Initialization & History ---
    if (!m_initialized) { Print("EhlersCASFWrapper Error: Not initialized!"); return EhlersCASFWrapper::DefaultPeriod; }
    if (bar < EhlersCASFWrapper::MinHistory) { return EhlersCASFWrapper::DefaultPeriod; }
    if (bar >= rates_total || bar - EhlersCASFWrapper::MinHistory < 0) { Print("CASFWrapper Error [Bar ",bar,"]: Bar index or history index out of bounds."); return m_Period_prev1; }

    // --- Intermediate Variables ---
    double Smooth_curr = 0.0;
    double Cycle_curr = 0.0;
    double Q1_curr = 0.0;
    double I1_curr = 0.0;
    double DC_curr = EhlersCASFWrapper::DefaultPeriod;
    double DeltaPhase_curr = EhlersCASFWrapper::DefaultDeltaPhase;
    double MedianDelta_curr = EhlersCASFWrapper::DefaultDeltaPhase;
    double InstPeriod_curr = EhlersCASFWrapper::DefaultPeriod;
    double Period_curr = EhlersCASFWrapper::DefaultPeriod; // Initialize with default

    bool error_occurred = false; // Flag to skip subsequent steps if one fails
    double periodPrev = m_Period_prev1; // *** Capture previous period state *before* calculations ***
    double instPeriodPrev = m_InstPeriod_prev1; // Capture previous InstPeriod

    // 1. Smooth
    if (bar >= 3 && bar < rates_total && bar-1 >= 0 && bar-2 >= 0 && bar-3 >= 0) {
        Smooth_curr = (price[bar] + 2 * price[bar - 1] + 2 * price[bar - 2] + price[bar - 3]) / 6.0;
    } else { error_occurred = true; Smooth_curr = m_Smooth_prev1; }
    if (!MathIsValidNumber(Smooth_curr)) { Smooth_curr = m_Smooth_prev1; if(!MathIsValidNumber(Smooth_curr)) error_occurred=true; }

    // 2. Cycle
    if (!error_occurred && MathIsValidNumber(m_Smooth_prev1) && MathIsValidNumber(m_Smooth_prev2) && MathIsValidNumber(m_Cycle_prev1) && MathIsValidNumber(m_Cycle_prev2)) {
        double alpha_term = (1.0 - 0.5 * m_alpha); double alpha_sq = alpha_term * alpha_term; double alpha_m1 = (1.0 - m_alpha); double alpha_m1_sq = alpha_m1 * alpha_m1;
        Cycle_curr = alpha_sq * (Smooth_curr - 2 * m_Smooth_prev1 + m_Smooth_prev2) + 2 * alpha_m1 * m_Cycle_prev1 - alpha_m1_sq * m_Cycle_prev2;
    } else { error_occurred = true; }
    if (!MathIsValidNumber(Cycle_curr)) { Cycle_curr = m_Cycle_prev1; if(!MathIsValidNumber(Cycle_curr)) error_occurred=true; }

    // 3. Q1
    if (!error_occurred && MathIsValidNumber(m_Cycle_prev2) && MathIsValidNumber(m_Cycle_prev4) && MathIsValidNumber(m_Cycle_prev6) && MathIsValidNumber(instPeriodPrev)) {
        Q1_curr = (0.0962*Cycle_curr + 0.5769*m_Cycle_prev2 - 0.5769*m_Cycle_prev4 - 0.0962*m_Cycle_prev6) * (0.5 + 0.08*instPeriodPrev);
    } else { error_occurred = true; }
    if (!MathIsValidNumber(Q1_curr)) { Q1_curr = m_Q1_prev1; if(!MathIsValidNumber(Q1_curr)) error_occurred=true; }

    // 4. I1
    if (!error_occurred && MathIsValidNumber(m_Cycle_prev3)) {
        I1_curr = m_Cycle_prev3;
    } else { error_occurred = true; }
     if (!MathIsValidNumber(I1_curr)) { I1_curr = m_I1_prev1; if(!MathIsValidNumber(I1_curr)) error_occurred=true; }

    // 5. DeltaPhase
    DeltaPhase_curr = EhlersCASFWrapper::DefaultDeltaPhase; // Start with default
    if (!error_occurred && MathIsValidNumber(Q1_curr) && MathIsValidNumber(m_Q1_prev1) && MathIsValidNumber(I1_curr) && MathIsValidNumber(m_I1_prev1)) {
        if (MathAbs(Q1_curr) > 1e-10 && MathAbs(m_Q1_prev1) > 1e-10) {
            double term1 = I1_curr / Q1_curr; double term2 = m_I1_prev1 / m_Q1_prev1;
            if (MathIsValidNumber(term1) && MathIsValidNumber(term2)) {
                double denominator = (1.0 + term1 * term2);
                if (MathAbs(denominator) > 1e-10) {
                    double dp_calc = (term1 - term2) / denominator;
                    if (MathIsValidNumber(dp_calc)) DeltaPhase_curr = dp_calc; // Only assign if calculation valid
                }
            }
        }
    } // If inputs invalid, keeps default

    // 6. Clamp DeltaPhase
    DeltaPhase_curr = fmax(0.1, fmin(DeltaPhase_curr, 1.1));

    // Update DeltaPhase history
    m_DeltaPhase_hist[m_dp_hist_idx] = DeltaPhase_curr;
    m_dp_hist_idx = (m_dp_hist_idx + 1) % CASF_MEDIAN_LOOKBACK;

    // 7. MedianDelta
    MedianDelta_curr = CASF_MedianHelper(m_DeltaPhase_hist, 0, CASF_MEDIAN_LOOKBACK);
     if (!MathIsValidNumber(MedianDelta_curr)) MedianDelta_curr = EhlersCASFWrapper::DefaultDeltaPhase;

    // 8. DC Period
    DC_curr = EhlersCASFWrapper::DefaultPeriod; // Default
    if (MathAbs(MedianDelta_curr) > 1e-10) {
        DC_curr = (2.0 * M_PI) / MedianDelta_curr + 0.5;
        if (!MathIsValidNumber(DC_curr)) DC_curr = EhlersCASFWrapper::DefaultPeriod;
    }

    // 9. InstPeriod
    if (MathIsValidNumber(DC_curr) && MathIsValidNumber(instPeriodPrev)) {
        InstPeriod_curr = 0.33 * DC_curr + 0.67 * instPeriodPrev;
    } else { InstPeriod_curr = instPeriodPrev; }
    if (!MathIsValidNumber(InstPeriod_curr)) InstPeriod_curr = instPeriodPrev;

    // 10. Period
    if (MathIsValidNumber(InstPeriod_curr) && MathIsValidNumber(periodPrev)) {
        Period_curr = 0.15 * InstPeriod_curr + 0.85 * periodPrev;
    } else { Period_curr = periodPrev; }
    if (!MathIsValidNumber(Period_curr)) Period_curr = periodPrev;

    // --- Update State Variables for the *NEXT* calculation ---
    m_Smooth_prev2 = m_Smooth_prev1;
    m_Smooth_prev1 = Smooth_curr;
    m_Cycle_prev6 = m_Cycle_prev4;
    m_Cycle_prev4 = m_Cycle_prev3;
    m_Cycle_prev3 = m_Cycle_prev2;
    m_Cycle_prev2 = m_Cycle_prev1;
    m_Cycle_prev1 = Cycle_curr;
    m_Q1_prev1 = Q1_curr;
    m_I1_prev1 = I1_curr;
    m_InstPeriod_prev1 = InstPeriod_curr;
    m_Period_prev1 = Period_curr; // Update the final period state
    
    return Period_curr;
}
#endif // EHLERS_CASF_WRAPPER_MQH