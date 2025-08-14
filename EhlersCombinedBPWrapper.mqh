// EhlersCombinedBPWrapper.mqh (Correcting Power Calc)
#include <Arrays\ArrayDouble.mqh>
#include <Math\Stat\Math.mqh>

#ifndef EHLERS_COMBINEDBP_WRAPPER_MQH
#define EHLERS_COMBINEDBP_WRAPPER_MQH

#define COMBINEDBP_MAX_PERIOD_CLASS 48
#define COMBINEDBP_MIN_PERIOD_CLASS 10
#define COMBINEDBP_POWER_SMOOTH_PERIOD 5 // EMA Period for Power Smoothing (adjust as needed)

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class EhlersCombinedBPWrapper
  {
private:
   // Filter State
   double m_BP_Prev1[COMBINEDBP_MAX_PERIOD_CLASS + 1];
   double m_BP_Prev2[COMBINEDBP_MAX_PERIOD_CLASS + 1];
   // Power State (Smoothed Power)
   double m_Pwr[COMBINEDBP_MAX_PERIOD_CLASS + 1]; // Now stores SMOOTHED power
   double            m_MaxPwr;
   bool              m_initialized;
   // Coefficients
   double m_g_alpha[COMBINEDBP_MAX_PERIOD_CLASS + 1];
   double m_g_beta[COMBINEDBP_MAX_PERIOD_CLASS + 1];
   double            m_delta;
   // Power Smoothing Alpha
   double            m_power_alpha;

public:
                     EhlersCombinedBPWrapper(void);
                    ~EhlersCombinedBPWrapper(void);
   bool              Init(const double delta);
   double            CalculateDominantCycle(const double &price[], const int bar, const int rates_total, const int spec_dil_comp_param);
private:
   void              InitializeArraysAndCoeffs(void);
  };

// --- Constructor ---
EhlersCombinedBPWrapper::EhlersCombinedBPWrapper(void) : m_initialized(false), m_MaxPwr(0.0), m_delta(0.0)
  {
   ArrayInitialize(m_BP_Prev1, 0.0);
   ArrayInitialize(m_BP_Prev2, 0.0);
   ArrayInitialize(m_Pwr, 0.0); // Initialize smoothed power state
   ArrayInitialize(m_g_alpha, 0.0);
   ArrayInitialize(m_g_beta, 0.0);
   m_power_alpha = 2.0 / (COMBINEDBP_POWER_SMOOTH_PERIOD + 1.0); // EMA alpha for power
  }
// --- Destructor ---
EhlersCombinedBPWrapper::~EhlersCombinedBPWrapper(void) { }
// --- Init ---
bool EhlersCombinedBPWrapper::Init(const double delta) { if(delta <= 1e-9 || delta >= 0.5) { return false; } if(m_initialized) Print("Warning: BP Wrapper Re-initializing"); m_delta = delta; InitializeArraysAndCoeffs(); m_initialized = true; Print("EhlersCombinedBPWrapper Initialized with Delta: ", m_delta); return true; }
// --- InitializeArraysAndCoeffs ---
void EhlersCombinedBPWrapper::InitializeArraysAndCoeffs(void) { ArrayInitialize(m_Pwr, 0.0); m_MaxPwr = 0.0; ArrayInitialize(m_BP_Prev1, 0.0); ArrayInitialize(m_BP_Prev2, 0.0); ArrayInitialize(m_g_alpha, 0.0); ArrayInitialize(m_g_beta, 0.0); Print("CombinedBP Wrapper: Pre-calculating coefficients..."); for(int N = COMBINEDBP_MIN_PERIOD_CLASS; N <= COMBINEDBP_MAX_PERIOD_CLASS; N++) { if(N <= 0) continue; double gamma = 1.0; double cosArg1 = 4.0*M_PI*m_delta/N; double cosVal1 = MathCos(cosArg1); if(MathAbs(cosVal1) > 1.0) { gamma = 1.0; } else if(MathAbs(cosVal1) < 1e-9) { gamma = 1e9; } else { gamma = 1.0 / cosVal1; } double beta_val = MathCos(2.0*M_PI/N); double alpha_val; double gamma_sq = gamma*gamma; if(gamma_sq - 1.0 <= 1e-9 || gamma <=1.0) { alpha_val = 1.0; } else { alpha_val = gamma - MathSqrt(gamma_sq - 1.0); } alpha_val = MathMax(1e-9, MathMin(alpha_val, 0.9999)); m_g_beta[N] = beta_val; m_g_alpha[N] = alpha_val; } Print("CombinedBP Wrapper: Coefficients pre-calculated."); }

// --- Calculation Method Implementation with Power Smoothing ---
double EhlersCombinedBPWrapper::CalculateDominantCycle(
   const double &price[],
   const int bar,
   const int rates_total,
   const int spec_dil_comp_param)
  {
   bool print_debug = (bar >= rates_total - 5 || bar % 100 == 0);

// ... (Initial checks remain the same) ...
   if(!m_initialized)
      return 0.0;
   if(bar < 2)
      return (COMBINEDBP_MIN_PERIOD_CLASS+COMBINEDBP_MAX_PERIOD_CLASS)/2.0;
   if(bar>=rates_total || bar-2<0)
      return (COMBINEDBP_MIN_PERIOD_CLASS+COMBINEDBP_MAX_PERIOD_CLASS)/2.0;

   double Filt = price[bar];
   double Filt_prev2 = price[bar-2];
   double Comp;

// Calculate new values for each period N
   for(int N = COMBINEDBP_MIN_PERIOD_CLASS; N <= COMBINEDBP_MAX_PERIOD_CLASS; N++)
     {
      if(N <= 0 || N > COMBINEDBP_MAX_PERIOD_CLASS)
         continue;

      // Get coefficients
      double alpha = m_g_alpha[N];
      double beta  = m_g_beta[N];

      // Calculate BP[N] using Mladen's filter
      double bp = 0.0;
      double bp_prev1 = m_BP_Prev1[N];
      double bp_prev2 = m_BP_Prev2[N];
      if(!MathIsValidNumber(bp_prev1))
         bp_prev1 = 0.0;
      if(!MathIsValidNumber(bp_prev2))
         bp_prev2 = 0.0;
      bp = 0.5*(1.0-alpha)*(Filt - Filt_prev2) + beta*(1.0+alpha)*bp_prev1 - alpha*bp_prev2;
      if(!MathIsValidNumber(bp))
         bp = 0.0;

      // Update State
      m_BP_Prev2[N] = bp_prev1;
      m_BP_Prev1[N] = bp;

      // --- Calculate Instantaneous Power & Smooth it ---
      Comp = (spec_dil_comp_param == 1) ? 1.0 : (double)N;
      if(Comp == 0)
         Comp = 1.0;
      double instant_pwr = (bp / Comp) * (bp / Comp);
      if(!MathIsValidNumber(instant_pwr) || instant_pwr < 0)
         instant_pwr = 0.0;

      // EMA Smoothing of Power
      double prev_smooth_pwr = m_Pwr[N]; // Get previous smoothed power from state
      if(!MathIsValidNumber(prev_smooth_pwr))
         prev_smooth_pwr = 0.0;
      double current_smooth_pwr = m_power_alpha * instant_pwr + (1.0 - m_power_alpha) * prev_smooth_pwr;

      // Store smoothed power back into state
      m_Pwr[N] = current_smooth_pwr;
      // --- End Power Calculation ---


      // --- DEBUG PRINT ---
      if(print_debug && (N == 10 || N == 20 || N == 40 || N == 48))
        {
         string bp_str = MathIsValidNumber(bp) ? StringFormat("%.4f", bp) : "INVALID";
         string pwr_str = (N>=0 && N<=COMBINEDBP_MAX_PERIOD_CLASS && MathIsValidNumber(m_Pwr[N])) ? StringFormat("%.6f", m_Pwr[N]) : "INVALID";
         PrintFormat("BPWrapper Calc (Mladen Filt+PwrEMA) [Bar %d, N=%d]: alpha=%.4f, beta=%.3f, BP= %s, SmPwr= %s", bar, N, alpha, beta, bp_str, pwr_str);
        }
     } // End N loop

// --- Post-processing (Uses the SMOOTHED m_Pwr now) ---
// ... Max Power calc ...
   double prevMaxPwr = m_MaxPwr;
   if(!MathIsValidNumber(m_MaxPwr) || m_MaxPwr < 0)
      m_MaxPwr = 0.0;
   m_MaxPwr *= 0.995;
   for(int Period = COMBINEDBP_MIN_PERIOD_CLASS; Period <= COMBINEDBP_MAX_PERIOD_CLASS; Period++)
     {
      if(Period >= 0 && Period <= COMBINEDBP_MAX_PERIOD_CLASS)
        {
         if(MathIsValidNumber(m_Pwr[Period]) && m_Pwr[Period] > m_MaxPwr)
           {
            m_MaxPwr = m_Pwr[Period];
           }
        }
     }
   if(!MathIsValidNumber(m_MaxPwr) || m_MaxPwr < 0)
      m_MaxPwr = 0.0;
//if (print_debug) { PrintFormat("BPWrapper Calc [Bar %d]: PrevMaxPwr=%.6g, DecayedMaxPwr=%.6g, NewMaxPwr=%.6g", bar, prevMaxPwr, prevMaxPwr * 0.995, m_MaxPwr); }
// ... Normalize Power ...
   double PwrNorm[COMBINEDBP_MAX_PERIOD_CLASS + 1];
   ArrayInitialize(PwrNorm, 0.0);
   if(m_MaxPwr > 1e-12)
     {
      for(int Period = COMBINEDBP_MIN_PERIOD_CLASS; Period <= COMBINEDBP_MAX_PERIOD_CLASS; Period++)
        {
         if(Period >= 0 && Period <= COMBINEDBP_MAX_PERIOD_CLASS)
           {
            if(MathIsValidNumber(m_Pwr[Period]))
              {
               PwrNorm[Period] = m_Pwr[Period] / m_MaxPwr;
              }
            else
              {
               PwrNorm[Period] = 0.0;
              }
           }
        }
     }
//if (print_debug) { /* PwrNorm printing */ }
// --- Compute Center of Gravity ---
double Spx = 0.0, Sp = 0.0;
double cg_threshold = 0.1; // <<< LOWER THRESHOLD (e.g., 10%)
for(int Period = COMBINEDBP_MIN_PERIOD_CLASS; Period <= COMBINEDBP_MAX_PERIOD_CLASS; Period++) {
    if (Period >= 0 && Period <= COMBINEDBP_MAX_PERIOD_CLASS) {
        if(MathIsValidNumber(PwrNorm[Period]) && PwrNorm[Period] >= cg_threshold) { // <<< USE LOWER THRESHOLD
            Spx += (double)Period * PwrNorm[Period];
            Sp += PwrNorm[Period];
        }
    }
}
//if (print_debug) { PrintFormat("BPWrapper Calc [Bar %d]: Spx=%.4f, Sp=%.4f", bar, Spx, Sp); }
// ... Final DC ...
   double dc = (Sp > 1e-12) ? Spx / Sp : (COMBINEDBP_MIN_PERIOD_CLASS + COMBINEDBP_MAX_PERIOD_CLASS) / 2.0;
   double dc_unclamped = dc;
   if(dc < COMBINEDBP_MIN_PERIOD_CLASS)
      dc = COMBINEDBP_MIN_PERIOD_CLASS;
   if(dc > COMBINEDBP_MAX_PERIOD_CLASS)
      dc = COMBINEDBP_MAX_PERIOD_CLASS;
   if(print_debug)
     {
      PrintFormat("BPWrapper Calc [Bar %d] End: UnclampedDC=%.2f, FinalDC=%.2f", bar, dc_unclamped, dc);
     }

   return dc;
  }

#endif // EHLERS_COMBINEDBP_WRAPPER_MQH
//+------------------------------------------------------------------+
