//+------------------------------------------------------------------+
//|                                           EhlersTASCWrapper.mqh |
//|                        Copyright 2024, Your Name/Company        |
//+------------------------------------------------------------------+
//| Based on John F. Ehlers' "Measuring Cycle Periods" from         |
//| TASC March 2008 - Channelized Receiver Approach                 |
//| Compatible with JE_Adaptive_Norm_CG_DZ indicator                |
//+------------------------------------------------------------------+

#include <Arrays\ArrayDouble.mqh>
#include <Math\Stat\Math.mqh>

#ifndef EHLERS_TASC_WRAPPER_MQH
#define EHLERS_TASC_WRAPPER_MQH

// Define constants globally needed for array sizing and loop bounds
#define TASC_WRAPPER_NUM_PERIODS 51 // Size 0-50 for indexing
#define TASC_WRAPPER_MINN_CALC 8
#define TASC_WRAPPER_MAXN_CALC 50
#define TASC_WRAPPER_MINN_CG 10
#define TASC_WRAPPER_MEDIAN_LOOKBACK 10

class EhlersTASCWrapper
{
private:
   // --- Constants ---
   static const int NumPeriods;
   static const int MinN_Calc;
   static const int MaxN_Calc;
   static const int MinN_CG;
   static const int MedianLookback;
   static const double DefaultPeriod;

   // --- One-time Init Values ---
   double m_alpha1_hpf; // HPF alpha (uses fixed 40 period equivalent)
   double m_log10_val;

   // --- State Variables (History per Period N) ---
   double m_EhlersI_Prev1[TASC_WRAPPER_NUM_PERIODS];
   double m_EhlersI_Prev2[TASC_WRAPPER_NUM_PERIODS];
   double m_Q_Prev1[TASC_WRAPPER_NUM_PERIODS];
   double m_Q_Prev2[TASC_WRAPPER_NUM_PERIODS];
   double m_Real_Prev1[TASC_WRAPPER_NUM_PERIODS];
   double m_Real_Prev2[TASC_WRAPPER_NUM_PERIODS];
   double m_Imag_Prev1[TASC_WRAPPER_NUM_PERIODS];
   double m_Imag_Prev2[TASC_WRAPPER_NUM_PERIODS];

   // --- Internal Buffers (Self-Contained) ---
   double m_HP_Buffer[];
   double m_SmoothHP_Buffer[];
   double m_DC_Buffer[];
   
   // --- Buffer Management ---
   int m_buffer_size;
   bool m_buffers_initialized;

   // --- Initialization Flag ---
   bool m_initialized;

   // --- Private Median Helper ---
   double CalculateMedian(const double &arr[], int bar_index, int count);
   
   // --- Internal Buffer Management ---
   bool ResizeInternalBuffers(int new_size);

public:
   // --- Constructor ---
   EhlersTASCWrapper(void);
   // --- Destructor ---
   ~EhlersTASCWrapper(void);

   // --- Initialization ---
   bool Init(void);

   // --- Original Method (for JE_DC2 compatibility) ---
   double CalculateDominantCycle(const double &price[], // Price array (non-series)
                                const int bar,         // Current bar index 'i'
                                const int rates_total, // Total bars
                                // Indicator Buffers for HP, SmoothHP, DC history:
                                double &HP[],
                                double &SmoothHP[],
                                double &DC[]           // Also used for output DC[bar]
                                );

   // --- Simplified Method (for JE_Adaptive_Norm_CG_DZ compatibility) ---
   double CalculateDominantCycle(const double &price[], // Price array (non-series)
                                const int bar,         // Current bar index 'i'
                                const int rates_total  // Total bars
                                );

   // --- Reset Method ---
   void Reset();
};

// --- Define and Initialize static const members OUTSIDE the class ---
const int EhlersTASCWrapper::NumPeriods     = TASC_WRAPPER_NUM_PERIODS;
const int EhlersTASCWrapper::MinN_Calc      = TASC_WRAPPER_MINN_CALC;
const int EhlersTASCWrapper::MaxN_Calc      = TASC_WRAPPER_MAXN_CALC;
const int EhlersTASCWrapper::MinN_CG        = TASC_WRAPPER_MINN_CG;
const int EhlersTASCWrapper::MedianLookback = TASC_WRAPPER_MEDIAN_LOOKBACK;
const double EhlersTASCWrapper::DefaultPeriod = 20.0;

//+------------------------------------------------------------------+
//| Constructor                                                      |
//+------------------------------------------------------------------+
EhlersTASCWrapper::EhlersTASCWrapper(void) : 
   m_initialized(false), 
   m_alpha1_hpf(0.0), 
   m_log10_val(0.0),
   m_buffer_size(0),
   m_buffers_initialized(false)
{
   // Initialize state arrays
   ArrayInitialize(m_EhlersI_Prev1, 0.0);
   ArrayInitialize(m_EhlersI_Prev2, 0.0);
   ArrayInitialize(m_Q_Prev1, 0.0);
   ArrayInitialize(m_Q_Prev2, 0.0);
   ArrayInitialize(m_Real_Prev1, 0.0);
   ArrayInitialize(m_Real_Prev2, 0.0);
   ArrayInitialize(m_Imag_Prev1, 0.0);
   ArrayInitialize(m_Imag_Prev2, 0.0);
}

//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
EhlersTASCWrapper::~EhlersTASCWrapper(void)
{
   // ArrayFree will be called automatically for dynamic arrays
}

//+------------------------------------------------------------------+
//| Initialize the wrapper                                           |
//+------------------------------------------------------------------+
bool EhlersTASCWrapper::Init(void)
{
   if(m_initialized) return true; // Already initialized

   // Calculate HPF alpha (fixed 40 period / 9 degrees)
   double cos9 = MathCos(9.0 * M_PI / 180.0);
   if(MathAbs(cos9) > 1e-10) 
   {
      m_alpha1_hpf = (1.0 - MathSin(9.0 * M_PI / 180.0)) / cos9;
   } 
   else 
   {
      m_alpha1_hpf = 0.0;
      Print("Warning: EhlersTASCWrapper Cosine(9) is zero in initialization.");
   }
   
   // Calculate Log10
   m_log10_val = MathLog(10.0);
   if (MathAbs(m_log10_val) < 1e-10) 
   {
       Print("Error: EhlersTASCWrapper MathLog(10.0) is zero or too small!");
       m_log10_val = 2.302585; // Approx value if log fails
   }

   // Reset state arrays
   ArrayInitialize(m_EhlersI_Prev1, 0.0);
   ArrayInitialize(m_EhlersI_Prev2, 0.0);
   ArrayInitialize(m_Q_Prev1, 0.0);
   ArrayInitialize(m_Q_Prev2, 0.0);
   ArrayInitialize(m_Real_Prev1, 0.0);
   ArrayInitialize(m_Real_Prev2, 0.0);
   ArrayInitialize(m_Imag_Prev1, 0.0);
   ArrayInitialize(m_Imag_Prev2, 0.0);

   m_initialized = true;
   Print("EhlersTASCWrapper Initialized. HPF Alpha=", m_alpha1_hpf);
   return true;
}

//+------------------------------------------------------------------+
//| Internal Buffer Management                                       |
//+------------------------------------------------------------------+
bool EhlersTASCWrapper::ResizeInternalBuffers(int new_size)
{
   if(new_size <= 0) return false;
   
   if(!m_buffers_initialized || new_size != m_buffer_size)
   {
      if(ArrayResize(m_HP_Buffer, new_size) < 0) return false;
      if(ArrayResize(m_SmoothHP_Buffer, new_size) < 0) return false;
      if(ArrayResize(m_DC_Buffer, new_size) < 0) return false;
      
      // Initialize with safe values
      ArrayInitialize(m_HP_Buffer, 0.0);
      ArrayInitialize(m_SmoothHP_Buffer, 0.0);
      ArrayInitialize(m_DC_Buffer, DefaultPeriod);
      
      m_buffer_size = new_size;
      m_buffers_initialized = true;
   }
   
   return true;
}

//+------------------------------------------------------------------+
//| Private Median Helper Implementation                             |
//+------------------------------------------------------------------+
double EhlersTASCWrapper::CalculateMedian(const double &arr[], int bar_index, int count)
{
    int start_idx = bar_index - count + 1;
    // Robust bounds checking
    if(count <= 0 || start_idx < 0 || bar_index < 0 || start_idx >= ArraySize(arr) || 
       bar_index >= ArraySize(arr) || count > ArraySize(arr))
    {
        return (bar_index > 0 && bar_index < ArraySize(arr)) ? arr[bar_index] : DefaultPeriod;
    }

    double temp_arr[];
    if(ArrayResize(temp_arr, count) != count) 
        return (bar_index > 0 && bar_index < ArraySize(arr)) ? arr[bar_index] : DefaultPeriod;
    if(ArrayCopy(temp_arr, arr, 0, start_idx, count) != count) 
        return (bar_index > 0 && bar_index < ArraySize(arr)) ? arr[bar_index] : DefaultPeriod;

    ArraySort(temp_arr);

    if(count % 2 == 1) 
        return temp_arr[count / 2];
    else if(count >= 2) 
        return (temp_arr[count / 2 - 1] + temp_arr[count / 2]) / 2.0;
    else 
        return DefaultPeriod;
}

//+------------------------------------------------------------------+
//| Original Calculation Method (for JE_DC2 compatibility)          |
//+------------------------------------------------------------------+
double EhlersTASCWrapper::CalculateDominantCycle(
   const double &price[],
   const int bar,         // Current bar index 'i'
   const int rates_total, // Total bars
   // Indicator Buffers for HP, SmoothHP, DC history:
   double &HP[],
   double &SmoothHP[],
   double &DC[]           // Also used for output DC[bar]
   )
{
   // --- Check Initialization & History ---
   if(!m_initialized) 
   { 
      Print("EhlersTASCWrapper Error: Not initialized!"); 
      return 20.0; 
   }
   
   int min_hist_needed = MedianLookback - 1;
   min_hist_needed = MathMax(min_hist_needed, 5);

   if (bar < min_hist_needed) 
   {
       if(bar >= 0) 
       { 
           if(bar < ArraySize(HP)) HP[bar] = 0; 
           if(bar < ArraySize(SmoothHP)) SmoothHP[bar] = 0; 
           if(bar < ArraySize(DC)) DC[bar] = 20.0;
       } 
       return 20.0;
   }
   
   // Bounds check for required history access
   if (bar >= rates_total || bar-5 < 0 || bar >= ArraySize(HP) || bar-5 >= ArraySize(HP) || 
       bar >= ArraySize(SmoothHP) || bar-1 >= ArraySize(SmoothHP) || 
       bar >= ArraySize(DC) || bar-(MedianLookback-1) >= ArraySize(DC) || 
       bar >= ArraySize(price) || bar-1 >= ArraySize(price))
   {
        if(bar >= 0 && bar < ArraySize(DC)) 
            DC[bar] = (bar > 0 && bar-1 >= 0 && bar-1 < ArraySize(DC)) ? DC[bar-1] : 20.0;
        return (bar > 0 && bar < ArraySize(DC)) ? DC[bar] : 20.0;
   }

   // --- Calculate HP and SmoothHP for current bar 'bar' ---
   double hp_prev = HP[bar-1];
   HP[bar] = 0.5 * (m_alpha1_hpf + 1.0) * (price[bar] - price[bar-1]) + m_alpha1_hpf * hp_prev;
   SmoothHP[bar] = (HP[bar] + 2*HP[bar-1] + 3*HP[bar-2] + 3*HP[bar-3] + 2*HP[bar-4] + HP[bar-5]) / 12.0;

   // --- Delta (Bandwidth Factor) ---
   double EhlersDelta = -0.015 * (bar + 1) + 0.5; 
   EhlersDelta = MathMax(0.15, EhlersDelta);

   // --- Spectral Analysis Loop ---
   double MaxAmpl = 0.0; 
   double TempDB[TASC_WRAPPER_NUM_PERIODS];
   ArrayInitialize(TempDB, 60.0);

   for(int N = MinN_Calc; N <= MaxN_Calc; N++)
   {
      if (N < 0 || N >= NumPeriods) continue;

      // --- Coefficients ---
      double angle_rad_N = (N != 0) ? (2.0*M_PI / N) : 0.0; 
      double EhlersBeta = MathCos(angle_rad_N);
      double angle_rad_720_delta_N = (N != 0) ? (720.0*EhlersDelta*M_PI / (180.0*N)) : 0.0; 
      double Cos720Delta = MathCos(angle_rad_720_delta_N);
      double EhlersGamma = 1.0; 
      if(MathAbs(Cos720Delta) > 1e-10) 
         EhlersGamma = 1.0 / Cos720Delta;
      double gamma_sq = EhlersGamma*EhlersGamma; 
      double alpha_bp = 0.0; 
      if(gamma_sq >= 1.0) 
         alpha_bp = EhlersGamma - MathSqrt(gamma_sq - 1.0);
      alpha_bp = MathMax(0.0, MathMin(alpha_bp, 1.0)); 
      double OneMinusAlpha = 1.0-alpha_bp; 
      double OnePlusAlpha = 1.0+alpha_bp;

      // --- I / Q Calculation ---
      double smoothHP_curr = SmoothHP[bar]; 
      double smoothHP_prev = SmoothHP[bar-1];
      double q_val = (N != 0) ? (N/(2.0*M_PI))*(smoothHP_curr - smoothHP_prev) : 0.0; 
      double i_val = smoothHP_curr;

      // --- Real / Imag Calculation using internal state ---
      double i_val_bar2 = m_EhlersI_Prev2[N]; 
      double q_val_bar2 = m_Q_Prev2[N];
      double real_prev1 = m_Real_Prev1[N]; 
      double real_prev2 = m_Real_Prev2[N];
      double imag_prev1 = m_Imag_Prev1[N]; 
      double imag_prev2 = m_Imag_Prev2[N];
      double real_curr = 0.5*OneMinusAlpha*(i_val - i_val_bar2) + EhlersBeta*OnePlusAlpha*real_prev1 - alpha_bp*real_prev2;
      double imag_curr = 0.5*OneMinusAlpha*(q_val - q_val_bar2) + EhlersBeta*OnePlusAlpha*imag_prev1 - alpha_bp*imag_prev2;

      // --- Update State for next bar ---
      m_EhlersI_Prev2[N] = m_EhlersI_Prev1[N]; 
      m_EhlersI_Prev1[N] = i_val;
      m_Q_Prev2[N] = m_Q_Prev1[N]; 
      m_Q_Prev1[N] = q_val;
      m_Real_Prev2[N] = m_Real_Prev1[N]; 
      m_Real_Prev1[N] = real_curr;
      m_Imag_Prev2[N] = m_Imag_Prev1[N]; 
      m_Imag_Prev1[N] = imag_curr;

      // --- Amplitude & DB ---
      double ampl_sq = real_curr*real_curr + imag_curr*imag_curr; 
      if(ampl_sq > MaxAmpl) MaxAmpl = ampl_sq;
      double db_val = 60.0; 
      if(MaxAmpl > 1e-12) 
      { 
         double amplRatio = (MaxAmpl > 0) ? (ampl_sq / MaxAmpl) : 0.0; 
         if(amplRatio > 1e-10 && amplRatio < 1.0) 
         { 
            double logArg = 0.01/(1.0 - 0.99*amplRatio); 
            if(logArg > 1e-10 && MathAbs(m_log10_val) > 1e-10) 
            { 
               db_val = -10.0*MathLog(logArg)/m_log10_val; 
            } 
         } 
         else if (amplRatio >= 1.0) 
         { 
            db_val = 0.0; 
         } 
      }
      if(db_val > 20.0) db_val = 20.0; 
      if(db_val < 0.0) db_val = 0.0;
      
      // Bounds check for TempDB access
      if (N >= 0 && N < NumPeriods) 
      {
          TempDB[N] = db_val;
      }
   } // End N loop

   // --- Center of Gravity ---
   double Num_CG = 0.0; 
   double Denom_CG = 0.0;
   for(int N = MinN_CG; N <= MaxN_Calc; N++) 
   {
       if (N < 0 || N >= NumPeriods) continue;
       if(TempDB[N] <= 3.0) 
       { 
          Num_CG += N * (20.0 - TempDB[N]); 
          Denom_CG += (20.0 - TempDB[N]); 
       }
   }

   // --- Calculate DC[bar] ---
   double dc_prev = DC[bar-1];
   if (Denom_CG > 1e-10) 
   { 
      DC[bar] = Num_CG / Denom_CG; 
   }
   else 
   { 
      DC[bar] = dc_prev; 
   }

   // --- Apply median filter ---
   double finalDomCyc = CalculateMedian(DC, bar, MedianLookback);

   // --- Final Clamp ---
   if(finalDomCyc < 8.0) finalDomCyc = 8.0; 
   if(finalDomCyc > 50.0) finalDomCyc = 50.0;

   return finalDomCyc;
}

//+------------------------------------------------------------------+
//| Simplified Calculation Method (for JE_Adaptive_Norm_CG_DZ)      |
//+------------------------------------------------------------------+
double EhlersTASCWrapper::CalculateDominantCycle(
   const double &price[],   // Price array (non-series)
   const int bar,          // Current bar index 'i'
   const int rates_total   // Total bars
   )
{
   // --- Check Initialization & History ---
   if(!m_initialized) 
   { 
      Print("EhlersTASCWrapper Error: Not initialized!"); 
      return DefaultPeriod; 
   }
   
   // Ensure internal buffers are properly sized
   if(!ResizeInternalBuffers(rates_total))
   {
      Print("EhlersTASCWrapper Error: Failed to resize internal buffers!");
      return DefaultPeriod;
   }
   
   int min_hist_needed = MedianLookback - 1;
   min_hist_needed = MathMax(min_hist_needed, 6); // Need at least 6 bars for SmoothHP

   if (bar < min_hist_needed) 
   {
       if(bar >= 0 && bar < ArraySize(m_DC_Buffer)) 
       { 
           m_DC_Buffer[bar] = DefaultPeriod;
       } 
       return DefaultPeriod;
   }
   
   // Bounds check for required history access
   if (bar >= rates_total || bar-5 < 0 || bar >= ArraySize(price) || bar-1 >= ArraySize(price))
   {
        if(bar >= 0 && bar < ArraySize(m_DC_Buffer)) 
        {
            m_DC_Buffer[bar] = (bar > 0 && bar-1 >= 0 && bar-1 < ArraySize(m_DC_Buffer)) ? m_DC_Buffer[bar-1] : DefaultPeriod;
        }
        return (bar > 0 && bar < ArraySize(m_DC_Buffer)) ? m_DC_Buffer[bar] : DefaultPeriod;
   }

   // Call the original method using internal buffers
   return CalculateDominantCycle(price, bar, rates_total, m_HP_Buffer, m_SmoothHP_Buffer, m_DC_Buffer);
}

//+------------------------------------------------------------------+
//| Reset the wrapper state                                          |
//+------------------------------------------------------------------+
void EhlersTASCWrapper::Reset()
{
   ArrayInitialize(m_EhlersI_Prev1, 0.0);
   ArrayInitialize(m_EhlersI_Prev2, 0.0);
   ArrayInitialize(m_Q_Prev1, 0.0);
   ArrayInitialize(m_Q_Prev2, 0.0);
   ArrayInitialize(m_Real_Prev1, 0.0);
   ArrayInitialize(m_Real_Prev2, 0.0);
   ArrayInitialize(m_Imag_Prev1, 0.0);
   ArrayInitialize(m_Imag_Prev2, 0.0);
   
   // Reset internal buffers
   if(m_buffers_initialized)
   {
      ArrayInitialize(m_HP_Buffer, 0.0);
      ArrayInitialize(m_SmoothHP_Buffer, 0.0);
      ArrayInitialize(m_DC_Buffer, DefaultPeriod);
   }
}

#endif // EHLERS_TASC_WRAPPER_MQH