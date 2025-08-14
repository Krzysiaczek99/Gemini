// EhlersMESAWrapper.mqh (Using Standard Arrays)
#include <Arrays\ArrayDouble.mqh> // Still needed for ArraySort maybe, but not core arrays
#include <Math\Stat\Math.mqh>

#ifndef EHLERS_MESA_WRAPPER_MQH_STD_ARRAYS // Changed guard name slightly
#define EHLERS_MESA_WRAPPER_MQH_STD_ARRAYS

// Define constants if needed - defaults are often passed via Init/Calculate
#define MESA_WRAPPER_MAX_LEN 500
#define MESA_WRAPPER_MAX_COEF 500 // Should be < MAX_LEN
#define MESA_WRAPPER_MAX_PERIOD 500
#define MESA_WRAPPER_HANN_LEN 12

//+------------------------------------------------------------------+
//| EhlersMESAWrapper Class Definition                               |
//+------------------------------------------------------------------+
class EhlersMESAWrapper
  {
private:
   // --- State Variables (Static) ---
   double m_coef2_flat[MESA_WRAPPER_MAX_COEF * (MESA_WRAPPER_HANN_LEN + 1)]; // Flat [NumCoef+1][12+1]
   int               m_numCoefStored;
   int               m_hannLenStored;

   // --- Configuration ---
   int               m_length;
   int               m_numCoef;
   int               m_lowerBound;
   int               m_upperBound;
   bool              m_initialized;

   // --- Internal Arrays (Standard Dynamic Arrays) ---
   double            m_P[];       // Power array
   double            m_bb1[];     // Burg backward array
   double            m_bb2[];     // Burg forward array
   double            m_coef[];    // MESA coefficients
   double            m_coefA[];   // Temp coefficients for Burg/Levinson
   double            m_Hcoef[];   // Hann-windowed coefficients
   double            m_Wave[];    // Spectrum wave power
   double            m_NWave[];   // Normalized wave power
   double            m_temp_bb1[];// Used for Burg update

   // --- Internal Flags for Array Sizing ---
   bool              m_arrays_sized; // Flag to check if arrays have been sized


public:
   // --- Constructor ---
                     EhlersMESAWrapper(void);
   // --- Destructor ---
                    ~EhlersMESAWrapper(void); // Still simple

   // --- Initialization ---
   bool              Init(int length, int numCoef, int lowerBound, int upperBound);

   // --- Calculation Method ---
   double            CalculateDominantCycle(const double &price[], // Price array (non-series)
         const int bar,         // Current bar index 'i'
         const int rates_total);
private:
   // --- Helper ---
   bool              SizeInternalArrays(void); // Renamed from Resize...

  };

// --- Constructor Implementation ---
EhlersMESAWrapper::EhlersMESAWrapper(void) :
   m_numCoefStored(0),
   m_hannLenStored(MESA_WRAPPER_HANN_LEN + 1),
   m_length(0), m_numCoef(0), m_lowerBound(0), m_upperBound(0),
   m_initialized(false),
   m_arrays_sized(false) // Arrays are not sized initially
  {
   ArrayInitialize(m_coef2_flat, 0.0);
// Standard dynamic arrays are initialized empty by default
  }

// --- Destructor Implementation ---
EhlersMESAWrapper::~EhlersMESAWrapper(void) {}

// --- Initialization Implementation ---
bool EhlersMESAWrapper::Init(int length, int numCoef, int lowerBound, int upperBound)
  {
// --- Input Validation (same as before) ---
   if(length < 2)
     {
      Print("Ehlers MESA Wrapper Error: length must be >= 2.");
      m_initialized = false;
      return false;
     }
   if(numCoef < 1)
     {
      Print("Ehlers MESA Wrapper Error: numCoef must be >= 1.");
      m_initialized = false;
      return false;
     }
   if(numCoef >= length)
     {
      Print("Ehlers MESA Wrapper Error: numCoef must be < length.");
      m_initialized = false;
      return false;
     }
   if(lowerBound < 1)
     {
      Print("Ehlers MESA Wrapper Error: lowerBound must be >= 1.");
      m_initialized = false;
      return false;
     }
   if(upperBound <= lowerBound)
     {
      Print("Ehlers MESA Wrapper Error: upperBound must be > lowerBound.");
      m_initialized = false;
      return false;
     }
   if(length >= MESA_WRAPPER_MAX_LEN)
     {
      PrintFormat("Ehlers MESA Wrapper Error: length (%d) >= MAX_LEN (%d).", length, MESA_WRAPPER_MAX_LEN);
      m_initialized = false;
      return false;
     }
   if(numCoef >= MESA_WRAPPER_MAX_COEF)
     {
      PrintFormat("Ehlers MESA Wrapper Error: numCoef (%d) >= MAX_COEF (%d).", numCoef, MESA_WRAPPER_MAX_COEF);
      m_initialized = false;
      return false;
     }
   if(upperBound >= MESA_WRAPPER_MAX_PERIOD)
     {
      PrintFormat("Ehlers MESA Wrapper Error: upperBound (%d) >= MAX_PERIOD (%d).", upperBound, MESA_WRAPPER_MAX_PERIOD);
      m_initialized = false;
      return false;
     }
// --- Store Parameters ---
   m_length = length;
   m_numCoef = numCoef;
   m_lowerBound = lowerBound;
   m_upperBound = upperBound;
// --- Flat Array Size Calc (same as before) ---
   m_numCoefStored = numCoef + 1;
   m_hannLenStored = MESA_WRAPPER_HANN_LEN + 1;
   int required_flat_size = m_numCoefStored * m_hannLenStored;
   if(required_flat_size > ArraySize(m_coef2_flat))
     {
      PrintFormat("Ehlers MESA Wrapper Error: Calculated flat size (%d) exceeds max allowed (%d).", required_flat_size, ArraySize(m_coef2_flat));
      m_initialized = false;
      return false;
     }
   ArrayInitialize(m_coef2_flat, 0.0);

// --- Reset Flags ---
   m_arrays_sized = false; // Force resizing on first Calculate call
   m_initialized = true;
   PrintFormat("EhlersMESAWrapper Initialized. L=%d, NumC=%d, Bounds[%d-%d]", m_length, m_numCoef, m_lowerBound, m_upperBound);
   return true;
  }

// --- Helper to Size Internal Arrays (Using ArrayResize) ---
// This should ideally only be called once or when parameters change significantly
bool EhlersMESAWrapper::SizeInternalArrays(void)
  {
// Check if already sized with current parameters (optional optimization)
// if (m_arrays_sized) return true; // If parameters don't change, maybe skip?

   PrintFormat("SizeInternalArrays: Attempting to size arrays. Target L=%d, NumC=%d, UB=%d",
               m_length, m_numCoef, m_upperBound);

   if(m_length < 1 || m_numCoef < 1 || m_upperBound < 1)
     {
      Print("SizeInternalArrays Error: Invalid target dimensions (length, numCoef, or upperBound < 1).");
      return false;
     }

   bool ok = true;
   string fail_msg = "";

// Use +1 for size because MQL arrays are 0-based but EL often uses 1-based logic
// Using ArrayResize - returns the new size, or -1 on error
   if(ArrayResize(m_P, m_numCoef + 1) < 0)
     {
      PrintFormat("SizeInternalArrays: FAILED m_P to %d", m_numCoef + 1);
      fail_msg +=" P ";
      ok = false;
     }
   if(ArrayResize(m_bb1, m_length + 1) < 0)
     {
      PrintFormat("SizeInternalArrays: FAILED m_bb1 to %d", m_length + 1);
      fail_msg +=" bb1 ";
      ok = false;
     }
   if(ArrayResize(m_bb2, m_length + 1) < 0)
     {
      PrintFormat("SizeInternalArrays: FAILED m_bb2 to %d", m_length + 1);
      fail_msg +=" bb2 ";
      ok = false;
     }
   if(ArrayResize(m_coef, m_numCoef + 1) < 0)
     {
      PrintFormat("SizeInternalArrays: FAILED m_coef to %d", m_numCoef + 1);
      fail_msg +=" coef ";
      ok = false;
     }
   if(ArrayResize(m_coefA, m_numCoef + 1) < 0)
     {
      PrintFormat("SizeInternalArrays: FAILED m_coefA to %d", m_numCoef + 1);
      fail_msg +=" coefA ";
      ok = false;
     }
   if(ArrayResize(m_Hcoef, m_numCoef + 1) < 0)
     {
      PrintFormat("SizeInternalArrays: FAILED m_Hcoef to %d", m_numCoef + 1);
      fail_msg +=" Hcoef ";
      ok = false;
     }
   if(ArrayResize(m_Wave, m_upperBound + 1) < 0)
     {
      PrintFormat("SizeInternalArrays: FAILED m_Wave to %d", m_upperBound + 1);
      fail_msg +=" Wave ";
      ok = false;
     }
   if(ArrayResize(m_NWave, m_upperBound + 1) < 0)
     {
      PrintFormat("SizeInternalArrays: FAILED m_NWave to %d", m_upperBound + 1);
      fail_msg +=" NWave ";
      ok = false;
     }
   if(ArrayResize(m_temp_bb1, m_length + 1)< 0)
     {
      PrintFormat("SizeInternalArrays: FAILED m_temp_bb1 to %d", m_length + 1);
      fail_msg +=" temp_bb1 ";
      ok = false;
     }

   if(!ok)
     {
      Print("Ehlers MESA Wrapper Error: Failed to size internal standard arrays:" + fail_msg);
      m_arrays_sized = false; // Mark as failed
     }
   else
     {
      PrintFormat("SizeInternalArrays: Success. New sizes: P=%d, bb1=%d, bb2=%d, coef=%d, coefA=%d, Hcoef=%d, Wave=%d, NWave=%d, temp_bb1=%d",
                  ArraySize(m_P), ArraySize(m_bb1), ArraySize(m_bb2), ArraySize(m_coef), ArraySize(m_coefA), ArraySize(m_Hcoef), ArraySize(m_Wave), ArraySize(m_NWave), ArraySize(m_temp_bb1));
      m_arrays_sized = true; // Mark as successfully sized
     }
   return ok;
  }

// --- Calculation Method Implementation (Using Standard Arrays) ---
double EhlersMESAWrapper::CalculateDominantCycle(
   const double &price[],
   const int bar,
   const int rates_total)
  {
   bool print_debug = (bar <= 5 || bar >= rates_total - 5 || bar % 100 == 0);

   if(print_debug)
      PrintFormat("MESA Calc [Bar %d] Start", bar);

// --- Initial Checks ---
   if(!m_initialized)
     {
      Print("Ehlers MESA Wrapper Error: Not initialized!");   // Return 0 or default if not init
      return 0.0;
     }
   if(bar < m_length - 1)
     {
      return 0.0;   // Insufficient history, return 0 or default
     }
   if(bar >= rates_total)
     {
      PrintFormat("MESA Calc [Bar %d] Error: bar index >= rates_total (%d)", bar, rates_total);   // Return 0 or default
      return 0.0;
     }

// --- Variable Initialization (Define default_cycle EARLY) ---
   double default_cycle = (m_lowerBound + m_upperBound) / 2.0; // Define default here
   double DominantCycle = default_cycle;
   bool skip_this_bar = false;
   double Pwr = 0.0;
   double Num = 0.0;
   double Denom = 0.0;
   double MaxWave = 0.0;
   double max_val_check = 1e100;

// --- Ensure arrays are sized ---
   if(!m_arrays_sized)
     {
      if(!SizeInternalArrays())
        {
         PrintFormat("MESA Calc [Bar %d] Exit: Failed to size internal arrays.", bar);
         return default_cycle; // Use the defined default_cycle
        }
     }
// *** Add Check AFTER sizing call ***
   if(ArraySize(m_bb1) <= m_length || ArraySize(m_bb2) <= m_length)
     {
      PrintFormat("MESA Calc [Bar %d] Error: bb1/bb2 still not sized correctly! bb1.Size=%d, bb2.Size=%d, m_length=%d", bar, ArraySize(m_bb1), ArraySize(m_bb2), m_length);
      // *** CORRECTED RETURN VALUE ***
      return default_cycle; // Use the MESA default value
     }


// --- 1. Initialize Power (Pwr) (same logic) ---
   Pwr = 0.0;
   for(int tt = 0; tt < m_length; tt++)
     {
      int idx = bar - tt;
      if(idx < 0 || idx >= rates_total)
        {
         PrintFormat("MESA DEBUG [Bar %d]: Power loop ERROR idx out of bounds (idx=%d, tt=%d, rates_total=%d). Skipping.", bar, idx, tt, rates_total);
         skip_this_bar = true;
         break;
        }
      if(!MathIsValidNumber(price[idx]))
        {
         PrintFormat("MESA DEBUG [Bar %d]: Invalid price[idx=%d] in Power loop. Skipping bar.", bar, idx);
         skip_this_bar = true;
         break;
        }
      Pwr += price[idx] * price[idx];
     }
   if(skip_this_bar || m_length <= 0)
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Exit: Skip in Power Init (m_length=%d)", bar, m_length);
      return default_cycle;
     }
   Pwr /= m_length;
   if(!MathIsValidNumber(Pwr))
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Warning: Pwr became invalid (%.4g) after averaging. Setting to 0.", bar, Pwr);
      Pwr = 0.0;
     }


// --- 2. Initialize Burg Arrays (bb1, bb2) (Use standard array access) ---
   double price_bar = price[bar];
   if(!MathIsValidNumber(price_bar))
     {
      PrintFormat("MESA DEBUG [Bar %d]: Invalid price[bar=%d] for bb1[1] init. Skipping.", bar, bar);
      return default_cycle;
     }
   m_bb1[1] = price_bar;

   int idx_end = bar - (m_length - 1);
   if(idx_end < 0 || idx_end >= rates_total)
     {
      PrintFormat("MESA DEBUG [Bar %d]: Invalid price index idx_end=%d for bb2[%d] init. Skipping.", bar, idx_end, m_length - 1);
      return default_cycle;
     }
   double price_end = price[idx_end];
   if(!MathIsValidNumber(price_end))
     {
      PrintFormat("MESA DEBUG [Bar %d]: Invalid price value at idx_end=%d for bb2[%d] init. Skipping.", bar, idx_end, m_length - 1);
      return default_cycle;
     }
   m_bb2[m_length - 1] = price_end;

   for(int tt = 2; tt <= m_length - 1; tt++)
     {
      int idx = bar - (tt - 1);
      if(idx < 0 || idx >= rates_total)
        {
         PrintFormat("MESA DEBUG [Bar %d]: bb1/bb2 loop price index out of bounds (idx=%d, tt=%d). Skipping.", bar, idx, tt);
         skip_this_bar = true;
         break;
        }
      double price_idx = price[idx];
      if(!MathIsValidNumber(price_idx))
        {
         PrintFormat("MESA DEBUG [Bar %d]: Invalid price[idx=%d] = %.4g in bb1/bb2 loop (tt=%d). Setting skip_this_bar.", bar, idx, price_idx, tt);
         skip_this_bar = true;
         continue;
        }
      if(tt >= ArraySize(m_bb1) || (tt - 1) >= ArraySize(m_bb2) || (tt-1) < 1)
        {
         PrintFormat("MESA DEBUG [Bar %d]: bb1/bb2 loop internal array index out of bounds (tt=%d, Size bb1=%d, Size bb2=%d). Skipping.", bar, tt, ArraySize(m_bb1), ArraySize(m_bb2));
         skip_this_bar = true;
         break;
        }
      m_bb1[tt] = price_idx;
      m_bb2[tt - 1] = price_idx;
     }

   if(skip_this_bar)
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Exit: Skip flag set during bb1/bb2 Initialization loop.", bar);
      return default_cycle;
     }

   if(print_debug)   // Debug print after init loop
     {
      string bb1_str = StringFormat("bb1[1]=%.4f", m_bb1[1]);
      if(m_length >= 2)
         bb1_str += StringFormat(", bb1[2]=%.4f", m_bb1[2]);
      if(m_length > 2)
         bb1_str += StringFormat(" ... bb1[%d]=%.4f", m_length-1, m_bb1[m_length-1]);
      PrintFormat("MESA DEBUG [Bar %d] AFTER bb init: %s", bar, bb1_str);
      string bb2_str = StringFormat("bb2[1]=%.4f", m_bb2[1]);
      if(m_length >= 3)
         bb2_str += StringFormat(", bb2[2]=%.4f", m_bb2[2]);
      if(m_length > 2)
         bb2_str += StringFormat(" ... bb2[%d]=%.4f", m_length-1, m_bb2[m_length-1]);
      PrintFormat("MESA DEBUG [Bar %d] AFTER bb init: %s", bar, bb2_str);
      PrintFormat("MESA DEBUG [Bar %d] AFTER bb init: bb1.Size=%d, bb2.Size=%d, m_length=%d", bar, ArraySize(m_bb1), ArraySize(m_bb2), m_length);
     }


// --- 3. Calculate Initial Coefficient (coef[1]) (Use standard array access) ---
   Num = 0.0;
   Denom = 0.0;
   for(int tt = 1; tt <= m_length - 1; tt++)
     {
      if(tt < 1 || tt >= ArraySize(m_bb1) || tt >= ArraySize(m_bb2))
        {
         PrintFormat("MESA Calc [Bar %d] Error: Index tt=%d out of bounds for bb1/bb2 in Coef1 Sum. bb1.Size=%d, bb2.Size=%d", bar, tt, ArraySize(m_bb1), ArraySize(m_bb2));
         skip_this_bar = true;
         break;
        }
      double b1 = m_bb1[tt];
      double b2 = m_bb2[tt];
      if(!MathIsValidNumber(b1) || !MathIsValidNumber(b2))
        {
         if(print_debug)
            PrintFormat("MESA Calc [Bar %d] Warning: Invalid number in Coef1 Sum at tt=%d (b1=%.4g, b2=%.4g). Skipping bar.", bar, tt, b1, b2);
         skip_this_bar = true;
         break;
        }
      double term_num = b1 * b2;
      double term_den = b1 * b1 + b2 * b2;
      if(!MathIsValidNumber(term_num) || !MathIsValidNumber(term_den))
        {
         if(print_debug)
            PrintFormat("MESA Calc [Bar %d] Warning: Invalid intermediate term in Coef1 Sum at tt=%d. Skipping bar.", bar, tt);
         skip_this_bar = true;
         break;
        }
      if(MathAbs(Num + term_num) > max_val_check || MathAbs(Denom + term_den) > max_val_check)
        {
         if(print_debug)
            PrintFormat("MESA Calc [Bar %d] Warning: Potential overflow in Coef1 Sum at tt=%d. Skipping bar.", bar, tt);
         skip_this_bar = true;
         break;
        }
      Num += term_num;
      Denom += term_den;
     }

   if(skip_this_bar)
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Exit: Skip flag set during Coef1 Sum loop.", bar);
      return default_cycle;
     }

   double coef1_temp = 0.0;
   if(MathAbs(Denom) < 1e-10)
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Warning: Coef[1] Denom is near zero (%.4g). Setting coef[1] to 0.", bar, Denom);
      coef1_temp = 0.0;
     }
   else
     {
      coef1_temp = 2.0 * Num / Denom;
      if(!MathIsValidNumber(coef1_temp) || MathAbs(coef1_temp) > 10.0)
        {
         if(print_debug)
            PrintFormat("MESA Calc [Bar %d] Warning: Coef[1] unstable (%.4g), clamping.", bar, coef1_temp);
         coef1_temp = 0.0;
        }
     }
   m_coef[1] = coef1_temp;


// --- 4. Calculate P[1] (Use standard array access) ---
   double coef1_sq = m_coef[1] * m_coef[1];
   if(!MathIsValidNumber(Pwr) || !MathIsValidNumber(coef1_sq))
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Warning: Invalid Pwr or coef1_sq for P[1] calc. Setting P[1] = Pwr.", bar);
      m_P[1] = Pwr;
     }
   else
     {
      m_P[1] = Pwr * (1.0 - coef1_sq);
      if(!MathIsValidNumber(m_P[1]))
        {
         if(print_debug)
            PrintFormat("MESA Calc [Bar %d] Warning: P[1] became invalid after calculation. Setting P[1] = Pwr.", bar);
         m_P[1] = Pwr;
        }
     }
   if(print_debug)
      PrintFormat("MESA Calc [Bar %d]: Pwr=%.5g, Coef[1]=%.4f, P[1]=%.5g", bar, Pwr, m_coef[1], m_P[1]);


// --- 5. Burg Algorithm Loop (Use standard array access) ---
// ... (Logic is identical, just uses standard array [] access) ...
   for(int mm = 2; mm <= m_numCoef; mm++)
     {
      // ... (rest of Burg loop using standard array access) ...
      if(mm-1 >= ArraySize(m_coef) || mm-1 < 1)
        {
         PrintFormat("MESA Calc [Bar %d, mm=%d] Error: Index out of bounds accessing m_coef for coefA. Skipping bar.", bar, mm);
         skip_this_bar = true;
         break;
        }
      for(int tt = 1; tt <= mm - 1; tt++)
        {
         if(tt >= ArraySize(m_coef) || tt < 1 || tt >= ArraySize(m_coefA))
           {
            PrintFormat("MESA Calc [Bar %d, mm=%d] Error: Index tt=%d out of bounds for coefA update. Skipping bar.", bar, mm, tt);
            skip_this_bar = true;
            break;
           }
         m_coefA[tt] = m_coef[tt];
        }
      if(skip_this_bar)
         break;
      double coefA_mm_minus_1 = m_coefA[mm - 1];
      if(!MathIsValidNumber(coefA_mm_minus_1))
        {
         PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: coefA[mm-1] is invalid. Skipping bar.", bar, mm);
         skip_this_bar = true;
         break;
        }
      for(int tt = 1; tt <= m_length - mm; tt++)
        {
         if(tt < 1 || tt >= ArraySize(m_bb1) || tt >= ArraySize(m_bb2) || tt >= ArraySize(m_temp_bb1) || tt+1 >= ArraySize(m_bb1) || tt+1 >= ArraySize(m_bb2))
           {
            PrintFormat("MESA Calc [Bar %d, mm=%d] Error: Index tt=%d out of bounds in Burg update prep. Skipping bar.", bar, mm, tt);
            skip_this_bar = true;
            break;
           }
         double b1_tt = m_bb1[tt];
         double b2_tt = m_bb2[tt];
         if(!MathIsValidNumber(b1_tt) || !MathIsValidNumber(b2_tt))
           {
            PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Invalid b1/b2 at tt=%d in Burg update prep. Skipping bar.", bar, mm, tt);
            skip_this_bar = true;
            break;
           }
         double temp_val = b1_tt - coefA_mm_minus_1 * b2_tt;
         if(!MathIsValidNumber(temp_val))
           {
            PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Invalid temp_bb1 calc at tt=%d. Skipping bar.", bar, mm, tt);
            skip_this_bar = true;
            break;
           }
         m_temp_bb1[tt] = temp_val;
        }
      if(skip_this_bar)
         break;
      for(int tt = 1; tt <= m_length - mm; tt++)
        {
         double b2_tt_plus_1 = m_bb2[tt + 1];
         double b1_tt_plus_1 = m_bb1[tt + 1];
         if(!MathIsValidNumber(b1_tt_plus_1) || !MathIsValidNumber(b2_tt_plus_1))
           {
            PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Invalid b1/b2 at tt+1=%d in Burg bb2 update. Skipping bar.", bar, mm, tt+1);
            skip_this_bar = true;
            break;
           }
         double temp_val = b2_tt_plus_1 - coefA_mm_minus_1 * b1_tt_plus_1;
         if(!MathIsValidNumber(temp_val))
           {
            PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Invalid bb2 calc at tt=%d. Skipping bar.", bar, mm, tt);
            skip_this_bar = true;
            break;
           }
         m_bb2[tt] = temp_val;
        }
      if(skip_this_bar)
         break;
      for(int tt = 1; tt <= m_length - mm; tt++)
        {
         m_bb1[tt] = m_temp_bb1[tt];
        }
      Num = 0.0;
      Denom = 0.0;
      for(int tt = 1; tt <= m_length - mm; tt++)
        {
         double b1 = m_bb1[tt];
         double b2 = m_bb2[tt];
         if(!MathIsValidNumber(b1) || !MathIsValidNumber(b2))
           {
            if(print_debug)
               PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Invalid b1/b2 in Coef[mm] Sum at tt=%d. Skipping bar.", bar, mm, tt);
            skip_this_bar = true;
            break;
           }
         double term_num = b1 * b2;
         double term_den = b1 * b1 + b2 * b2;
         if(!MathIsValidNumber(term_num) || !MathIsValidNumber(term_den))
           {
            if(print_debug)
               PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Invalid intermediate term in Coef[mm] Sum at tt=%d. Skipping bar.", bar, mm, tt);
            skip_this_bar = true;
            break;
           }
         if(MathAbs(Num + term_num) > max_val_check || MathAbs(Denom + term_den) > max_val_check)
           {
            if(print_debug)
               PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Potential overflow in Coef[mm] Sum at tt=%d. Skipping bar.", bar, mm, tt);
            skip_this_bar = true;
            break;
           }
         Num += term_num;
         Denom += term_den;
        }
      if(skip_this_bar)
         break;
      double coef_mm_temp = 0.0;
      if(MathAbs(Denom) < 1e-10)
        {
         if(print_debug)
            PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Coef[mm] Denom is near zero (%.4g). Setting coef[mm] to 0.", bar, mm, Denom);
         coef_mm_temp = 0.0;
        }
      else
        {
         coef_mm_temp = 2.0 * Num / Denom;
         if(!MathIsValidNumber(coef_mm_temp) || MathAbs(coef_mm_temp) > 10.0)
           {
            if(print_debug)
               PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Coef[mm] unstable (%.4g), clamping.", bar, mm, coef_mm_temp);
            coef_mm_temp = 0.0;
           }
        }
      if(mm < 1 || mm >= ArraySize(m_coef))
        {
         PrintFormat("MESA Calc [Bar %d, mm=%d] Error: Index out of bounds updating m_coef. Skipping bar.", bar, mm);
         skip_this_bar = true;
         break;
        }
      m_coef[mm] = coef_mm_temp;
      double p_prev = m_P[mm - 1];
      double coef_mm = m_coef[mm];
      double coef_mm_sq = coef_mm * coef_mm;
      if(!MathIsValidNumber(p_prev) || !MathIsValidNumber(coef_mm_sq))
        {
         if(print_debug)
            PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Invalid p_prev or coef_mm_sq for P[mm] calc. Setting P[mm] = p_prev.", bar, mm);
         m_P[mm] = p_prev;
        }
      else
        {
         m_P[mm] = p_prev * (1.0 - coef_mm_sq);
         if(!MathIsValidNumber(m_P[mm]))
           {
            if(print_debug)
               PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: P[mm] became invalid after calculation. Setting P[mm] = p_prev.", bar, mm);
            m_P[mm] = p_prev;
           }
        }
      if(print_debug && (mm == 2 || mm == m_numCoef || mm == m_numCoef / 2))
        {
         PrintFormat("MESA Calc [Bar %d, mm=%d]: Coef[mm]=%.4f, P[mm]=%.5g", bar, mm, m_coef[mm], m_P[mm]);
        }
      for(int tt = 1; tt <= mm - 1; tt++)
        {
         int idx_mm_minus_tt = mm - tt;
         if(tt < 1 || tt >= ArraySize(m_coefA) || tt >= ArraySize(m_coef) || idx_mm_minus_tt < 1 || idx_mm_minus_tt >= ArraySize(m_coefA))
           {
            PrintFormat("MESA Calc [Bar %d, mm=%d] Error: Index tt=%d or mm-tt=%d out of bounds in Levinson recursion. Skipping bar.", bar, mm, tt, idx_mm_minus_tt);
            skip_this_bar = true;
            break;
           }
         double coefA_tt = m_coefA[tt];
         double coef_mm_val = m_coef[mm];
         double coefA_mm_tt = m_coefA[idx_mm_minus_tt];
         if(!MathIsValidNumber(coefA_tt) || !MathIsValidNumber(coef_mm_val) || !MathIsValidNumber(coefA_mm_tt))
           {
            if(print_debug)
               PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Invalid value in Levinson recursion for coef[%d]. Skipping update.", bar, mm, tt);
            continue;
           }
         double new_coef_tt = coefA_tt - coef_mm_val * coefA_mm_tt;
         if(!MathIsValidNumber(new_coef_tt))
           {
            if(print_debug)
               PrintFormat("MESA Calc [Bar %d, mm=%d] Warning: Invalid result in Levinson recursion for coef[%d]. Skipping update.", bar, mm, tt);
            continue;
           }
         m_coef[tt] = new_coef_tt;
        }
      if(skip_this_bar)
         break;
     } // End mm loop

   if(skip_this_bar)
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Exit: Skip occurred within Burg Loop (mm loop)", bar);
      return default_cycle;
     }

// --- 6. Hann Windowing ---
// ... (Logic is identical, just uses standard m_coef[] access) ...
   for(int tt = 1; tt <= m_numCoef; tt++)
     {
      int base_idx_tt = tt * m_hannLenStored;
      if(base_idx_tt + m_hannLenStored > ArraySize(m_coef2_flat))
        {
         skip_this_bar = true;
         break;
        }
      for(int count_idx = m_hannLenStored - 2; count_idx >= 1; count_idx--)
        {
         int target_idx = base_idx_tt + count_idx + 1;
         int source_idx = base_idx_tt + count_idx;
         if(target_idx < 0 || target_idx >= ArraySize(m_coef2_flat) || source_idx < 0 || source_idx >= ArraySize(m_coef2_flat))
           {
            skip_this_bar = true;
            break;
           }
         m_coef2_flat[target_idx] = m_coef2_flat[source_idx];
        }
      if(skip_this_bar)
         break;
      int newest_flat_idx = base_idx_tt + 1;
      if(tt < 1 || tt >= ArraySize(m_coef) || newest_flat_idx < 0 || newest_flat_idx >= ArraySize(m_coef2_flat))
        {
         skip_this_bar = true;
         break;
        }
      double current_coef_tt = m_coef[tt];
      if(!MathIsValidNumber(current_coef_tt))
         current_coef_tt = 0.0;
      m_coef2_flat[newest_flat_idx] = current_coef_tt;
     }
   if(skip_this_bar)
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Exit: Skip occurred during Hann Coef2 Update", bar);
      return default_cycle;
     }

// --- 7. Calculate Hann Weighted Coefficients (Hcoef) ---
// ... (Logic is identical, just uses standard m_Hcoef[] access) ...
   for(int tt = 1; tt <= m_numCoef; tt++)
     {
      if(tt < 1 || tt >= ArraySize(m_Hcoef))
        {
         skip_this_bar = true;
         break;
        }
      m_Hcoef[tt] = 0.0;
      double hann_sum_weights = 0.0;
      int base_idx_tt = tt * m_hannLenStored;
      if(base_idx_tt + m_hannLenStored > ArraySize(m_coef2_flat))
        {
         skip_this_bar = true;
         break;
        }
      for(int count = 1; count <= MESA_WRAPPER_HANN_LEN; count++)
        {
         double hann_weight = (1.0 - MathCos(2.0 * M_PI * count / (MESA_WRAPPER_HANN_LEN + 1.0)));
         int flat_idx = base_idx_tt + count;
         if(flat_idx < 0 || flat_idx >= ArraySize(m_coef2_flat))
           {
            skip_this_bar = true;
            break;
           }
         double coef2_val = m_coef2_flat[flat_idx];
         if(!MathIsValidNumber(coef2_val))
            coef2_val = 0.0;
         m_Hcoef[tt] += hann_weight * coef2_val;
         hann_sum_weights += hann_weight;
        }
      if(skip_this_bar)
         break;
      if(MathAbs(hann_sum_weights) > 1e-10)
        {
         double hcoef_val = m_Hcoef[tt] / hann_sum_weights;
         if(!MathIsValidNumber(hcoef_val))
            hcoef_val = 0.0;
         m_Hcoef[tt] = hcoef_val;
        }
      else
        {
         m_Hcoef[tt] = 0.0;
        }
     }
   if(skip_this_bar)
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Exit: Skip occurred during HCoef Calculation", bar);
      return default_cycle;
     }
   if(print_debug)
     {
      string hcoef_str = "";
      for(int c=1; c<=m_numCoef; c+=MathMax(1,m_numCoef/5))
        {
         if(c >= 1 && c < ArraySize(m_Hcoef))
           {
            hcoef_str += StringFormat(" H%d=%.3f", c, m_Hcoef[c]);
           }
        }
      PrintFormat("MESA Calc [Bar %d]: HCoeffs Sampled:%s", bar, hcoef_str);
     }

// --- 8. Compute Spectrum (Wave) ---
// ... (Logic is identical, just uses standard array [] access) ...
   if(m_numCoef < 1 || m_numCoef >= ArraySize(m_P))
     {
      PrintFormat("MESA Calc [Bar %d] Error: Index m_numCoef=%d out of bounds for m_P. Skipping bar.", bar, m_numCoef);
      return default_cycle;
     }
   double P_NumCoef = m_P[m_numCoef];
   if(!MathIsValidNumber(P_NumCoef))
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Warning: P[NumCoef=%d] is invalid (%.4g). Using 0 for Spectrum calc.", bar, m_numCoef, P_NumCoef);
      P_NumCoef = 0.0;
     }
   for(int Period = m_lowerBound; Period <= m_upperBound; Period++)
     {
      if(Period <= 0 || Period >= ArraySize(m_Wave))
        {
         PrintFormat("MESA Calc [Bar %d] Warning: Period=%d out of bounds [1..%d] for Wave array. Skipping Period.", bar, Period, ArraySize(m_Wave)-1);
         continue;
        }
      double DenReal = 0.0;
      double DenImag = 0.0;
      for(int mm = 1; mm <= m_numCoef; mm++)
        {
         if(mm < 1 || mm >= ArraySize(m_Hcoef))
           {
            PrintFormat("MESA Calc [Bar %d] Error: Index mm=%d out of bounds for m_Hcoef in Spectrum calc. Skipping bar.", bar, mm);
            skip_this_bar = true;
            break;
           }
         double hcoef_mm = m_Hcoef[mm];
         if(!MathIsValidNumber(hcoef_mm))
           {
            if(print_debug)
               PrintFormat("MESA Calc [Bar %d] Warning: Invalid Hcoef[%d] in Spectrum calc for Period=%d. Ignoring term.", bar, mm, Period);
            continue;
           }
         double angle_rad = 2.0 * M_PI * mm / Period;
         DenReal += hcoef_mm * MathCos(angle_rad);
         DenImag -= hcoef_mm * MathSin(angle_rad);
        }
      if(skip_this_bar)
         break;
      double one_minus_DenReal = 1.0 - DenReal;
      double wave_denom = one_minus_DenReal * one_minus_DenReal + DenImag * DenImag;
      double wave_val = 0.0;
      if(MathAbs(wave_denom) > 1e-10 && MathIsValidNumber(wave_denom))
        {
         wave_val = P_NumCoef / wave_denom;
        }
      else
         if(print_debug)
           {
            PrintFormat("MESA Calc [Bar %d] Warning: Wave[%d] denominator near zero or invalid (%.4g). Setting Wave to 0.", bar, Period, wave_denom);
           }
      if(!MathIsValidNumber(wave_val))
        {
         if(print_debug)
            PrintFormat("MESA Calc [Bar %d] Warning: Wave[%d] became invalid (%.4g). Setting to 0.", bar, Period, wave_val);
         wave_val = 0.0;
        }
      m_Wave[Period] = wave_val;
      if(print_debug && (Period == m_lowerBound || Period == m_upperBound || Period % (MathMax(1,(m_upperBound-m_lowerBound+1)/5)) == 0))
        {
         PrintFormat("MESA Calc [Bar %d]: Wave[%d]=%.5g (DenR=%.3f, DenI=%.3f, Denom=%.4g)", bar, Period, m_Wave[Period], DenReal, DenImag, wave_denom);
        }
     }
   if(skip_this_bar)
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Exit: Skip occurred during Spectrum Calculation", bar);
      return default_cycle;
     }


// --- 9. Normalize Spectrum & Find Peak (Use standard array access) ---
// ... (Logic is identical, just uses standard array [] access) ...
   MaxWave = 0.0;
   if(m_upperBound >= m_lowerBound)
     {
      if(m_lowerBound >= 0 && m_lowerBound < ArraySize(m_Wave) && MathIsValidNumber(m_Wave[m_lowerBound]))
        {
         MaxWave = m_Wave[m_lowerBound];
        }
      else
        {
         for(int Period = m_lowerBound; Period <= m_upperBound; Period++)
           {
            if(Period >= 0 && Period < ArraySize(m_Wave) && MathIsValidNumber(m_Wave[Period]))
              {
               MaxWave = m_Wave[Period];
               break;
              }
           }
         if(MaxWave == 0.0 && print_debug)
            PrintFormat("MESA Calc [Bar %d] Warning: Could not find a valid starting MaxWave. All Waves might be zero or invalid.", bar);
        }
      for(int Period = m_lowerBound + 1; Period <= m_upperBound; Period++)
        {
         if(Period >= 0 && Period < ArraySize(m_Wave) && MathIsValidNumber(m_Wave[Period]) && m_Wave[Period] > MaxWave)
           {
            MaxWave = m_Wave[Period];
           }
        }
     }
   else
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Warning: Invalid bounds (lower=%d, upper=%d) for finding MaxWave.", bar, m_lowerBound, m_upperBound);
      return default_cycle;
     }
   DominantCycle = default_cycle;
   double MaxNWave = -1.0;
   for(int Period = m_lowerBound; Period <= m_upperBound; Period++)
     {
      if(Period < 0 || Period >= ArraySize(m_NWave) || Period >= ArraySize(m_Wave))
        {
         if(print_debug)
            PrintFormat("MESA Calc [Bar %d] Warning: Period=%d out of bounds for NWave/Wave arrays during normalization. Skipping Period.", bar, Period);
         continue;
        }
      double nwave_val = 0.0;
      double current_wave = m_Wave[Period];
      if(MathAbs(MaxWave) > 1e-10 && MathIsValidNumber(current_wave))
        {
         nwave_val = current_wave / MaxWave;
         if(!MathIsValidNumber(nwave_val) || nwave_val < 0)
           {
            if(print_debug)
               PrintFormat("MESA Calc [Bar %d] Warning: NWave[%d] became invalid or negative (%.4g). Setting to 0.", bar, Period, nwave_val);
            nwave_val = 0.0;
           }
        }
      else
        {
         nwave_val = 0.0;
        }
      m_NWave[Period] = nwave_val;
      if(m_NWave[Period] > MaxNWave)
        {
         MaxNWave = m_NWave[Period];
         DominantCycle = Period;
        }
     }


// --- Final Debug Print & Return ---
   if(print_debug)
     {
      PrintFormat("MESA Calc [Bar %d] End: MaxWave=%.5g, MaxNWave=%.3f, DominantCycle=%d", bar, MaxWave, MaxNWave, (int)DominantCycle);
     }
   if(DominantCycle < m_lowerBound || DominantCycle > m_upperBound)
     {
      if(print_debug)
         PrintFormat("MESA Calc [Bar %d] Warning: Final DominantCycle (%d) outside bounds [%d-%d]. Returning default.", bar, (int)DominantCycle, m_lowerBound, m_upperBound);
      return default_cycle;
     }
   return DominantCycle;

  } // End CalculateDominantCycle

// End of class definition moved outside the function accidentally in previous edits, ensure it's here

#endif // EHLERS_MESA_WRAPPER_MQH_STD_ARRAYS
//+------------------------------------------------------------------+
