// EhlersGriffithsWrapper.mqh
// Wrapper for John Ehlers' Griffiths Dominant Cycle (TASC Jan 2005) (v1.1 - Pwr Fix)

#ifndef EHLERS_GRIFFITHS_WRAPPER_MQH
#define EHLERS_GRIFFITHS_WRAPPER_MQH

#include <Math\Stat\Math.mqh>
#include <Arrays/Array.mqh> // For ArrayInitialize

// --- Constants ---
#define GRIFFITHS_WRAPPER_MAX_LEN 200 // Max allowed Length parameter

//+------------------------------------------------------------------+
//| EhlersGriffithsWrapper Class                                     |
//+------------------------------------------------------------------+
class EhlersGriffithsWrapper
{
private:
   // --- Configuration ---
   int    m_Length;
   int    m_LowerBound;
   int    m_UpperBound;
   bool   m_initialized;
   double m_Mu;

   // --- State Variables ---
   double m_XX[];
   double m_coef[];
   double m_Pwr[];           // <<< Renamed from Pwr in previous static version
   double m_CyclePrev;

public:
   // --- Constructor ---
                     EhlersGriffithsWrapper(void);
   // --- Destructor ---
                    ~EhlersGriffithsWrapper(void);

   // --- Initialization ---
   bool              Init(int length, int lowerBound, int upperBound);

   // --- Calculation Method ---
   double            CalculateDominantCycle(const double &signal[],
                                            const int bar,
                                            const int rates_total);
private:
    bool             ResizeArrays();

};

// --- Constructor Implementation ---
EhlersGriffithsWrapper::EhlersGriffithsWrapper(void) :
   m_Length(0), m_LowerBound(0), m_UpperBound(0), m_initialized(false), m_Mu(0.0), m_CyclePrev(0.0) {}

// --- Destructor Implementation ---
EhlersGriffithsWrapper::~EhlersGriffithsWrapper(void) {}


// --- Helper to Resize Dynamic Arrays ---
bool EhlersGriffithsWrapper::ResizeArrays()
{
    if (m_Length < 2 || m_UpperBound < m_LowerBound || m_UpperBound <= 0) { Print("GriffithsWrapper ResizeArrays Error: Invalid Length or Bounds."); return false; }
    bool ok = true; string fail_msg = "";
    if(ArrayResize(m_XX, m_Length + 1) < 0)    { fail_msg += " XX "; ok = false; }
    if(ArrayResize(m_coef, m_Length + 1) < 0)  { fail_msg += " coef "; ok = false; }
    if(ArrayResize(m_Pwr, m_UpperBound + 1) < 0) { fail_msg += " Pwr "; ok = false; } // Check m_Pwr resize
    if (!ok) { Print("GriffithsWrapper ResizeArrays Error: Failed to resize:" + fail_msg); }
    else {
        ArrayInitialize(m_XX, 0.0); ArrayInitialize(m_coef, 0.0); ArrayInitialize(m_Pwr, 0.0); // Initialize m_Pwr
        PrintFormat("GriffithsWrapper ResizeArrays Success: XX=%d, coef=%d, Pwr=%d", ArraySize(m_XX), ArraySize(m_coef), ArraySize(m_Pwr));
    }
    return ok;
}


// --- Initialization Implementation ---
bool EhlersGriffithsWrapper::Init(int length, int lowerBound, int upperBound)
{
   if (length < 2 || length > GRIFFITHS_WRAPPER_MAX_LEN) { Print("EhlersGriffithsWrapper Init Error: Length (", length, ") must be >= 2 and <= ", GRIFFITHS_WRAPPER_MAX_LEN); m_initialized = false; return false; }
   if (lowerBound < 1) { Print("EhlersGriffithsWrapper Init Error: LowerBound (", lowerBound, ") must be >= 1."); m_initialized = false; return false; }
   if (upperBound <= lowerBound) { Print("EhlersGriffithsWrapper Init Error: UpperBound (", upperBound, ") must be > LowerBound (", lowerBound, ")."); m_initialized = false; return false; }
   // Remove pre-check for Pwr size, ResizeArrays handles it
   // if (upperBound >= ArraySize(m_Pwr) && ArrayResize(m_Pwr, upperBound + 1) < 0) { ... }

   m_Length = length; m_LowerBound = lowerBound; m_UpperBound = upperBound;
   m_Mu = (m_Length > 0) ? (1.0 / m_Length) : 0.01;

   if (!ResizeArrays()) { m_initialized = false; return false; }

   m_CyclePrev = (m_LowerBound + m_UpperBound) / 2.0;
   m_initialized = true;
   Print("EhlersGriffithsWrapper Initialized. L=", m_Length, ", Bounds[", m_LowerBound, "-", m_UpperBound, "], Mu=", m_Mu);
   return true;
}


//+------------------------------------------------------------------+
//| CalculateDominantCycle Method                                    |
//+------------------------------------------------------------------+
double EhlersGriffithsWrapper::CalculateDominantCycle(
   const double &signal[],
   const int bar,
   const int rates_total)
{
    const double DefaultCycle = (m_LowerBound + m_UpperBound) / 2.0;
    bool print_griff_debug = (bar % 100 == 0 || bar >= rates_total - 10);

    // --- Checks ---
    if (!m_initialized) { Print("EhlersGriffithsWrapper Error: Not initialized!"); return DefaultCycle; }
    if (bar < m_Length) { return m_CyclePrev; }
    if (bar >= rates_total || bar - m_Length < 0) { Print("GriffithsWrapper Error [Bar ",bar,"]: Bar index or history index out of bounds."); return m_CyclePrev; }
    if (ArraySize(m_XX) <= m_Length || ArraySize(m_coef) <= m_Length || ArraySize(m_Pwr) <= m_UpperBound) { Print("GriffithsWrapper Error [Bar ",bar,"]: Internal arrays not sized correctly!"); return m_CyclePrev; }

    bool calculation_error = false;

    // 1. Fill XX array
    for(int Count=1; Count <= m_Length; Count++) {
       int signal_index = bar - m_Length + Count;
       if(signal_index < 0 || signal_index >= rates_total) { calculation_error=true; break; }
       if(!MathIsValidNumber(signal[signal_index])) {
             if(print_griff_debug) PrintFormat("Griffiths Warning [Bar %d]: Invalid signal value at index %d. Using 0.", bar, signal_index);
             m_XX[Count] = 0.0;
             // calculation_error = true; break; // Option: Treat as error
       } else { m_XX[Count] = signal[signal_index]; }
    }
    if(calculation_error) { if(print_griff_debug)Print("Griffiths Error [Bar %d]: History index error filling XX.", bar); return m_CyclePrev; }

    // 2. Calculate XBar
    double XBar=0.;
    for(int Count=1; Count <= m_Length; Count++) {
       int xx_index = m_Length - Count + 1;
       XBar += m_XX[xx_index] * m_coef[Count];
    }
    if (!MathIsValidNumber(XBar)) { PrintFormat("Griffiths Warning [Bar %d]: XBar became invalid. Using 0.", bar); XBar = 0.0; }

    // 3. Calculate Error and Update Coefficients
    double error_term = m_XX[m_Length] - XBar;
    if (!MathIsValidNumber(error_term)) { if(print_griff_debug) PrintFormat("Griffiths Warning [Bar %d]: error_term became invalid. Skipping coef update.", bar); }
    else {
        for(int Count=1; Count <= m_Length; Count++) {
           int xx_index = m_Length - Count + 1;
           double coef_update = m_Mu * error_term * m_XX[xx_index];
           if(MathIsValidNumber(coef_update) && MathIsValidNumber(m_coef[Count])) {
                m_coef[Count] = m_coef[Count] + coef_update;
           }
        }
    }

    // 4. Calculate Power Spectrum using *updated* coefficients
    double Real_pwr, Imag_pwr, Denom_pwr;
    for(int Period=m_LowerBound; Period <= m_UpperBound; Period++) {
       if(Period <= 0 || Period >= ArraySize(m_Pwr)) continue; // Use ArraySize(m_Pwr)
       Real_pwr=0.; Imag_pwr=0.;
       for(int Count=1; Count <= m_Length; Count++) {
          if(!MathIsValidNumber(m_coef[Count])) { calculation_error=true; break;} // Check coef validity
          double angle_rad=2.*M_PI*Count/Period;
          Real_pwr += m_coef[Count]*MathCos(angle_rad);
          Imag_pwr += m_coef[Count]*MathSin(angle_rad);
       }
       if(calculation_error) break;

       if(!MathIsValidNumber(Real_pwr) || !MathIsValidNumber(Imag_pwr)) {
           m_Pwr[Period] = 0.0; continue; // Use m_Pwr
       }

       Denom_pwr=(1.-Real_pwr)*(1.-Real_pwr)+Imag_pwr*Imag_pwr;
       m_Pwr[Period]=(MathAbs(Denom_pwr)>1e-12)?(0.1/Denom_pwr):0.; // Use m_Pwr
       if (!MathIsValidNumber(m_Pwr[Period])) m_Pwr[Period] = 0.0; // Use m_Pwr
    }
    if(calculation_error) { if(print_griff_debug)Print("Griffiths Error [Bar %d]: Calculation error during spectrum calc.", bar); return m_CyclePrev; }


    // 5. Find Dominant Cycle from Spectrum
    double MaxPwr=-1.;
    double Cycle=DefaultCycle; // Start with default
    for(int Period=m_LowerBound; Period <= m_UpperBound; Period++) {
       if(Period >= ArraySize(m_Pwr)) continue; // Use ArraySize(m_Pwr)
       if(m_Pwr[Period]>MaxPwr) { // Use m_Pwr
          MaxPwr=m_Pwr[Period]; // Use m_Pwr
          Cycle=Period;
       }
    }

    // 6. Smooth/Clamp Cycle
    if(Cycle > m_CyclePrev+2.) Cycle=m_CyclePrev+2.;
    if(Cycle < m_CyclePrev-2.) Cycle=m_CyclePrev-2.;
    Cycle=fmax(m_LowerBound,fmin(Cycle,m_UpperBound));

    // 7. Update State for *NEXT* bar
    m_CyclePrev = Cycle;

    // --- Debug Print ---
    if(print_griff_debug) PrintFormat("Griffiths [Bar %d]: Period=%.1f", bar, Cycle);

    // 8. Return cycle calculated for *THIS* bar
    return Cycle;
}


#endif // EHLERS_GRIFFITHS_WRAPPER_MQH