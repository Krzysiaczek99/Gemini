// EhlersAutocorrWrapper.mqh
// Wrapper for John Ehlers' Autocorrelation Dominant Cycle calculation

#ifndef EHLERS_AUTOCORR_WRAPPER_MQH
#define EHLERS_AUTOCORR_WRAPPER_MQH

#include <Math\Stat\Math.mqh>
#include <Arrays/ArrayDouble.mqh> // For ArrayInitialize

// --- Constants from original function ---
#define AUTOCORR_MAX_LAG 48
#define AUTOCORR_MIN_PERIOD 10
#define AUTOCORR_CORR_START_LAG 3
#define AUTOCORR_ARRAY_SIZE (AUTOCORR_MAX_LAG + 1) // Size 49 (Indices 0-48)

//+------------------------------------------------------------------+
//| EhlersAutocorrWrapper Class                                      |
//+------------------------------------------------------------------+
class EhlersAutocorrWrapper
{
private:
   // --- Configuration ---
   int  m_AvgLength;
   bool m_initialized;

   // --- State Variables ---
   double m_MaxPwr_prev;                  // Stores Max Power from the previous bar
   double m_R1_state[AUTOCORR_ARRAY_SIZE]; // Stores R[Period, 1] (smoothed power) from the previous bar

public:
   // --- Constructor ---
                     EhlersAutocorrWrapper(void);
   // --- Destructor ---
                    ~EhlersAutocorrWrapper(void); // Simple destructor

   // --- Initialization ---
   bool              Init(int avgLength);

   // --- Calculation Method ---
   double            CalculateDominantCycle(const double &price[], // Price array (non-series)
                                            const int bar,         // Current bar index 'i'
                                            const int rates_total); // Total bars available
};

// --- Constructor Implementation ---
EhlersAutocorrWrapper::EhlersAutocorrWrapper(void) :
   m_AvgLength(0),
   m_initialized(false),
   m_MaxPwr_prev(0.0)
{
   ArrayInitialize(m_R1_state, 0.0); // Initialize state array
}

// --- Destructor Implementation ---
EhlersAutocorrWrapper::~EhlersAutocorrWrapper(void)
{
   // Nothing dynamic to delete
}

// --- Initialization Implementation ---
bool EhlersAutocorrWrapper::Init(int avgLength)
{
   m_AvgLength = avgLength; // Store averaging length (0 means use Lag)

   // Reset state variables
   m_MaxPwr_prev = 0.0;
   ArrayInitialize(m_R1_state, 0.0);

   m_initialized = true;
   Print("EhlersAutocorrWrapper Initialized. AvgLength=", m_AvgLength);
   return true;
}


//+------------------------------------------------------------------+
//| CalculateDominantCycle Method                                    |
//+------------------------------------------------------------------+
double EhlersAutocorrWrapper::CalculateDominantCycle(
   const double &price[],
   const int bar,
   const int rates_total)
{
    const double DefaultCycle = (AUTOCORR_MIN_PERIOD + AUTOCORR_MAX_LAG) / 2.0; // Default return value

    // --- Check Initialization & History ---
    if (!m_initialized) { Print("EhlersAutocorrWrapper Error: Not initialized!"); return DefaultCycle; }
    if (bar < AUTOCORR_MAX_LAG) { return DefaultCycle; } // Not enough history for max lag
    if (bar >= rates_total) { Print("AutocorrWrapper Error [Bar ",bar,"]: Bar index out of bounds."); return DefaultCycle;}


    // --- Local arrays for intermediate calculations ---
    double Corr[AUTOCORR_ARRAY_SIZE];       ArrayInitialize(Corr, 0.0);
    double CosinePart[AUTOCORR_ARRAY_SIZE]; ArrayInitialize(CosinePart, 0.0);
    double SinePart[AUTOCORR_ARRAY_SIZE];   ArrayInitialize(SinePart, 0.0);
    double SqSum[AUTOCORR_ARRAY_SIZE];      ArrayInitialize(SqSum, 0.0);
    double PwrNorm[AUTOCORR_ARRAY_SIZE];    ArrayInitialize(PwrNorm, 0.0);
    double R1_current[AUTOCORR_ARRAY_SIZE]; ArrayInitialize(R1_current, 0.0); // Store current bar's smoothed power

    bool history_ok = true; // Flag for price data validity

    // 1. Autocorrelation Loop
    for (int Lag = 0; Lag <= AUTOCORR_MAX_LAG; Lag++) {
        // Determine averaging length M
        int M = m_AvgLength;
        if (m_AvgLength <= 0) M = Lag; // If AvgLength is 0, use Lag
        M = MathMin(M, bar - Lag + 1); // Ensure M doesn't exceed available data
        if (M <= 1) continue;

        double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0;
        for (int count = 0; count < M; count++) {
            int idxX = bar - count;
            int idxY = bar - Lag - count;
            if (idxX < 0 || idxY < 0 || idxX >= rates_total || idxY >= rates_total) { history_ok = false; break; } // Check index validity

            double X = price[idxX]; double Y = price[idxY];
            if (!MathIsValidNumber(X) || !MathIsValidNumber(Y)) { history_ok = false; break; } // Check value validity

            Sx += X; Sy += Y; Sxx += X * X; Sxy += X * Y; Syy += Y * Y;
        }
        if (!history_ok || M <= 1) { Corr[Lag] = 0.0; continue; } // Skip if history bad or M too small

        double VarX = M * Sxx - Sx * Sx;
        double VarY = M * Syy - Sy * Sy;
        double denom = VarX * VarY;
        Corr[Lag] = (denom > 1e-12) ? (M * Sxy - Sx * Sy) / MathSqrt(denom) : 0.0;
    }
    if (!history_ok) { return DefaultCycle; } // Exit if price data was bad

    // 2. Fourier Transform Loop & Raw Power (SqSum)
    for (int Period = AUTOCORR_MIN_PERIOD; Period <= AUTOCORR_MAX_LAG; Period++) {
        if (Period >= AUTOCORR_ARRAY_SIZE) continue; // Safety check
        for (int N = AUTOCORR_CORR_START_LAG; N <= AUTOCORR_MAX_LAG; N++) {
            if (N >= 0 && N < AUTOCORR_ARRAY_SIZE) { // Check Corr index
                CosinePart[Period] += Corr[N] * MathCos(2 * M_PI * N / Period);
                SinePart[Period] += Corr[N] * MathSin(2 * M_PI * N / Period);
            }
        }
        SqSum[Period] = CosinePart[Period] * CosinePart[Period] + SinePart[Period] * SinePart[Period];
    }

    // 3. Power Smoothing (R array calculation) & Store Current Smoothed Power
    for (int Period = AUTOCORR_MIN_PERIOD; Period <= AUTOCORR_MAX_LAG; Period++) {
        if (Period < 0 || Period >= AUTOCORR_ARRAY_SIZE) continue; // Safety check index

        double r1_prev = m_R1_state[Period]; // Get R[Period, 1] from the PREVIOUS bar's state
        if(!MathIsValidNumber(r1_prev)) r1_prev = 0.0; // Handle potential NaN from state

        double sqsum_sq = SqSum[Period] * SqSum[Period];
        if(!MathIsValidNumber(sqsum_sq)) sqsum_sq = 0.0;

        // R[Period, 1] = .2*SqSum[Period]*SqSum[Period] + .8*R[Period, 2];
        double r1_new = 0.2 * sqsum_sq + 0.8 * r1_prev; // Calculate R[Period, 1] for CURRENT bar
        R1_current[Period] = r1_new; // Store locally for normalization & CoG THIS bar
    }

    // 4. Find Maximum Smoothed Power Level
    double CurrentMaxPwr = 0.991 * m_MaxPwr_prev; // Apply decay to previous max
    for (int Period = AUTOCORR_MIN_PERIOD; Period <= AUTOCORR_MAX_LAG; Period++) {
        if (Period >= 0 && Period < AUTOCORR_ARRAY_SIZE) {
             double r1_curr = R1_current[Period]; // Use the value calculated THIS bar
             if (MathIsValidNumber(r1_curr) && r1_curr > CurrentMaxPwr)
                 CurrentMaxPwr = r1_curr;
        }
    }
    // Update state AFTER using it for calculations on this bar
    m_MaxPwr_prev = CurrentMaxPwr;

    // 5. Normalize Smoothed Power
    if (MathAbs(CurrentMaxPwr) > 1e-12) {
        for (int Period = AUTOCORR_MIN_PERIOD; Period <= AUTOCORR_MAX_LAG; Period++) {
            if (Period >= 0 && Period < AUTOCORR_ARRAY_SIZE) {
                 double r1_curr = R1_current[Period]; // Use current bar's smoothed power
                 if (MathIsValidNumber(r1_curr)) PwrNorm[Period] = r1_curr / CurrentMaxPwr;
                 else PwrNorm[Period] = 0.0;
            }
        }
    } else {
        ArrayInitialize(PwrNorm, 0.0); // Set all to zero if MaxPower is zero
    }

    // 6. Compute the dominant cycle using the CG of the Normalized Smoothed Spectrum
    double Spx = 0., Sp = 0.;
    double DominantCycle = DefaultCycle; // Start with default
    for (int Period = AUTOCORR_MIN_PERIOD; Period <= AUTOCORR_MAX_LAG; Period++) {
        if (Period >= 0 && Period < AUTOCORR_ARRAY_SIZE) {
            if (PwrNorm[Period] >= 0.5) { // <<< EL Threshold is 0.5
                Spx += (double)Period * PwrNorm[Period];
                Sp += PwrNorm[Period];
            }
        }
    }

    if (MathAbs(Sp) > 1e-12) DominantCycle = Spx / Sp; // Calculate CoG

    // 7. Final Clamp
    DominantCycle = fmax(AUTOCORR_MIN_PERIOD, fmin(DominantCycle, AUTOCORR_MAX_LAG));

    // 8. *** Update the R1 state AFTER all calculations for this bar are done ***
    for (int Period = AUTOCORR_MIN_PERIOD; Period <= AUTOCORR_MAX_LAG; Period++) {
         if (Period >= 0 && Period < AUTOCORR_ARRAY_SIZE) {
             m_R1_state[Period] = R1_current[Period]; // Store this bar's R1 for the *next* bar
         }
    }

    return DominantCycle;
}


#endif // EHLERS_AUTOCORR_WRAPPER_MQH