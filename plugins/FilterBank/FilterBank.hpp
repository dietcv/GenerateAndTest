#pragma once
#include "Utils.hpp"

class FilterBank : public SCUnit {
public:
    FilterBank();
    ~FilterBank() = default;
private:
    void next(int nSamples);
   
    // Constants
    const float m_sampleRate;
    static constexpr int NUM_BANDS = 24;

    // Core processing
    Utils::FilterBank<NUM_BANDS> filterBank;
   
    // Cache for SlopeSignal state
    float freqPast, spreadPast, warpPast, resonancePast;
   
    enum InputParams {
        Input,      // Audio input
        Freq,       // Base frequency
        Spread,     // Frequency spread
        Warp,       // Frequency warping
        Resonance   // Filter resonance (0-1)
    };
   
    enum Outputs {
        Out
    };
};