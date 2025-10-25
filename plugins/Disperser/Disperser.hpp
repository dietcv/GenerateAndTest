#pragma once
#include "Utils.hpp"

class Disperser : public SCUnit {
public:
    Disperser();
    ~Disperser() = default;
private:
    void next(int nSamples);
   
    // Constants
    const float m_sampleRate;
    static constexpr int NUM_ALLPASSES = 8;
   
    // Core processing
    Utils::Disperser<NUM_ALLPASSES> disperser;
    Utils::OnePoleHz m_dcBlocker;

    // Feedback state
    float m_feedbackState{0.0f};
    
    // Cache for SlopeSignal state
    float freqPast, resonancePast, mixPast, feedbackPast;
   
    enum InputParams {
        Input,      // Audio input
        Freq,       // Allpass frequency
        Resonance,  // Filter resonance (0-1)
        Mix,        // Dry/wet mix
        Feedback    // Feedback amount (0-0.999)
    };
   
    enum Outputs {
        Out
    };
};