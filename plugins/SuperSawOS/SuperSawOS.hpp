#pragma once
#include "Utils.hpp"
#include "VariableOversampling.hpp"

class SuperSawOS : public SCUnit {
public:
    SuperSawOS();
    ~SuperSawOS() = default;
private:
    void next_aa(int nSamples);

    // Constants
    const float m_sampleRate;
    
    // Core processing
    Utils::SuperSawOsc m_superSaw;
    Utils::OnePoleSlope m_pitchTrackingHPF;
    VariableOversampling<4> m_oversampling;
    
    enum InputParams {
        Freq,           // Base frequency
        Mix,            // SuperSaw mix (0-1): 0=center only, 1=full spread
        Detune,         // Detune amount (0-1)
        Oversample      // Oversampling factor (0=1x, 1=2x, 2=4x, 3=8x, 4=16x)
    };
    
    enum Outputs {
        Out
    };
};