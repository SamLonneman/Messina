#include "portaudio.h"
#include "RtMidi.h"

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>


namespace {

    // User-configurable constants
    constexpr int MIN_FREQUENCY = 100;
    constexpr int MAX_FREQUENCY = 800;
    constexpr double WINDOW_SIZE_FACTOR = 2.0;
    constexpr double ABSOLUTE_THRESHOLD = 0.1;
    constexpr int SAMPLE_RATE = 44100;
    constexpr int HOP_SIZE = 256;
    constexpr double SUGGESTED_INPUT_LATENCY = 0;
    constexpr double SUGGESTED_OUTPUT_LATENCY = 0.025;
    constexpr bool DO_LOW_PASS_FILTER = false;

    // Derived constants
    constexpr int MAX_LAG = (SAMPLE_RATE + MIN_FREQUENCY - 1) / MIN_FREQUENCY;
    constexpr int WINDOW_SIZE = MAX_LAG * WINDOW_SIZE_FACTOR;
    constexpr int YIN_REQUIRED_SIZE = WINDOW_SIZE + MAX_LAG;
    constexpr int BUFFER_SIZE = YIN_REQUIRED_SIZE + HOP_SIZE;
    constexpr int YIN_STRIDE = 2;

    // Global constants
    constexpr std::array<const char*, 12> NOTES = {" A", "A#", " B", " C", "C#", " D", "D#", " E", " F", "F#", " G", "G#"};

    // Scales
    constexpr std::array<int, 8> MAJOR = {0, 2, 4, 5, 7, 9, 11, 12};
    constexpr std::array<int, 8> MINOR = {0, 2, 3, 5, 7, 8, 10, 12};

    // Global variables
    RtMidiIn midiIn;

}

// Modulo operation for doubles since fmod is technically a remainder operator
inline double mod(double a, double b)
{
    return a - b * std::floor(a / b);
}

// Given three points with x values separated by 1, calculate the x value of the vertex of the interpolated parabola
inline double parabolicInterpolation(double x1, double y1, double x2, double y2, double x3, double y3)
{
    return x2 + (y1 - y3) / (2.0 * (y1 - 2.0 * y2 + y3));
}

// Difference function
long long df(const int16_t* samples, int lag)
{
    long long sum = 0;
    for (int j = 1; j <= WINDOW_SIZE; j++)
    {
        int diff = samples[j] - samples[j + lag];
        sum += diff * diff;
    }
    return sum;
}

// Given a waveform at some time t, estimate the fundamental frequency using the YIN algorithm
double estimateF_0(const int16_t* samples)
{
    // For each lag value up to MAX_LAG,
    int optimalLag = -1;
    double minCMNDF = std::numeric_limits<double>::max();
    long long sumDF = 0;
    bool thresholdReached = false;
    for (int currentLag = 0; currentLag < MAX_LAG; currentLag++)
    {

        // Calculate the cumulative mean normalized difference function of that lag at time t
        long long currentDF = df(samples, currentLag);
        sumDF += currentDF;
        double currentCMNDF = currentLag ? (double)currentDF * currentLag / sumDF : 1;

        // Keep track of the lag value which minimizes the CMNDF
        if (currentCMNDF < minCMNDF)
        {
            minCMNDF = currentCMNDF;
            optimalLag = currentLag;
        }

        // If CMNDF falls below and then exceeds the absoluteThreshold, immediately select the best lag found thusfar
        if (currentCMNDF < ABSOLUTE_THRESHOLD)
        {
            thresholdReached = true;
        }
        else if (thresholdReached)
        {
            break;
        }
    }

    // Use parabolic interpolation of the neighbors of the selected lag value to capture the true period if it lies between samples
    long long priorDF = df(samples, optimalLag - 1);
    long long optimalDF = df(samples, optimalLag);
    long long nextDF = df(samples, optimalLag + 1);
    double periodInSamples = parabolicInterpolation(optimalLag - 1, priorDF, optimalLag, optimalDF, optimalLag + 1, nextDF);

    // Calculate and return final fundamental frequency estimate from interpolated period
    return SAMPLE_RATE / periodInSamples;

}

// Given a fundamental frequency F_0 in Hz, return the note name
std::string frequencyToNote(double F_0)
{
    int centsFromA440 = std::round(std::log2(F_0 / 440.0) * 1200);
    int semitonesFromA440 = (centsFromA440 + 50) / 100;
    std::string note = NOTES[((semitonesFromA440 % 12) + 12) % 12];
    int octave = ((semitonesFromA440 + 9) / 12) + 4;
    int errorInCents = ((((centsFromA440 + 50) % 100) + 100) % 100) - 50;
    return note + std::to_string(octave) + " (" + (errorInCents >= 0 ? "+" : "") + std::to_string(errorInCents) + " cents)";
}

// Quantize the given frequency to nearest note of the 12 tone equal tempered system
inline double quantizeChromatic(double F_0)
{
    return std::exp2(std::round(std::log2(F_0 / 440) * 12) / 12) * 440;
}

// Quantize the given frequency to nearest note in a given key (tonic and scale)
double quantizeToScale(double F_0, int tonic, const std::array<int, 8>& scale)
{
    double semitonesFromA440 = std::log2(F_0 / 440) * 12;
    double semitonesFromTonicMod12 = mod(semitonesFromA440 - tonic, 12);
    double smallestDiff = scale[0] - semitonesFromTonicMod12;
    double minAbsDiff = std::abs(smallestDiff);
    for (int note : scale)
    {
        double diff = note - semitonesFromTonicMod12;
        double absDiff = std::abs(diff);
        if (absDiff < minAbsDiff)
        {
            smallestDiff = diff;
            minAbsDiff = absDiff;
        }
    }
    return std::exp2((semitonesFromA440 + smallestDiff) / 12) * 440;
}

// Structure used to represent a voice to resample
struct Voice
{
    double frequency;
    int volume;
    double sourceIndex;
    bool sustained;
};

// Callback function which takes input from the mic and returns output to the speakers
static int audioCallback(const void* input, void* output, unsigned long frameCount, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void* userData)
{
    // Cast input, output, and buffer to appropriate types
    const int16_t* in = static_cast<const int16_t*>(input);
    int16_t* out = static_cast<int16_t*>(output);
    std::vector<int16_t>& buffer = *static_cast<std::vector<int16_t>*>(userData);

    // If no input or output, just continue
    if (!in || !out)
    {
        return paContinue;
    }

    // Write into buffer
    buffer.insert(buffer.end(), in, in + HOP_SIZE);

    // Keep collecting data until we have enough to run YIN
    if (buffer.size() < YIN_REQUIRED_SIZE)
    {
        return paContinue;
    }

    // Estimate pitch once every YIN_STRIDE frames
    double F_0 = 0;
    static int yinCounter = 0;
    if (++yinCounter == YIN_STRIDE)
    {
        F_0 = estimateF_0(buffer.data());
        yinCounter = 0;
    }

    // Take prior pitch if current pitch is out of bounds
    static double priorF_0 = 440;
    if (F_0 > MIN_FREQUENCY && F_0 < MAX_FREQUENCY)
    {
        priorF_0 = F_0;
    }
    else
    {
        F_0 = priorF_0;
    }

    // Determine pitch period from frequency
    double pitchPeriod = SAMPLE_RATE / F_0;

    // Quantize frequency to nearest semitone
    double quantizedF_0 = quantizeToScale(F_0, 4, MAJOR);

    // Persistent list of active voices; first is the main autotuned passthrough.
    static std::vector<Voice> voices = {{quantizedF_0, 50, 0}};

    // Update the frequency of the main voice
    voices[0].frequency = quantizedF_0;

    // Read all messages currently in MIDI queue
    std::vector<unsigned char> message;
    while (midiIn.getMessage(&message), !message.empty())
    {
        // If the message is about the sustain pedal
        static bool sustain = false;
        if (message[0] == 176)
        {
            // If it is a sustain release
            if (!(sustain = message[2]))
            {
                // Erase all voices flagged as sustained
                for (int i = voices.size() - 1; i >= 1; i--)
                {
                    if (voices[i].sustained)
                    {
                        voices.erase(voices.begin() + i);
                    }
                }
            }
            continue;
        }

        // Otherwise, we have a standard note, so calculate frequency and volume
        double frequency = std::exp2((message[1] - 69) / 12.0) * 440;
        int volume = message[2];

        // If the volume is non-zero, we have the start of a note
        if (volume)
        {
            // Iterate through list of voices looking for thise note
            bool alreadyHere = false;
            for (int i = 1; i < voices.size(); i++)
            {
                // If already in voices, update the volume and remove sustained flag
                if (voices[i].frequency == frequency)
                {
                    voices[i].volume = volume;
                    voices[i].sustained = false;
                    alreadyHere = true;
                    break;
                }
            }

            // If not found in voices, add it
            if (!alreadyHere)
            {
                voices.push_back({frequency, volume, 0, false});
            }
            continue;
        }
        
        // Otherwise, we have the end of a note. Find it in the list
        for (int i = 1; i < voices.size(); i++)
        {
            if (voices[i].frequency == frequency)
            {
                // If sustain is on, flag it as sustained, otherwise erase it
                if (sustain)
                {
                    voices[i].sustained = true;
                }
                else
                {
                    voices.erase(voices.begin() + i);
                }
                break;
            }
        }
    }

    // Precompute uniform gain (lower all levels based on the number of voices)
    double uniformGain = 1 / std::sqrt(voices.size());
    
    // Clear the output so we can layer in voices
    memset(out, 0, HOP_SIZE * sizeof(int16_t));

    // Resample each target frequency into out
    for (Voice& voice : voices)
    {
        // Get target pitch ratio
        double pitchRatio = voice.frequency / F_0;

        // Dynamic contrast: give more weight to notes with more volume
        double dynamicGain = std::pow(voice.volume / 80.0 + 0.3, 5);

        // EQ Adjustment: give more weight to notes with lower frequencies
        double equalizerGain = std::pow(voice.frequency / 440, -0.5);

        // Final gain: all aspects of gain combined
        double gain = uniformGain * dynamicGain * equalizerGain;

        // Perform resampling
        for (int targetIndex = 0; targetIndex < HOP_SIZE; targetIndex++)
        {
            // Linearly interpolate values lying between samples
            int x1 = static_cast<int>(voice.sourceIndex);
            int x2 = x1 + 1;
            double y1 = buffer[x1];
            double y2 = buffer[x2];
            out[targetIndex] += ((y2 - y1) * (voice.sourceIndex - x1) + y1) * gain;

            // Increment sourceIndex by pitchRatio, backing up one period if we run out of data
            voice.sourceIndex += pitchRatio;
            if (voice.sourceIndex >= buffer.size())
            {
                voice.sourceIndex -= pitchPeriod;
            }
        }
    }

    // Optionally apply a low pass filter
    if (DO_LOW_PASS_FILTER)
    {
        for (int i = 0; i < HOP_SIZE; i++)
        {
            static float prev = 0;
            out[i] = prev = prev + 0.15f * (out[i] - prev);
        }
    }

    // Consume this hop from the buffer
    buffer.erase(buffer.begin(), buffer.begin() + HOP_SIZE);

    // Reset source indices carefully
    for (Voice& voice : voices)
    {
        // Account for the buffer shift before the next hop
        voice.sourceIndex -= HOP_SIZE;

        // Shift the sourceindex back to the start (mod pitchPeriod to ensure a continuous waveform)
        voice.sourceIndex = mod(voice.sourceIndex, pitchPeriod);   
    }

    // Continue processing
    return paContinue;

}

// Entry point
int main()
{
    // Initialize RtMidi
    midiIn.openPort(0);

    // Initialize PortAudio
    Pa_Initialize();
    PaStream* stream;

    // Set up input parameters
    PaStreamParameters inputParams;
    inputParams.device = Pa_GetDefaultInputDevice();
    inputParams.channelCount = 1;
    inputParams.sampleFormat = paInt16;
    inputParams.suggestedLatency = SUGGESTED_INPUT_LATENCY;
    inputParams.hostApiSpecificStreamInfo = nullptr;

    // Set up output parameters
    PaStreamParameters outputParams;
    outputParams.device = Pa_GetDefaultOutputDevice();
    outputParams.channelCount = 1;
    outputParams.sampleFormat = paInt16;
    outputParams.suggestedLatency = SUGGESTED_OUTPUT_LATENCY;
    outputParams.hostApiSpecificStreamInfo = nullptr;

    // Initialize a buffer to persist between callbacks
    std::vector<int16_t> buffer(BUFFER_SIZE);

    // Open and start stream
    Pa_OpenStream(&stream, &inputParams, &outputParams, SAMPLE_RATE, HOP_SIZE, paClipOff, audioCallback, &buffer);
    Pa_StartStream(stream);

    // Pause the main thread until user input
    std::cin.get();

    // Stop and close stream, terminate PortAudio
    Pa_StopStream(stream);
    Pa_CloseStream(stream);
    Pa_Terminate();

}