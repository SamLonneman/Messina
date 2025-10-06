#include "portaudio.h"
#include "CircularBuffer.h"

#include <array>
#include <climits>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>


namespace {
    // User-configurable constants
    constexpr int MIN_FREQUENCY = 100;
    constexpr float WINDOW_SIZE_SCALAR = 1;
    constexpr double ABSOLUTE_THRESHOLD = 0.1;
    constexpr int SAMPLE_RATE = 44100;
    constexpr int HOP_SIZE = 256;
    constexpr double SUGGESTED_INPUT_LATENCY = 0;
    constexpr double SUGGESTED_OUTPUT_LATENCY = 0.02;

    // Derived constants
    constexpr int MAX_LAG = (SAMPLE_RATE + MIN_FREQUENCY - 1) / MIN_FREQUENCY;
    constexpr int WINDOW_SIZE = MAX_LAG * WINDOW_SIZE_SCALAR;
    constexpr int YIN_REQUIRED_SIZE = WINDOW_SIZE + MAX_LAG;
    constexpr int BUFFER_SIZE = YIN_REQUIRED_SIZE + HOP_SIZE;

    // Global constants
    constexpr std::array<const char*, 12> NOTES = {" A", "A#", " B", " C", "C#", " D", "D#", " E", " F", "F#", " G", "G#"};

    // Global state
    double priorF_0 = 0;

}

// Given three points (x values must be separated by 1), calculate x value of vertex of interpolated parabola via magic
double parabolicInterpolation(double x1, double y1, double x2, double y2, double x3, double y3)
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
    int differenceInCents = std::round(std::log2(F_0 / 440.0) * 1200);
    int differenceInSemitones = (differenceInCents + 50) / 100;
    std::string note = NOTES[((differenceInSemitones % 12) + 12) % 12];
    int octave = ((differenceInSemitones + 9) / 12) + 4;
    int errorInCents = ((((differenceInCents + 50) % 100) + 100) % 100) - 50;
    return note + std::to_string(octave) + " (" + (errorInCents >= 0 ? "+" : "") + std::to_string(errorInCents) + " cents)";
}

static int audioCallback(const void* input, void* output, unsigned long frameCount, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void* userData)
{
    // Cast input, output, and buffer to appropriate types
    const int16_t* in = static_cast<const int16_t*>(input);
    int16_t* out = static_cast<int16_t*>(output);
    CircularBuffer<int16_t>& circularBuffer = *static_cast<CircularBuffer<int16_t>*>(userData);

    // If no input or output, just continue
    if (!in || !out)
    {
        return paContinue;
    }

    // Write into circular buffer
    circularBuffer.write(in, frameCount);
    
    // Keep collecting data until we have enough to run YIN
    if (circularBuffer.size() < YIN_REQUIRED_SIZE)
    {
        // Until then, just send input directly to output
        std::memcpy(out, in, frameCount);
        return paContinue;
    }

    // Read circular buffer into yinBuffer
    std::array<int16_t, YIN_REQUIRED_SIZE> yinBuffer;
    circularBuffer.read(&yinBuffer[0], YIN_REQUIRED_SIZE);
    
    // Estimate pitch
    double F_0 = estimateF_0(&yinBuffer[0]);
    
    // Smooth pitch output
    if (F_0 > 95 && F_0 < 800) {
        priorF_0 = F_0;
    }
    else
    {
        F_0 = priorF_0;
    }
    
    // Print pitch
    std::cout << "Pitch: " << frequencyToNote(F_0) << std::endl;
    
    // Consume from the circular buffer now that we have run YIN
    circularBuffer.read(out, frameCount);
    circularBuffer.consume(frameCount);

    // Continue processing
    return paContinue;
    
}

// Entry point
int main()
{
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

    // Initialize a circular buffer to persist between callbacks
    CircularBuffer<int16_t> circularBuffer(BUFFER_SIZE);

    // Open and start stream
    Pa_OpenStream(&stream, &inputParams, &outputParams, SAMPLE_RATE, HOP_SIZE, paClipOff, audioCallback, &circularBuffer);
    Pa_StartStream(stream);

    // Pause the main thread until user input
    std::cin.get();

    // Stop and close stream, terminate PortAudio
    Pa_StopStream(stream);
    Pa_CloseStream(stream);
    Pa_Terminate();

}