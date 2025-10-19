#include "portaudio.h"

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
    constexpr double WINDOW_SIZE_FACTOR = 2.0;
    constexpr double ABSOLUTE_THRESHOLD = 0.1;
    constexpr int SAMPLE_RATE = 44100;
    constexpr int HOP_SIZE = 256;
    constexpr double SUGGESTED_INPUT_LATENCY = 0;
    constexpr double SUGGESTED_OUTPUT_LATENCY = 0.025;

    // Derived constants
    constexpr int MAX_LAG = (SAMPLE_RATE + MIN_FREQUENCY - 1) / MIN_FREQUENCY;
    constexpr int WINDOW_SIZE = MAX_LAG * WINDOW_SIZE_FACTOR;
    constexpr int YIN_REQUIRED_SIZE = WINDOW_SIZE + MAX_LAG;
    constexpr int BUFFER_SIZE = YIN_REQUIRED_SIZE + HOP_SIZE;

    // Global constants
    constexpr std::array<const char*, 12> NOTES = {" A", "A#", " B", " C", "C#", " D", "D#", " E", " F", "F#", " G", "G#"};

}

// Modulo operation for doubles since fmod is technically a remainder operator
double mod(double a, double b)
{
    return a - b * std::floor(a / b);
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

    // Estimate pitch
    double F_0 = estimateF_0(buffer.data());

    // Take prior pitch if current pitch is out of bounds
    static double priorF_0 = 440;
    if (F_0 > 95 && F_0 < 800)
    {
        priorF_0 = F_0;
    }
    else
    {
        F_0 = priorF_0;
    }

    // Determine pitch period from frequency
    double pitchPeriod = SAMPLE_RATE / F_0;

    // Print pitch
    std::cout << "Pitch: " << frequencyToNote(F_0) << ", i.e. " << F_0 << " Hz" << std::endl;

    // Clear the output so we can layer in tones
    for (int i = 0; i < HOP_SIZE; i++)
    {
        out[i] = 0;
    }

    // Determine target pitch and pitch ratio
    double targetF_0 = 220;
    double pitchRatio = targetF_0 / F_0;

    // Resample the hop according to the pitch ratio
    static double sourceIndex = 0;
    for (int targetIndex = 0; targetIndex < HOP_SIZE; targetIndex++)
    {
        // Linearly interpolate values lying between samples
        int x1 = static_cast<int>(sourceIndex);
        int x2 = x1 + 1;
        double y1 = buffer[x1];
        double y2 = buffer[x2];
        out[targetIndex] += ((y2 - y1) * (sourceIndex - x1) + y1) / 2;

        // Increment sourceIndex by pitchRatio, backing up one period if we run out of data
        sourceIndex += pitchRatio;
        if (sourceIndex >= buffer.size())
        {
            sourceIndex -= pitchPeriod;
        }
    }

    double targetF_0_2 = 277.18;
    double pitchRatio2 = targetF_0_2 / F_0;

    // Resample the hop according to the pitch ratio
    static double sourceIndex2 = 0;
    for (int targetIndex = 0; targetIndex < HOP_SIZE; targetIndex++)
    {
        // Linearly interpolate values lying between samples
        int x1 = static_cast<int>(sourceIndex2);
        int x2 = x1 + 1;
        double y1 = buffer[x1];
        double y2 = buffer[x2];
        out[targetIndex] += ((y2 - y1) * (sourceIndex2 - x1) + y1) / 2;

        // Increment sourceIndex by pitchRatio, backing up one period if we run out of data
        sourceIndex2 += pitchRatio2;
        if (sourceIndex2 >= buffer.size())
        {
            sourceIndex2 -= pitchPeriod;
        }
    }

    // Consume this hop from the buffer
    buffer.erase(buffer.begin(), buffer.begin() + HOP_SIZE);

    // Account for the buffer shift before the next hop
    sourceIndex -= HOP_SIZE;
    sourceIndex2 -= HOP_SIZE;

    // Shift the sourceindex back to the start (mod pitchPeriod) to ensure a continuous waveform
    sourceIndex = mod(sourceIndex, pitchPeriod);
    sourceIndex2 = mod(sourceIndex2, pitchPeriod);

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