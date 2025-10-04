#include "portaudio.h"

#include <climits>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>


namespace {
    // Yin algorithm parameters
    const int minFrequency = 100;
    const float windowSizeScalar = 2;
    const double absoluteThreshold = 0.1;

    // Global variables
    const int sampleRate = 44100;
    int maxLag;
    int windowSize;

}

// Given three points (x values must be separated by 1), calculate x value of vertex of interpolated parabola via magic
double parabolicInterpolation(double x1, double y1, double x2, double y2, double x3, double y3)
{
    return x2 + (y1 - y3) / (2.0 * (y1 - 2.0 * y2 + y3));
}

// Difference function
long long df(int16_t* samples, int lag)
{
    long long sum = 0;
    for (int j = 1; j <= windowSize; j++)
    {
        int diff = samples[j] - samples[j + lag];
        sum += diff * diff;
    }
    return sum;
}

// Given a waveform at some time t, estimate the fundamental frequency using the YIN algorithm
double estimateF_0(int16_t* samples)
{
    // For each lag value up to maxLag,
    int optimalLag = -1;
    double minCMNDF = std::numeric_limits<double>::max();
    long long sumDF = 0;
    bool thresholdReached = false;
    for (int currentLag = 0; currentLag < maxLag; currentLag++)
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
        if (currentCMNDF < absoluteThreshold)
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
    return sampleRate / periodInSamples;

}

static int audioCallback(const void* input, void* output, unsigned long frameCount, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void* userData)
{
    // Cast input and userData to appropriate types
    const int16_t* in = (const int16_t*)input;
    std::vector<int16_t>& buffer = *(std::vector<int16_t>*)userData;

    // Copy input to buffer
    buffer.assign(in, in + frameCount);

    // Run YIN
    double pitch = estimateF_0(buffer.data());
    std::cout << "Pitch: " << pitch << " Hz" << std::endl;

    // Continue processing
    return paContinue;

}

// Entry point
int main()
{
    // Determine maxLag and window size
    maxLag = (sampleRate + minFrequency - 1) / minFrequency;
    windowSize = maxLag * windowSizeScalar;

    // Real time pitch detection
    Pa_Initialize();
    PaStream* stream;
    std::vector<uint16_t> audioBuffer(2048);
    Pa_OpenDefaultStream(&stream, 1, 0, paInt16, sampleRate, 1024, audioCallback, &audioBuffer);
    Pa_StartStream(stream);

    std::cin.get();

    Pa_StopStream(stream);
    Pa_CloseStream(stream);
    Pa_Terminate();

}