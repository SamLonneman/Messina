#include <iostream>
#include <vector>
#include <climits>
#include "AudioFile.h"


long long df(int16_t* samples, int windowSize, int lag)
{
    long long sum = 0;
    for (int j = 1; j <= windowSize; j++)
    {
        sum += pow(samples[j] - samples[j + lag], 2);
    }
    return sum;
}

int main()
{
    // Load audio file and print summary
    AudioFile<int16_t> audioFile;
    audioFile.load("../audio/a440.wav");
    audioFile.printSummary();

    // Get stats
    int numSamples = audioFile.getNumSamplesPerChannel();
    int sampleRate = audioFile.getSampleRate();
    
    // Determine interesting lags and window size
    int minFrequency = 135;                                      // Lowest sung note (Hz)
    int maxFrequency = 1400;                                     // Highest sung note (Hz)
    int minLag = sampleRate / maxFrequency;                      // Smallest period to check (samples, rounded down)
    int maxLag = (sampleRate + minFrequency - 1) / minFrequency; // Largest period to check (samples, rounded up)
    int windowSize = maxLag * 2.5;                               // Rule of thumb is to use window size at least twice maxLag.
    
    // Display lag and window size info
    std::cout << "minLag: " << minLag << " samples, i.e. " << (double)minLag * 1000 / sampleRate << " ms, i.e. " << (double)sampleRate / minLag << " Hz" << std::endl;
    std::cout << "maxLag: " << maxLag << " samples, i.e. " << (double)maxLag * 1000 / sampleRate << " ms, i.e. " << (double)sampleRate / maxLag << " Hz" << std::endl;
    std::cout << "windowSize: " << windowSize << " samples, i.e. " << (double)windowSize * 1000 / sampleRate << "ms, i.e. 2.5x maxLag" << std::endl;
    
    // Estimate pitch at t = 0
    int candidate = 0;
    int channel = 0;
    int optimalLag = -1;
    long long minDF = LLONG_MAX;
    for (int currentLag = minLag; currentLag < maxLag; currentLag++)
    {
        long long currentDF = df(audioFile.samples[channel].data() + candidate, windowSize, currentLag);
        std::cout << "lag: " << currentLag << ", df: " << currentDF << std::endl;
        if (currentDF < minDF)
        {
            minDF = currentDF;
            optimalLag = currentLag;
        }
    }

    // Print result
    double F_0 = (double)sampleRate / optimalLag;
    std::cout << "Estimated F_0 at t = " << (double)candidate * 1000 / sampleRate << " ms: " << F_0 << std::endl;

}