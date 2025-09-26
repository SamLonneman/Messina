#include <iostream>
#include "AudioFile.h"


long long df(int16_t* samples, int window_size, int lag)
{
    long long sum = 0;
    for (int j = 0; j < window_size; j++)
    {
        sum += (samples[j] - samples[j + lag]) * (samples[j] - samples[j + lag]);
    }
    return sum;
}

int main()
{
    std::cout << "Hello, World!" << std::endl;

    AudioFile<int16_t> audioFile;
    audioFile.load("../audio/a440.wav");
    audioFile.printSummary();

    int channel = 0;
    int numSamples = audioFile.getNumSamplesPerChannel();
    int sampleRate = audioFile.getSampleRate();
    for (int i = 0; i < 20; i++)
    {
        int16_t currentSample = audioFile.samples[channel][i];
        std::cout << i << ": " << currentSample << std::endl;
    }


    int candidate = 0;
    int windowSize = 20000;
    for (int lag = 0; lag < 500; lag++)
    {
        std::cout << "lag: " << lag << ", df: " << df(audioFile.samples[channel].data() + candidate, windowSize, lag) << std::endl;
    }

    int period = 100;
    std::cout << sampleRate/period << std::endl;
}