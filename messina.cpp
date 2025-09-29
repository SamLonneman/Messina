#include <iostream>
#include <iomanip>
#include <chrono>
#include <climits>
#include <vector>

#include "AudioFile.h"


void printVector(const std::vector<double>& vec)
{
    std::cout << "[";
    for (double value : vec)
    {
        std::cout << "  " << std::setw(10) << value;
    }
    std::cout << "  ]" << std::endl;
}

void printMatrix(const std::vector<std::vector<double>>& mat)
{
    for (const std::vector<double>& row : mat)
    {
        printVector(row);
    }
}

std::vector<double> gaussianElimination(std::vector<std::vector<double>>& mat)
{
    // Get matrix dimensions
    int m = mat.size();
    int n = mat[0].size();

    // Forward elimination: for each column...
    int currentRow = 0;
    for (int currentColumn = 0; currentColumn < n; currentColumn++)
    {
        // Break if currentRow is out of bounds
        if (currentRow >= m)
        {
            break;
        }

        // Beginning at the current row, find the value in this column with the largest absolute value
        int largestAbsValue = 0;
        int largestAbsValueRow = currentRow;
        for (int i = currentRow; i < m; i++)
        {
            int currentAbsValue = std::abs(mat[i][currentColumn]);
            if (currentAbsValue > largestAbsValue)
            {
                largestAbsValue = currentAbsValue;
                largestAbsValueRow = i;
            }
        }

        // Swap current row with the found row if applicable
        if (largestAbsValueRow != currentRow)
        {
            std::swap(mat[currentRow], mat[largestAbsValueRow]);
        }

        // If the current value is still 0, then this is not a pivot column, so proceed to the next column.
        if (mat[currentRow][currentColumn] == 0)
        {
            continue;
        }

        // Add a multiple of the current row to each other row to make their current column value 0.
        for (int i = currentRow + 1; i < m; i++)
        {
            double multiplier = -mat[i][currentColumn] / mat[currentRow][currentColumn];
            for (int j = currentColumn; j < n; j++)
            {
                mat[i][j] += mat[currentRow][j] * multiplier;
            }
        }

        // Advance to next row
        currentRow++;
    }

    // Verify that the rank is sufficient to yield a unique solution.
    if (currentRow < n - 1)
    {
        throw std::runtime_error("Back substitution failed: system does not have a unique solution.");
    }

    // Perform back substitution to calculate solution
    std::vector<double> solution(n - 1);
    for (int i = m - 1; i >= 0; i--)
    {
        solution[i] = mat[i][n-1];
        for (int j = n - 2; j > i; j--)
        {
            solution[i] -= mat[i][j] * solution[j];
        }
        solution[i] /= mat[i][i];
    }

    // Return solution vector
    return solution;

}

// Autocorrelation function
long long acf(int16_t* samples, int windowSize, int lag)
{
    long long sum = 0;
    for (int j = 1; j <= windowSize; j++)
    {
        sum += samples[j] * samples[j + lag];
    }
    return sum;
}

// Difference function v1 (intuitive form)
long long dfv1(int16_t* samples, int windowSize, int lag)
{
    long long sum = 0;
    for (int j = 1; j <= windowSize; j++)
    {
        int diff = samples[j] - samples[j + lag];
        sum += diff * diff;
    }
    return sum;
}

// Difference function v2 (in terms of ACF)
long long dfv2(int16_t* samples, int windowSize, int lag)
{
    return acf(samples, windowSize, 0) + acf(samples + lag, windowSize, 0) - 2 * acf(samples, windowSize, lag);
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

    // Determine lags of interest and window size
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
    int16_t* samples = audioFile.samples[channel].data();
    for (int currentLag = minLag; currentLag < maxLag; currentLag++)
    {
        long long currentDF = dfv1(samples + candidate, windowSize, currentLag);
        if (currentDF < minDF)
        {
            minDF = currentDF;
            optimalLag = currentLag;
        }
    }

    // Parabolic interpolation
    long long priorDF = dfv1(samples + candidate, windowSize, optimalLag - 1);
    long long optimalDF = dfv1(samples + candidate, windowSize, optimalLag);
    long long nextDF = dfv1(samples + candidate, windowSize, optimalLag + 1);
    std::vector<std::vector<double>> systemOfEquations =
    {
        { (double)((optimalLag - 1) * (optimalLag - 1)), (double)(optimalLag - 1), 1, (double)priorDF   },
        { (double)((optimalLag    ) * (optimalLag    )), (double)(optimalLag    ), 1, (double)optimalDF },
        { (double)((optimalLag + 1) * (optimalLag + 1)), (double)(optimalLag + 1), 1, (double)nextDF    }
    };
    std::vector<double> solution = gaussianElimination(systemOfEquations);
    double interpolatedOptimalLag = -solution[1] / (2 * solution[0]);

    // From interpolated optimal lag, calculate final F_0 estimate
    double F_0 = sampleRate / interpolatedOptimalLag;
    std::cout << "\nEstimated F_0 at t = " << (double)candidate * 1000 / sampleRate << " ms: " << F_0 << " Hz" << std::endl;

}