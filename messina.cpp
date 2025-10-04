#include "AudioFile.h"

#include <climits>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>


namespace {

    // Hyperparameters
    const int minFrequency = 135; // Lowest sung note (Hz)
    const int maxFrequency = 1400; // Highest sung note (Hz)
    const float windowSizeScalar = 2;
    const double absoluteThreshold = 0.1;

    // Global variables
    int numSamples;
    int sampleRate;
    int minLag;
    int maxLag;
    int windowSize;

}

// Pretty-print a vector
void printVector(const std::vector<double>& vec)
{
    std::cout << "[";
    for (double value : vec)
    {
        std::cout << "  " << std::setw(10) << value;
    }
    std::cout << "  ]" << std::endl;
}

// Pretty-print a matrix
void printMatrix(const std::vector<std::vector<double>>& mat)
{
    for (const std::vector<double>& row : mat)
    {
        printVector(row);
    }
}

// Convert a period expressed in samples into a period expressed in Milliseconds
double samplesToMs(double numSamples)
{
    return numSamples * 1000 / sampleRate;
}

// Convert a period expressed in samples into a frequency expressed in Hz
double samplesToHz(double numSamples)
{
    return sampleRate / numSamples;
}

// Given an augmented matrix representing a system of linear equations, returns a solution vector using Gaussian elimination with back substitution.
std::vector<double> solveSystemOfEquations(std::vector<std::vector<double>>& mat)
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

// Given three points, returns the x value of the vertex of the interpolated parabola
double parabolicInterpolation(double x1, double y1, double x2, double y2, double x3, double y3)
{
    // If the three points are colinear, the parabola is trivial.
    if (y1 == y2 == y3)
    {
        return x2;
    }

    // Otherwise, prepare and solve a system of linear equations for the parabolic coefficients
    std::vector<std::vector<double>> systemOfEquations =
    {
        { std::pow(x1, 2), x1, 1, y1 },
        { std::pow(x2, 2), x2, 1, y2 },
        { std::pow(x3, 2), x3, 1, y3 }
    };
    std::vector<double> coefficients = solveSystemOfEquations(systemOfEquations);

    // Return x position of vertex by calculating the value of x for which the derivative of the parabola is 0.
    return -coefficients[1] / (2 * coefficients[0]);

}

// Autocorrelation function
long long acf(int16_t* samples, int lag)
{
    long long sum = 0;
    for (int j = 1; j <= windowSize; j++)
    {
        sum += samples[j] * samples[j + lag];
    }
    return sum;
}

// Difference function v1 (intuitive form)
long long dfv1(int16_t* samples, int lag)
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
long long dfv2(int16_t* samples, int lag)
{
    return acf(samples, 0) + acf(samples + lag, 0) - 2 * acf(samples, lag);
}

// Given a waveform at some time t, estimate the fundamental frequency using the YIN algorithm
double estimateF_0(int16_t* samples)
{
    // For each lag value up to maxLag,
    int optimalLag = -1;
    double minCMNDF = std::numeric_limits<double>::max();
    long long sumDF = 0;
    for (int currentLag = 0; currentLag < maxLag; currentLag++)
    {
        
        // Calculate the cumulative mean normalized difference function of that lag at time t
        long long currentDF = dfv1(samples, currentLag);
        sumDF += currentDF;
        double currentCMNDF = currentLag ? (double)currentDF * currentLag / sumDF : 1;

        // Find the lag value which minimizes the CMNDF
        if (currentCMNDF < minCMNDF)
        {
            minCMNDF = currentCMNDF;
            optimalLag = currentLag;

            // If a lag value is found with a CMNDF under the absolute threshold, select it immediately
            if (currentCMNDF < absoluteThreshold)
            {
                break;
            }
        }
    }

    // Use parabolic interpolation of the neighbors of the selected lag value to capture the true period if it lies between samples
    long long priorDF = dfv1(samples, optimalLag - 1);
    long long optimalDF = dfv1(samples, optimalLag);
    long long nextDF = dfv1(samples, optimalLag + 1);
    double periodInSamples = parabolicInterpolation(optimalLag - 1, priorDF, optimalLag, optimalDF, optimalLag + 1, nextDF);

    // Calculate and return final fundamental frequency estimate from interpolated period
    return samplesToHz(periodInSamples);

}

// Entry point
int main()
{
    // Load audio file and print summary
    AudioFile<int16_t> audioFile;
    audioFile.load("../audio/a440_vocal.wav");
    int16_t* samples = audioFile.samples[0].data();
    
    // Get audio file properties
    numSamples = audioFile.getNumSamplesPerChannel();
    sampleRate = audioFile.getSampleRate();
    
    // Determine lags of interest and window size
    minLag = sampleRate / maxFrequency;                      // Smallest period to check (samples, rounded down)
    maxLag = (sampleRate + minFrequency - 1) / minFrequency; // Largest period to check (samples, rounded up)
    windowSize = maxLag * windowSizeScalar;                  // Rule of thumb is to use window size at least twice maxLag.
    
    // Print relevant debug info
    std::cout << std::endl;
    audioFile.printSummary();
    std::cout << "minLag: " << minLag << " samples, i.e. " << samplesToMs(minLag) << " ms, i.e. " << samplesToHz(minLag) << " Hz" << std::endl;
    std::cout << "maxLag: " << maxLag << " samples, i.e. " << samplesToMs(maxLag) << " ms, i.e. " << samplesToHz(maxLag) << " Hz" << std::endl;
    std::cout << "windowSize: " << windowSize << " samples, i.e. " << samplesToMs(windowSize) << " ms, i.e. 2.5x maxLag" << std::endl;
    std::cout << "|======================================|" << std::endl << std::endl;

    // Estimate pitch throughout file, every 1000 samples
    double sum = 0;
    int n = 0;
    for (int i = 0; i < numSamples - windowSize; i += 1000)
    {
        double F_0 = estimateF_0(samples + i);
        sum += F_0;
        n++;
        std::cout << "Estimated F_0 at sample " << i << " (t = " << samplesToMs(i) << " ms): " << F_0 << " Hz" << std::endl;
    }
    std::cout << "Average F_0: " << sum / n << std::endl;

}