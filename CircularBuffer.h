#pragma once
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <vector>

template <typename T>
class CircularBuffer
{   
private:
    // Internal working variables
    std::vector<T> buffer_;
    size_t capacity_;
    size_t writePtr_;
    size_t readPtr_;
    size_t size_;

public:
    // Construct a circular buffer with a given capacity
    CircularBuffer(size_t capacity) :
        buffer_(capacity, 0),
        capacity_(capacity),
        writePtr_(0),
        readPtr_(0),
        size_(0)
    { }

    // Write data to the circular buffer
    void write(const T *data, size_t numSamples)
    {
        if (numSamples > capacity_ - size_)
        {
            throw std::runtime_error("Not enough space available in buffer to write!");
        }
        size_t firstChunk = std::min(numSamples, capacity_ - writePtr_);
        size_t secondChunk = numSamples - firstChunk;
        std::memcpy(&buffer_[writePtr_], data, firstChunk * sizeof(T));
        std::memcpy(&buffer_[0], data + firstChunk, secondChunk * sizeof(T));
        writePtr_ = (writePtr_ + numSamples) % capacity_;
        size_ += numSamples;
    }

    // Read data from the circular buffer, with optional offset
    void read(T *output, size_t numSamples, size_t offset=0) const
    {
        if (offset + numSamples > size_)
        {
            throw std::runtime_error("Not enough data in buffer to read!");
        }
        size_t start = (readPtr_ + offset) % capacity_;
        size_t firstChunk = std::min(numSamples, capacity_ - start);
        size_t secondChunk = numSamples - firstChunk;
        std::memcpy(output, &buffer_[start], firstChunk * sizeof(T));
        std::memcpy(output + firstChunk, &buffer_[0], secondChunk * sizeof(T));
    }

    // Consume data from the circular buffer
    void consume(size_t numSamples)
    {
        if (numSamples > size_)
        {
            throw std::runtime_error("Not enough data in buffer to consume!");
        }
        readPtr_ = (readPtr_ + numSamples) % capacity_;
        size_ -= numSamples;
    }

    // Get the current size of the buffer
    size_t size() const
    {
        return size_;
    }

};