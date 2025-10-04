# Messina

## Installation and setup
```bash
git clone https://github.com/SamLonneman/Messina.git
mkdir external
cd external
git clone https://github.com/adamstark/AudioFile.git
git clone https://github.com/PortAudio/portaudio.git
cd ..
mkdir build
cd build
cmake ..
```

## Build and run
```bash
cmake --build .
./messina
```