# Messina

## Installation and setup
```bash
git clone https://github.com/SamLonneman/Messina.git
cd Messina
git clone https://github.com/PortAudio/portaudio.git external/portaudio
git clone https://github.com/thestk/rtmidi.git external/rtmidi
cmake -B build
```

## Build and run (development)
```bash
cmake --build build
./build/messina
```