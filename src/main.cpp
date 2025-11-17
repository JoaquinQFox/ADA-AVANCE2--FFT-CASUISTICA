#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"
#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// Lectura de audio
vector <double> cargar_normalizar_wav (const char* filename) {
    
    drwav wav;

    // Lectura de muestras de audios WAV
    vector<int16_t> samples(wav.totalPCMFrameCount);
    drwav_read_pcm_frames_s16(&wav, wav.totalPCMFrameCount, samples.data());
    drwav_uninit(&wav);
}


// FFT



// Filtrado de frecuencias cardiacas



// IFFT



// Extracción de BPM



// Detección de anomalias


// Main