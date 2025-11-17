#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"
#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
using namespace std;

// Instrucciones de uso:
// audio WAV
// audio MONO
// audio de 16 bits

// Lectura de audio
vector <double> cargar_normalizar_wav (const char* filename) {
    
    drwav wav;
    
    // Excepción de error al abrir
    if (!drwav_init_file(&wav, filename, NULL))
        throw runtime_error("No se abre audio WAV");

    // Excepción de audio MONO
    if (wav.channels != 1) {
        drwav_uninit(&wav);
        throw runtime_error("El audio solo permite de 1 canal (MONO)");
    }

    // Lectura de muestras de audios WAV
    vector<int16_t> samples(wav.totalPCMFrameCount);
    drwav_read_pcm_frames_s16(&wav, wav.totalPCMFrameCount, samples.data());
    drwav_uninit(&wav);

    // Normalización
    vector <double> normalizado;
    normalizado.reserve(samples.size());

    for (int16_t s : samples)
        normalizado.push_back(s / 32768.0);
    
    return normalizado;
}


// FFT



// Filtrado de frecuencias cardiacas
void filtrarFrecuencias(vector<complex<double>>& fft, double fs) { // fs = frecuencia de muestreo
    int n = fft.size();
    double minFreq = 0.5;
    double maxFreq = 3.5;

    for (int i = 0; i <= n/2; i++) {
        double freq = (i * fs) / n;

        if (freq < minFreq || freq > maxFreq) {
            fft[i] = {0.0, 0.0};

            if (i != 0 && i != n/2)
                fft[n - i] = {0.0, 0.0};
        }
    }
}


// IFFT



// Extracción de BPM



// Detección de anomalias


// Main