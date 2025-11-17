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
// solo lee desde la carpeta de audios

// Lectura de audio
vector <double> cargar_normalizar_wav (const char* filename) {
    
    drwav wav;
    string ruta = string("audios/") + filename;

    // Excepción de error al abrir
    if (!drwav_init_file(&wav, ruta.c_str(), NULL))
        throw runtime_error("No se encuentra el audio WAV");

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
vector<complex<double>> fft(const vector<complex<double>>& x) {
    size_t N = x.size();
    if (N <= 1) return x;
    
    // Verificar que N sea potencia de 2
    if ((N & (N - 1)) != 0) {
        throw runtime_error("FFT requiere que el tamaño sea potencia de 2");
    }

    // Separar muestras pares e impares
    vector<complex<double>> even(N / 2), odd(N / 2);
    for (size_t i = 0; i < N / 2; ++i) {
        even[i] = x[2 * i];
        odd[i]  = x[2 * i + 1];
    }

    // FFT recursiva de pares e impares
    vector<complex<double>> Fe = fft(even);
    vector<complex<double>> Fo = fft(odd);

    // Combinar resultados
    vector<complex<double>> F(N);
    const double PI = acos(-1.0);

    for (size_t k = 0; k < N / 2; ++k) {
        complex<double> wk = polar(1.0, -2.0 * PI * k / static_cast<double>(N));
        complex<double> t  = wk * Fo[k];

        F[k]         = Fe[k] + t;
        F[k + N / 2] = Fe[k] - t;
    }

    return F;
}

vector<complex<double>> fft_real(const vector<double>& x) {
    vector<complex<double>> xc(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        xc[i] = complex<double>(x[i], 0.0);   // parte imaginaria = 0
    }
    return fft(xc);
}

// Calcula la siguiente potencia de 2 mayor o igual a n (para zero-padding)
size_t siguiente_potencia2(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

// Recibe la señal normalizada en el tiempo y devuelve el espectro listo para filtrar
vector<complex<double>> obtenerEspectroParaFiltrado(const vector<double>& senalTiempo) {
    // Copiar y rellenar con ceros hasta la siguiente potencia de 2
    vector<double> senal = senalTiempo;
    size_t N_fft = siguiente_potencia2(senal.size());
    senal.resize(N_fft, 0.0);

    // Aplicar FFT
    vector<complex<double>> espectro = fft_real(senal);
    return espectro;
}


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
int main() {
    try {

        vector<double> senal = cargar_normalizar_wav("mi_audio.wav"); 

        // ========= Persona 1: FFT y entrega de espectro =========
        vector<complex<double>> espectro = obtenerEspectroParaFiltrado(senal);
        
        // Información para depuración
        cout << "Señal original: " << senal.size() << " muestras" << endl;
        cout << "Espectro FFT: " << espectro.size() << " componentes" << endl;

        // Obtener frecuencia de muestreo real del archivo
        drwav wav;
        string ruta = string("audios/") + "mi_audio.wav";
        double fs = 44100.0; // valor por defecto
        if (drwav_init_file(&wav, ruta.c_str(), NULL)) {
            fs = wav.sampleRate;
            drwav_uninit(&wav);
        }
        cout << "Frecuencia de muestreo: " << fs << " Hz" << endl;

        // ========= Persona 3: filtrado de frecuencias cardíacas =========
        cout << "Aplicando filtro de frecuencias cardíacas (0.5 - 3.5 Hz)..." << endl;
        filtrarFrecuencias(espectro, fs);
        cout << "Filtrado completado." << endl;

        // A partir de aquí siguen Persona 4 (IFFT), 5 (BPM), 6 (anomalías), etc.

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
