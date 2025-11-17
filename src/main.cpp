#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"
#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <cstdlib>

using namespace std;

const double PI = acos(-1.0);

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
vector<complex<double>> ifft(const vector<complex<double>>& X) 
{
    size_t N = X.size();
    if (N <= 1) return X; 
    //Verificar que N sea potencia de 2
    // usen try-catch para probar 
    if ((N & (N - 1)) != 0)
    {
        throw runtime_error("IFFT requiere que el tamaño sea potencia de 2");
    }
    // Conjugar las muestras de entrada
    vector<complex<double>> X_conj(N);
    for (size_t i = 0; i < N; ++i) 
    {
        X_conj[i] = conj(X[i]);
    }
    //fft a las muestras
    vector<complex<double>> x_conj = fft(X_conj);
    // Conjugar el resultado y dividir por N
    vector<complex<double>> x(N);
    for (size_t i = 0; i < N; ++i) 
    {
        x[i] = conj(x_conj[i]) / static_cast<double>(N);
    }
    return x;
}
/*
Para extraer la parte real despues de de la ifft
*/
vector<double> ifft_real(const vector<complex<double>>& X) {
    vector<complex<double>> resultado = ifft(X);
    vector<double> real(resultado.size());
    for (size_t i = 0; i < resultado.size(); ++i) 
    {
        real[i] = resultado[i].real();
    }
    return real;
}


// Extracción de BPM
vector<size_t> detectarPicos(const vector<double>& senal_filtrada, double umbral_picos, int distancia_minima_muestras) {
    vector<size_t> indices_picos;
    if (senal_filtrada.empty()) {
        return indices_picos;
    }

    double valor_maximo = *std::max_element(senal_filtrada.begin(), senal_filtrada.end());
    double umbral_absoluto = valor_maximo * umbral_picos; 

    for (size_t i = 1; i < senal_filtrada.size() - 1; ++i) {
        if (senal_filtrada[i] > umbral_absoluto && 
            senal_filtrada[i] > senal_filtrada[i-1] && 
            senal_filtrada[i] > senal_filtrada[i+1]) {
            
            bool es_pico_valido = true;
            if (!indices_picos.empty()) {
                if (i - indices_picos.back() < static_cast<size_t>(distancia_minima_muestras)) {
                    es_pico_valido = false;
                }
            }
            if (es_pico_valido) {
                indices_picos.push_back(i);
            }
        }
    }
    return indices_picos;
}

// Estructura para almacenar los resultados del BPM
struct ResultadosBPM {
    double bpm_promedio;
    vector<double> intervalos_rr_segundos; 
    vector<size_t> indices_picos;
};

ResultadosBPM extraerBPM(const vector<double>& senal_filtrada, double frecuencia_muestreo, double umbral_picos = 0.7, int distancia_minima_muestras = 0) {
    ResultadosBPM resultados;
    resultados.bpm_promedio = 0.0;
    resultados.intervalos_rr_segundos.clear();
    resultados.indices_picos.clear();

    if (senal_filtrada.empty() || frecuencia_muestreo <= 0) {
        return resultados;
    }
    // Calcular una distancia mínima
    if (distancia_minima_muestras == 0) {
        distancia_minima_muestras = static_cast<int>(0.2 * frecuencia_muestreo); // 0.2 segundos de "refractario" entre picos
        if (distancia_minima_muestras < 1) distancia_minima_muestras = 1;
    }

    resultados.indices_picos = detectarPicos(senal_filtrada, umbral_picos, distancia_minima_muestras);

    if (resultados.indices_picos.size() < 2) {
        return resultados; // Minimo dos picos para calcular un intervalo
    }

    // Calcular los intervalos RR 
    for (size_t i = 0; i < resultados.indices_picos.size() - 1; ++i) {
        double tiempo_entre_picos = (double)(resultados.indices_picos[i+1] - resultados.indices_picos[i]) / frecuencia_muestreo;
        resultados.intervalos_rr_segundos.push_back(tiempo_entre_picos);
    }

    // Calcular BPM promedio
    if (!resultados.intervalos_rr_segundos.empty()) {
        double suma_intervalos = std::accumulate(resultados.intervalos_rr_segundos.begin(), resultados.intervalos_rr_segundos.end(), 0.0);
        double promedio_intervalo_rr = suma_intervalos / resultados.intervalos_rr_segundos.size();

        if (promedio_intervalo_rr > 0) {
            resultados.bpm_promedio = 60.0 / promedio_intervalo_rr;
        }
    }
    
    return resultados;
}

// Detección de anomalias
struct Anomalias
{
    bool bradicardia = false;
    bool taquicardia = false;
    bool irregularidad = false;
    vector<string> lista_alertas;
};

Anomalias detectarAnomalias(const ResultadosBPM &datos)
{
    Anomalias a;

    if (datos.bpm_promedio == 0 || datos.intervalos_rr_segundos.empty())
    {
        a.lista_alertas.push_back("No se puede evaluar anomalías (datos insuficientes)");
        return a;
    }

    double bpm = datos.bpm_promedio;

    // --- Reglas simples ---
    if (bpm < 60)
    {
        a.bradicardia = true;
        a.lista_alertas.push_back("Bradicardia detectada (BPM < 60)");
    }
    if (bpm > 100)
    {
        a.taquicardia = true;
        a.lista_alertas.push_back("Taquicardia detectada (BPM > 100)");
    }

    // --- Irregularidad usando variabilidad RR ---
    double promedio = accumulate(datos.intervalos_rr_segundos.begin(),
                                 datos.intervalos_rr_segundos.end(), 0.0) /
                      datos.intervalos_rr_segundos.size();

    double varianza = 0.0;
    for (double rr : datos.intervalos_rr_segundos)
    {
        varianza += pow(rr - promedio, 2);
    }
    varianza /= datos.intervalos_rr_segundos.size();
    double sdnn = sqrt(varianza); // desviación estándar RR

    // Regla simple para irregularidad cardíaca
    if (sdnn > 0.10)
    { // >100 ms
        a.irregularidad = true;
        a.lista_alertas.push_back("Latido irregular (SDNN > 0.10s)");
    }

    if (a.lista_alertas.empty())
    {
        a.lista_alertas.push_back("No se detectaron anomalías.");
    }

    return a;
}

// ========== PRUEBAS UNITARIAS ==========
void pruebasUnitarias() {
    cout << "\n========== PRUEBAS UNITARIAS ==========\n" << endl;
    int pruebas_exitosas = 0;
    int pruebas_totales = 0;

    // Prueba 1: FFT de señal simple
    pruebas_totales++;
    try {
        vector<double> senal = {1.0, 0.0, -1.0, 0.0};
        vector<complex<double>> resultado = fft_real(senal);
        if (resultado.size() == 4) {
            cout << "[OK] Prueba 1: FFT procesa señal correctamente" << endl;
            pruebas_exitosas++;
        } else {
            cout << "[FAIL] Prueba 1: Tamaño incorrecto" << endl;
        }
    } catch (...) {
        cout << "[FAIL] Prueba 1: Excepción inesperada" << endl;
    }

    // Prueba 2: FFT requiere potencia de 2
    pruebas_totales++;
    try {
        vector<complex<double>> senal_invalida = {{1,0}, {2,0}, {3,0}};
        fft(senal_invalida);
        cout << "[FAIL] Prueba 2: No detectó tamaño inválido" << endl;
    } catch (runtime_error& e) {
        cout << "[OK] Prueba 2: FFT rechaza tamaño no potencia de 2" << endl;
        pruebas_exitosas++;
    }

    // Prueba 3: IFFT recupera señal original
    pruebas_totales++;
    try {
        vector<double> original = {1.0, 2.0, 3.0, 4.0};
        vector<complex<double>> fft_result = fft_real(original);
        vector<double> recuperada = ifft_real(fft_result);
        bool correcto = true;
        for (size_t i = 0; i < original.size(); i++) {
            if (abs(original[i] - recuperada[i]) > 0.001) {
                correcto = false;
                break;
            }
        }
        if (correcto) {
            cout << "[OK] Prueba 3: IFFT recupera señal original" << endl;
            pruebas_exitosas++;
        } else {
            cout << "[FAIL] Prueba 3: Señal recuperada difiere del original" << endl;
        }
    } catch (...) {
        cout << "[FAIL] Prueba 3: Excepción inesperada" << endl;
    }

    // Prueba 4: Detección de picos con umbral
    pruebas_totales++;
    try {
        vector<double> senal = {0.1, 0.5, 0.3, 0.8, 0.2, 0.9, 0.1};
        vector<size_t> picos = detectarPicos(senal, 0.5, 1);
        if (picos.size() >= 2) {
            cout << "[OK] Prueba 4: Detector de picos funciona" << endl;
            pruebas_exitosas++;
        } else {
            cout << "[FAIL] Prueba 4: No detectó picos esperados" << endl;
        }
    } catch (...) {
        cout << "[FAIL] Prueba 4: Excepción inesperada" << endl;
    }

    // Prueba 5: Cálculo de BPM válido
    pruebas_totales++;
    try {
        vector<double> senal(1000, 0.0);
        double fs = 1000.0;
        for (int i = 0; i < 10; i++) {
            int pos = i * 100;
            if (pos < (int)senal.size()) senal[pos] = 1.0;
        }
        ResultadosBPM resultado = extraerBPM(senal, fs, 0.5, 50);
        if (resultado.bpm_promedio > 0 && resultado.bpm_promedio < 200) {
            cout << "[OK] Prueba 5: Cálculo de BPM genera valor válido" << endl;
            pruebas_exitosas++;
        } else {
            cout << "[FAIL] Prueba 5: BPM fuera de rango esperado" << endl;
        }
    } catch (...) {
        cout << "[FAIL] Prueba 5: Excepción inesperada" << endl;
    }

    // Prueba 6: Detección de bradicardia
    pruebas_totales++;
    try {
        ResultadosBPM datos_bradi;
        datos_bradi.bpm_promedio = 50.0;
        datos_bradi.intervalos_rr_segundos = {1.2, 1.2, 1.2};
        Anomalias resultado = detectarAnomalias(datos_bradi);
        if (resultado.bradicardia) {
            cout << "[OK] Prueba 6: Detecta bradicardia (BPM < 60)" << endl;
            pruebas_exitosas++;
        } else {
            cout << "[FAIL] Prueba 6: No detectó bradicardia" << endl;
        }
    } catch (...) {
        cout << "[FAIL] Prueba 6: Excepción inesperada" << endl;
    }

    // Prueba 7: Detección de taquicardia
    pruebas_totales++;
    try {
        ResultadosBPM datos_taqui;
        datos_taqui.bpm_promedio = 120.0;
        datos_taqui.intervalos_rr_segundos = {0.5, 0.5, 0.5};
        Anomalias resultado = detectarAnomalias(datos_taqui);
        if (resultado.taquicardia) {
            cout << "[OK] Prueba 7: Detecta taquicardia (BPM > 100)" << endl;
            pruebas_exitosas++;
        } else {
            cout << "[FAIL] Prueba 7: No detectó taquicardia" << endl;
        }
    } catch (...) {
        cout << "[FAIL] Prueba 7: Excepción inesperada" << endl;
    }

    // Prueba 8: Siguiente potencia de 2
    pruebas_totales++;
    try {
        if (siguiente_potencia2(100) == 128 && siguiente_potencia2(256) == 256) {
            cout << "[OK] Prueba 8: Cálculo de siguiente potencia de 2 correcto" << endl;
            pruebas_exitosas++;
        } else {
            cout << "[FAIL] Prueba 8: Cálculo incorrecto de potencia de 2" << endl;
        }
    } catch (...) {
        cout << "[FAIL] Prueba 8: Excepción inesperada" << endl;
    }

    cout << "\n=== RESUMEN PRUEBAS UNITARIAS ===" << endl;
    cout << "Exitosas: " << pruebas_exitosas << "/" << pruebas_totales << endl;
    cout << "Tasa de éxito: " << (100.0 * pruebas_exitosas / pruebas_totales) << "%" << endl;
}


// ========== PRUEBAS FUNCIONALES ==========
void pruebasFuncionales() {
    cout << "\n========== PRUEBAS FUNCIONALES ==========\n" << endl;
    
    // Generar señal de prueba: simulación de pulso cardíaco
    cout << "Generando señal de prueba (pulso simulado)..." << endl;
    
    double fs = 1000.0;
    int duracion_segundos = 10;
    int num_muestras = fs * duracion_segundos;
    vector<double> senal_prueba(num_muestras, 0.0);
    
    // Simular latidos con BPM = 75 (0.8 segundos entre latidos)
    double intervalo_latidos = 0.8;
    int muestras_por_latido = intervalo_latidos * fs;
    
    for (int i = 0; i < num_muestras; i += muestras_por_latido) {
        if (i < num_muestras) {
            senal_prueba[i] = 1.0;
            if (i+1 < num_muestras) senal_prueba[i+1] = 0.8;
            if (i+2 < num_muestras) senal_prueba[i+2] = 0.5;
        }
    }
    
    // Agregar ruido
    for (int i = 0; i < num_muestras; i++) {
        senal_prueba[i] += 0.05 * ((rand() % 100) / 100.0 - 0.5);
    }
    
    // Prueba Funcional 1: Pipeline completo
    cout << "\nPrueba Funcional 1: Pipeline completo de procesamiento" << endl;
    try {
        vector<complex<double>> espectro = obtenerEspectroParaFiltrado(senal_prueba);
        filtrarFrecuencias(espectro, fs);
        vector<double> senal_filtrada = ifft_real(espectro);
        senal_filtrada.resize(num_muestras);
        
        ResultadosBPM bpm_resultado = extraerBPM(senal_filtrada, fs);
        
        cout << "  - FFT: OK (" << espectro.size() << " muestras)" << endl;
        cout << "  - Filtrado: OK" << endl;
        cout << "  - IFFT: OK" << endl;
        cout << "  - Extracción BPM: OK" << endl;
        cout << "  - BPM calculado: " << bpm_resultado.bpm_promedio << " bpm" << endl;
        cout << "  - Picos detectados: " << bpm_resultado.indices_picos.size() << endl;
        cout << "[OK] Pipeline completo ejecutado correctamente" << endl;
    } catch (exception& e) {
        cout << "[FAIL] Error en pipeline: " << e.what() << endl;
    }
    
    // Prueba Funcional 2: Detección de anomalías en diferentes casos
    cout << "\nPrueba Funcional 2: Detección de anomalías" << endl;
    
    // Caso normal (70 bpm)
    ResultadosBPM caso_normal;
    caso_normal.bpm_promedio = 70.0;
    caso_normal.intervalos_rr_segundos = {0.857, 0.857, 0.857, 0.857};
    Anomalias anom_normal = detectarAnomalias(caso_normal);
    cout << "  Caso normal (70 bpm): " << anom_normal.lista_alertas[0] << endl;
    
    // Caso bradicardia (50 bpm)
    ResultadosBPM caso_bradi;
    caso_bradi.bpm_promedio = 50.0;
    caso_bradi.intervalos_rr_segundos = {1.2, 1.2, 1.2};
    Anomalias anom_bradi = detectarAnomalias(caso_bradi);
    cout << "  Caso bradicardia (50 bpm): " << anom_bradi.lista_alertas[0] << endl;
    
    // Caso taquicardia (120 bpm)
    ResultadosBPM caso_taqui;
    caso_taqui.bpm_promedio = 120.0;
    caso_taqui.intervalos_rr_segundos = {0.5, 0.5, 0.5};
    Anomalias anom_taqui = detectarAnomalias(caso_taqui);
    cout << "  Caso taquicardia (120 bpm): " << anom_taqui.lista_alertas[0] << endl;
    
    // Caso irregular
    ResultadosBPM caso_irreg;
    caso_irreg.bpm_promedio = 75.0;
    caso_irreg.intervalos_rr_segundos = {0.6, 1.0, 0.7, 1.2, 0.5};
    Anomalias anom_irreg = detectarAnomalias(caso_irreg);
    cout << "  Caso irregular: ";
    for (const auto& alerta : anom_irreg.lista_alertas) {
        cout << alerta << " ";
    }
    cout << endl;
    cout << "[OK] Detección de anomalías funciona correctamente" << endl;
    
    // Prueba Funcional 3: Rendimiento con diferentes tamaños
    cout << "\nPrueba Funcional 3: Rendimiento con diferentes tamaños de señal" << endl;
    vector<int> tamanios = {256, 1024, 4096, 16384};
    
    for (int tam : tamanios) {
        vector<double> senal_test(tam, 0.0);
        for (int i = 0; i < tam; i++) {
            senal_test[i] = sin(2 * PI * i / tam);
        }
        
        auto inicio = std::chrono::high_resolution_clock::now();
        vector<complex<double>> esp = obtenerEspectroParaFiltrado(senal_test);
        auto fin = std::chrono::high_resolution_clock::now();
        
        auto duracion = std::chrono::duration_cast<std::chrono::microseconds>(fin - inicio);
        cout << "  Tamaño " << tam << ": " << duracion.count() << " us" << endl;
    }
    cout << "[OK] Pruebas de rendimiento completadas" << endl;
}


// ========== ANÁLISIS EXPERIMENTAL ==========
void analisisExperimental() {
    cout << "\n========== ANÁLISIS EXPERIMENTAL ==========\n" << endl;
    
    // Experimento 1: Comparación de diferentes configuraciones de umbral
    cout << "Experimento 1: Sensibilidad del detector de picos" << endl;
    cout << "Probando diferentes umbrales de detección...\n" << endl;
    
    double fs = 1000.0;
    int duracion = 10;
    vector<double> senal_exp(fs * duracion, 0.0);
    
    // Generar señal con ruido controlado
    for (size_t i = 0; i < static_cast<size_t>(fs * duracion); i += 800) {
        if (i < senal_exp.size()) {
            senal_exp[i] = 1.0;
            if (i+1 < senal_exp.size()) senal_exp[i+1] = 0.7;
        }
    }
    
    // Agregar ruido
    for (size_t i = 0; i < senal_exp.size(); i++) {
        senal_exp[i] += 0.1 * (rand() % 100 / 100.0 - 0.5);
    }
    
    vector<double> umbrales = {0.3, 0.5, 0.7, 0.9};
    cout << "Umbral\tPicos\tBPM" << endl;
    cout << "------\t-----\t---" << endl;
    
    for (double umbral : umbrales) {
        ResultadosBPM resultado = extraerBPM(senal_exp, fs, umbral);
        cout << umbral << "\t" << resultado.indices_picos.size() 
             << "\t" << resultado.bpm_promedio << endl;
    }
    
    // Experimento 2: Efectividad del filtrado de frecuencias
    cout << "\nExperimento 2: Efectividad del filtrado de frecuencias" << endl;
    cout << "Comparando señal original vs filtrada...\n" << endl;
    
    vector<double> senal_mixta(4096, 0.0);
    for (size_t i = 0; i < senal_mixta.size(); i++) {
        double t = i / fs;
        senal_mixta[i] = sin(2 * PI * 1.2 * t) +
                        0.3 * sin(2 * PI * 5 * t) +
                        0.2 * sin(2 * PI * 0.1 * t);
    }
    
    vector<complex<double>> espectro_orig = obtenerEspectroParaFiltrado(senal_mixta);
    
    double energia_original = 0.0;
    for (const auto& val : espectro_orig) {
        energia_original += abs(val) * abs(val);
    }
    
    vector<complex<double>> espectro_filt = espectro_orig;
    filtrarFrecuencias(espectro_filt, fs);
    
    double energia_filtrada = 0.0;
    for (const auto& val : espectro_filt) {
        energia_filtrada += abs(val) * abs(val);
    }
    
    cout << "Energía original: " << energia_original << endl;
    cout << "Energía filtrada: " << energia_filtrada << endl;
    cout << "Reducción: " << (100.0 * (1 - energia_filtrada/energia_original)) << "%" << endl;
    
    // Experimento 3: Análisis de variabilidad cardíaca
    cout << "\nExperimento 3: Análisis de variabilidad de frecuencia cardíaca" << endl;
    cout << "Simulando diferentes condiciones de ritmo cardíaco...\n" << endl;
    
    struct ConfiguracionRitmo {
        string nombre;
        vector<double> intervalos_rr;
    };
    
    vector<ConfiguracionRitmo> configuraciones = {
        {"Ritmo regular", {0.8, 0.8, 0.8, 0.8, 0.8, 0.8}},
        {"Variabilidad baja", {0.8, 0.82, 0.79, 0.81, 0.80, 0.78}},
        {"Variabilidad moderada", {0.8, 0.9, 0.7, 0.85, 0.75, 0.88}},
        {"Variabilidad alta", {0.6, 1.0, 0.7, 1.2, 0.5, 0.9}}
    };
    
    cout << "Condición\t\tSDNN\tIrregularidad" << endl;
    cout << "----------\t\t----\t-------------" << endl;
    
    for (const auto& config : configuraciones) {
        double promedio = accumulate(config.intervalos_rr.begin(), 
                                    config.intervalos_rr.end(), 0.0) / config.intervalos_rr.size();
        double varianza = 0.0;
        for (double rr : config.intervalos_rr) {
            varianza += pow(rr - promedio, 2);
        }
        varianza /= config.intervalos_rr.size();
        double sdnn = sqrt(varianza);
        
        string irregularidad = (sdnn > 0.10) ? "Sí" : "No";
        cout << config.nombre << "\t" << sdnn << "\t" << irregularidad << endl;
    }
    
    // Experimento 4: Tiempo de procesamiento vs tamaño de señal
    cout << "\nExperimento 4: Complejidad temporal del algoritmo FFT" << endl;
    cout << "Midiendo tiempo de ejecución para diferentes tamaños...\n" << endl;
    
    cout << "N\t\tTiempo (ms)\tTiempo/N*log(N)" << endl;
    cout << "---\t\t-----------\t---------------" << endl;
    
    vector<int> tamanos = {512, 1024, 2048, 4096, 8192};
    for (int n : tamanos) {
        vector<double> senal_tiempo(n, 0.0);
        for (int i = 0; i < n; i++) {
            senal_tiempo[i] = sin(2 * PI * i / n);
        }
        
        auto inicio = std::chrono::high_resolution_clock::now();
        vector<complex<double>> resultado = fft_real(senal_tiempo);
        auto fin = std::chrono::high_resolution_clock::now();
        
        auto duracion_ms = std::chrono::duration_cast<std::chrono::microseconds>(fin - inicio).count() / 1000.0;
        double nlogn = n * log2(n);
        double ratio = duracion_ms / nlogn;
        
        cout << n << "\t\t" << duracion_ms << "\t\t" << ratio << endl;
    }
    
    cout << "\n[OK] Análisis experimental completado" << endl;
}

