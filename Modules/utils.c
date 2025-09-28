/**
 * @file utils.c
 * @author David Ramírez Betancourth
 */

#include "utils.h"
#include <hackrf.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>
#include <unistd.h>
#include <string.h>
#include "utils.h"  // tu header con Paths_t, BackendParams_t, etc.

// Callback que se llama por cada bloque de datos recibido
static int rx_callback(hackrf_transfer* transfer) {
    FILE* fp = (FILE*)transfer->rx_ctx;
    if (fp) {
        size_t written = fwrite(transfer->buffer, 1, transfer->valid_length, fp);
        if (written != transfer->valid_length) {
            fprintf(stderr, "Error writing to file\n");
            return -1;
        }
    }
    return 0;
}

int instantaneous_capture(BackendParams_t* params, Paths_t* paths) {
    char outfile_path[2048];
    snprintf(outfile_path, sizeof(outfile_path), "%s/0.cs8", paths->samples_path);

    FILE* fp = fopen(outfile_path, "wb");
    if (!fp) {
        perror("Error opening output file");
        return 1;
    }

    hackrf_device* device = NULL;
    if (hackrf_init() != HACKRF_SUCCESS) {
        fprintf(stderr, "hackrf_init() failed\n");
        fclose(fp);
        return 1;
    }

    if (hackrf_open(&device) != HACKRF_SUCCESS) {
        fprintf(stderr, "hackrf_open() failed\n");
        hackrf_exit();
        fclose(fp);
        return 1;
    }

    // Configuración de sample rate
    if (hackrf_set_sample_rate(device, params->sample_rate) != HACKRF_SUCCESS) {
        fprintf(stderr, "Failed to set sample rate\n");
        hackrf_close(device);
        hackrf_exit();
        fclose(fp);
        return 1;
    }

    // Configuración de frecuencia central
    if (hackrf_set_freq(device, (uint64_t)params->frequency) != HACKRF_SUCCESS) {
        fprintf(stderr, "Failed to set frequency\n");
        hackrf_close(device);
        hackrf_exit();
        fclose(fp);
        return 1;
    }

    // Configuración de ancho de banda del filtro baseband
    if (hackrf_set_baseband_filter_bandwidth(device, params->filter_bw) != HACKRF_SUCCESS) {
        fprintf(stderr, "Failed to set baseband filter bandwidth\n");
        hackrf_close(device);
        hackrf_exit();
        fclose(fp);
        return 1;
    }

    // Ganancia opcional
    hackrf_set_amp_enable(device, 1);

    // Iniciar captura
    hackrf_start_rx(device, rx_callback, fp);

    // Calcular duración equivalente a num_samples
    double capture_time_s = (double)params->num_samples / params->sample_rate;

    printf("Capturing %lu samples (%.3f seconds at %.2f Msps, filter %.2f MHz)...\n",
           params->num_samples,
           capture_time_s,
           params->sample_rate / 1e6,
           params->filter_bw / 1e6);

    usleep((useconds_t)(capture_time_s * 1e6));

    hackrf_stop_rx(device);
    hackrf_close(device);
    hackrf_exit();
    fclose(fp);

    printf("Capture completed. Output saved to: %s\n", outfile_path);
    return 0;
}




signal_iq_t* load_cs8(const char* filename) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        perror("Error: No se pudo abrir el archivo de datos CS8");
        return NULL;
    }

    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    if (file_size % 2 != 0) {
        fprintf(stderr, "Error: Tamaño del archivo inválido. Debe ser múltiplo de 2.\n");
        fclose(file);
        return NULL;
    }

    
    signal_iq_t* signal_data = (signal_iq_t*)malloc(sizeof(signal_iq_t));
    if (!signal_data) {
        perror("Error: No se pudo reservar memoria para la estructura signal_iq_t");
        fclose(file);
        return NULL;
    }

    signal_data->n_signal = file_size / 2;
    int8_t* raw_data = (int8_t*)malloc(file_size);
    signal_data->signal_iq = (complex double*)malloc(signal_data->n_signal * sizeof(complex double));

    if (!raw_data || !signal_data->signal_iq) {
        perror("Error: No se pudo reservar memoria para los datos IQ");
        free(raw_data);
        free(signal_data->signal_iq);
        free(signal_data);
        fclose(file);
        return NULL;
    }

    if (fread(raw_data, 1, (size_t)file_size, file) != (size_t)file_size) {
        perror("Error: Lectura incompleta del archivo");
        free(raw_data);
        free(signal_data->signal_iq);
        free(signal_data);
        fclose(file);
        return NULL;
    }

    for (size_t i = 0; i < signal_data->n_signal; i++) {
        signal_data->signal_iq[i] = raw_data[2 * i] + raw_data[2 * i + 1] * I;
    }

    free(raw_data);
    fclose(file);

    return signal_data;
}

Paths_t get_paths(void) {
    Paths_t paths;
    // Get the current working directory
    char cwd[2048];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        snprintf(paths.samples_path, sizeof(paths.samples_path), "%s/Samples", cwd);
        snprintf(paths.json_path, sizeof(paths.json_path), "%s/JSON", cwd);
    }

    printf("Sample path: %s\n", paths.samples_path);
    printf("JSON path: %s\n", paths.json_path);
    return paths;
}



int hackrf_is_connected() {
    int rc;

    rc = hackrf_init();
    if (rc != HACKRF_SUCCESS) {
        fprintf(stderr, "[HackRF check] hackrf_init failed: %s\n", hackrf_error_name(rc));
        return 0;
    }

    hackrf_device *device = NULL;
    rc = hackrf_open(&device);
    if (rc != HACKRF_SUCCESS) {
        fprintf(stderr, "[HackRF check] hackrf_open failed: %s\n",
                rc == HACKRF_SUCCESS ? "unknown error" : hackrf_error_name(rc));
        hackrf_exit();
        return 0;
    }

    hackrf_close(device);
    hackrf_exit();
    return 1;
}

