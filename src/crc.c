#include <string.h>
#include "crc.h"

const uint16_t CRC16_POLY = 0x1021;
const uint16_t CRC16_INIT = 0xFFFF;

uint16_t crc16(const int *bits, int len) {
    uint16_t crc = CRC16_INIT;
    for (int i = 0; i < len; i++) {
        crc ^= (bits[i] << 15);  // meter el bit en el MSB
        for (int j = 0; j < 8; j++) {
            if (crc & 0x8000)
                crc = (crc << 1) ^ CRC16_POLY;
            else
                crc <<= 1;
        }
    }
    return crc;
}

void crc16_to_bits(uint16_t crc, int *crc_bits) {
    for (int i = 0; i < 16; i++) {
        crc_bits[i] = (crc >> (15 - i)) & 1;  // MSB primero
    }
}

uint16_t bits_to_crc16(const int *crc_bits) {
    uint16_t crc = 0;
    for (int i = 0; i < 16; i++) {
        crc = (crc << 1) | (crc_bits[i] & 1);  // MSB primero
    }
    return crc;
}

int crc_verify_embedded(const int *decoded, int total_len) {
    int payload_len = total_len - 16;
    uint16_t crc_calc = crc16(decoded, payload_len);
    uint16_t crc_recv = bits_to_crc16(decoded + payload_len);
    return crc_calc == crc_recv;
}
