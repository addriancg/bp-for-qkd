

#include <stdint.h>

uint16_t crc16(const int *bits, int len);
void crc16_to_bits(uint16_t crc, int *crc_bits);
uint16_t bits_to_crc16(const int *crc_bits);
int crc_verify_embedded(const int *decoded, int total_len);
