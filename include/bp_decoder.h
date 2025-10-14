
#include "matrixUtils.h"
#include "dvbs2ldpcShort.h"
#include "structs.h"

// Estructura para el decodificador BP
typedef struct DecoderConfig DecoderConfig;

struct DecoderConfig
{
	int max_iter;
    float alpha;        // Scaling factor (0.7-0.9)
    float beta;         // Offset (0.0-0.5)
};


// Funciones principales
int decode_minsum(const LDPCCode *code, const float *llr_in, 
                  int *bits_out, DecoderConfig config);

