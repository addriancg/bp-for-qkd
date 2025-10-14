#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "getchecknodetable.h"

Matrices getchecknodetable(double rate)
{
	Matrices m;

	if(rate == 1.0 / 2.0){

		int ct1_rows = 5, ct1_cols = 8;
        int ct2_rows = 15, ct2_cols = 3;

		m.ct1.rows = ct1_rows;
        m.ct1.cols = ct1_cols;
        m.ct1.data = malloc(ct1_rows * sizeof(int *));
        for (int i = 0; i < ct1_rows; i++) {
            m.ct1.data[i] = malloc(ct1_cols * sizeof(int));
        }

		static int ct1[5][8] = {
    		{20, 712, 2386, 6354, 4061, 1062, 5045, 5158},
    		{21, 2543, 5748, 4822, 2348, 3089, 6328, 5876},
    		{22, 926, 5701, 269, 3693, 2438, 3190, 3507},
    		{23, 2802, 4520, 3577, 5324, 1091, 4667, 4449},
    		{24, 5140, 2003, 1263, 4742, 6497, 1185, 6202}
		};	

		for (int i = 0; i < ct1_rows; i++) {
            for (int j = 0; j < ct1_cols; j++) {
                m.ct1.data[i][j] = ct1[i][j];
            }
        }

		m.ct2.rows = ct2_rows;
        m.ct2.cols = ct2_cols;
        m.ct2.data = malloc(ct2_rows * sizeof(int *));
        for (int i = 0; i < ct2_rows; i++) {
            m.ct2.data[i] = malloc(ct2_cols * sizeof(int));
        }

		static int ct2[15][3] = {
			{0, 4046, 6934},
    		{1, 2855, 66},
    		{2, 6694, 212},
    		{3, 3439, 1158},
    		{4, 3850, 4422},
    		{5, 5924, 290},
    		{6, 1467, 4049},
    		{7, 7820, 2242},
    		{8, 4606, 3080},
    		{9, 4633, 7877},
    		{10, 3884, 6868},
    		{11, 8935, 4996},
    		{12, 3028, 764},
    		{13, 5988, 1057},
    		{14, 7411, 3450}
		};

		for (int i = 0; i < ct2_rows; i++) {
            for (int j = 0; j < ct2_cols; j++) {
                m.ct2.data[i][j] = ct2[i][j];
            }
        }
	}
	return m;
}

void free_matrices(Matrices *m) {
    // Liberar ct1
    if (m->ct1.data) {
        for (int i = 0; i < m->ct1.rows; i++) {
            free(m->ct1.data[i]);
        }
        free(m->ct1.data);
    }
    
    // Liberar ct2
    if (m->ct2.data) {
        for (int i = 0; i < m->ct2.rows; i++) {
            free(m->ct2.data[i]);
        }
        free(m->ct2.data);
    }
}