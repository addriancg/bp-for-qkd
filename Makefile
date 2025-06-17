# Nombre del ejecutable
TARGET = prueba_h

# Compilador y opciones
CC = gcc
CFLAGS = -std=c11 -O3 -Iinclude -lm

# Archivos fuente
SRC = \
    prueba_h.c \
    src/dvbs2ldpcShort.c \
    src/matrixUtils.c \
    src/getchecknodetable.c \
    src/crc.c

# Regla principal
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

# Limpieza
clean:
	rm -f $(TARGET)

