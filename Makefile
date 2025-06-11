
CC = gcc
CFLAGS = -Wall -Wextra -O2
LDFLAGS = -lm

SOURCES = main.c Conversions.c externalfunctions.c key_gen.c key_sign.c randombytes.c fips202.c
OBJECTS = $(SOURCES:.c=.o)
TARGET = app

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	-del /Q $(OBJECTS) $(TARGET)

.PHONY: clean
