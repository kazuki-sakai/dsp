CC = gcc
RM = rm -f

CFLAGS  = -std=c99 -O2 -g \
			-Wall -Werror -Wmissing-prototypes -Wstrict-prototypes \
			-Wpointer-arith -Wcast-qual -Wwrite-strings \
			-fno-common -fno-strict-aliasing -Wcast-align -Wnested-externs \
			-Wno-format-zero-length -Wunused -Wunused-result
#			-DDEBUG \
#			-Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings
LDFLAGS = 
LDLIBS  = -lm

.SUFFIXES: .c.o

TARGETS = lowpass highpass bandpass bandstop highpass2 bandstop2
OBJS    = lowpass.o highpass.o bandpass.o bandstop.o highpass2.o bandpass2.o dsplib.o iir.o

all: $(TARGETS)

lowpass: lowpass.o dsplib.o iir.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

highpass: highpass.o dsplib.o iir.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

highpass2: highpass2.o dsplib.o iir.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

bandpass: bandpass.o dsplib.o iir.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

bandstop: bandstop.o dsplib.o iir.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

bandstop2: bandstop2.o dsplib.o iir.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

test_wave: test_wave.o dsplib.o iir.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

lowpass.o: lowpass.c dsplib.c dsplib.h iir.c iir.h
highpass.o: highpass.c dsplib.c dsplib.h iir.c iir.h
bandpass.o: bandpass.c dsplib.c dsplib.h iir.c iir.h
bandstop.o: bandstop.c dsplib.c dsplib.h iir.c iir.h
bandstop2.o: bandstop2.c dsplib.c dsplib.h iir.c iir.h
highpass2.o: highpass2.c dsplib.c dsplib.h iir.c iir.h
test_wave.o: test_wave.c dsplib.c dsplib.h iir.c iir.h

dsplib.o: dsplib.c dsplib.h iir.c iir.h

iir.o: iir.c iir.h

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	$(RM) $(TARGETS) $(OBJS)

all_clean:
	$(RM) $(TARGETS) *.o test_wave
