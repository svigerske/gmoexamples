all : addrow loadcntr

addrow : addrow.o loadgms.o gmomcc.o gevmcc.o

loadcntr : loadcntr.o gmomcc.o gevmcc.o

clean:
	rm -f *.o addrow loadcntr

%.c : gams/apifiles/C/api/%.c
	cp $< $@

LDFLAGS = -ldl -Wl,-rpath,\$$ORIGIN
CFLAGS = -Igams/apifiles/C/api -DGAMSDIR=\"gams\" -g