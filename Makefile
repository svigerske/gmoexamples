all : addrow

addrow : addrow.o loadgms.o gmomcc.o gevmcc.o

clean:
	rm -f *.o addrow

%.c : gams/apifiles/C/api/%.c
	cp $< $@

LDFLAGS = -ldl -Wl,-rpath,\$$ORIGIN
CFLAGS = -Igams/apifiles/C/api -DGAMSDIR=\"gams\" -g