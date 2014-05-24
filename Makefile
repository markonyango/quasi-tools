CC=gcc
CFLAGS=-O2 -c -Wall -pedantic -Werror

all: qa count

qa:  qa.o
	$(CC) qa.o -o qa
	@cp qa ~/bin/

count:	count.o fasthash.o output.o
	$(CC) count.o fasthash.o output.o -o count -lm
	@cp count ~/bin/

qa.o: ./Quality_Assessment/qa.c
	$(CC) $(CFLAGS) ./Quality_Assessment/qa.c

count.o: ./Quantification/count.c
	$(CC) $(CFLAGS) ./Quantification/count.c

fasthash.o: ./Quantification/fasthash.c
	$(CC) $(CFLAGS) ./Quantification/fasthash.c

output.o: ./Quantification/output.c
	$(CC) $(CFLAGS) ./Quantification/output.c

clean:
	@echo Removing files...
	@rm -rfv *.o qa count ~/bin/qa ~/bin/count
