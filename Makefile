sample: sample.c factor64.o
	gcc -O2 -o sample sample.c factor64.c

factor64.o: factor64.c factor.bin
	gcc -O2 -march=native -c factor64.c

factor.bin: makefactortable
	wget https://people.maths.bris.ac.uk/~maarb/code/factor64_data.xz
	unxz factor64_data.xz
	./makefactortable >> factor64_data
	mv factor64_data factor.bin

makefactortable: makefactortable.c
	gcc -O2 -o makefactortable makefactortable.c
