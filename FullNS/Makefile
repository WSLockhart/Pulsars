#For linux compilers
dgcc= gcc

ray: math.o jmath.o metric.o ns.o io.o image.o integrals.o ray.o time.o
	$(dgcc) math.o jmath.o metric.o ns.o io.o image.o integrals.o ray.o time.o -o ray -lm

ray_spin.o: ray.c definitions.h global.h
	$(dgcc) -c ray.c

time.o: time.c definitions.h global.h
	$(dgcc) -c time.c

integrals.o: integrals.c definitions.h global.h
	$(dgcc) -c integrals.c

image.o: image.c definitions.h global.h
	$(dgcc) -c image.c

io.o: io.c definitions.h global.h
	$(dgcc) -c io.c

ns.o: ns.c definitions.h global.h
	$(dgcc) -c ns.c

metric.o: metric.c definitions.h global.h
	$(dgcc) -c metric.c

math.o: math.c
	$(dgcc) -c math.c

jmath.o: jmath.c
	$(dgcc) -c jmath.c


clean:
	rm -f *.o *.trace *~
