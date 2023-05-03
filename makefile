all: mpi

sr:
	gcc -o sssp_serial src/util.h src/util.c  src/sssp_serial.c 
#	./sssp_serial input/sample_graph_paper_csc_8v_16e_weighted.txt 4 out.txt

mpi:
	mpicc -o sssp_parallel src/util.h src/util.c  src/sssp_parallel.c 
#	mpirun -np 4 ./sssp_parallel input/sample_graph_paper_csc_8v_16e_weighted.txt 4 out.txt

clear:
	@rm -f a.out out.txt sssp_serial

c: clear

z: zip

zip:
	@zip -r muhammed_can_kucukaslan.zip Makefile src/

out16: sr mpi
	mpirun -np 4 ./sssp_parallel input/output16_csc_65536v_247122e_weighted.txt 65474 outmpi4.txt
	mpirun --oversubscribe -np 64 ./sssp_parallel input/output16_csc_65536v_247122e_weighted.txt 65474 outmpi64.txt
	./sssp_serial input/output16_csc_65536v_247122e_weighted.txt 65474 outsr.txt
	diff -s outmpi4.txt outmpi64.txt
	diff -s outmpi4.txt outsr.txt
	diff -s outmpi64.txt outsr.txt

# out of 2^20=1048576 nodes
out20_100: sr mpi
	mpirun --oversubscribe -np 4 ./sssp_parallel input/g20_100_csc_1048576v_87078847e_weighted.txt 65474 outmpi4.txt
	mpirun --oversubscribe -np 64 ./sssp_parallel input/g20_100_csc_1048576v_87078847e_weighted.txt 65474 outmpi64.txt
	./sssp_serial input/g20_100_csc_1048576v_87078847e_weighted.txt 65474 outsr.txt
	diff -s outmpi4.txt outmpi64.txt
	diff -s outmpi4.txt outsr.txt
	diff -s outmpi64.txt outsr.txt

gen:
	cd input
	./generator_omp	10 -e 100 -o g10_100.txt
	./edgelist2csc.py g10_100.txt 1024
